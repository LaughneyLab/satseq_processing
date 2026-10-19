## ============================================================
##  plot_IREA_summary.R
##
##  Plot1_radar_polarization.pdf  – overlapping P13/K337 radar grid
##  Plot2_bubble_cytokine.pdf     – P13 | K337 side-by-side bubble matrix
##
##  Usage:  Rscript --vanilla plot_IREA_summary.R [data_dir] [out_dir]
## ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
  library(cluster)
  library(grid)
  library(grDevices)
})

args     <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(args) >= 1) args[1] else "excel"
OUT_DIR  <- if (length(args) >= 2) args[2] else "."
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COL_P13  <- "#FBB040"
COL_K337 <- "#E43B3F"
COL      <- c(P13 = COL_P13, K337 = COL_K337)

CT_LABELS <- c(
  Alveolar_Macrophage     = "Alveolar Macrophage",
  Interstitial_Macrophage = "Interstitial Macrophage",
  Inflammatory_Monocyte   = "Inflammatory Monocyte",
  Neutrophil              = "Neutrophil",
  neutrophil              = "Neutrophil",
  B_cell                  = "B cell",
  NK_cell                 = "NK cell",
  pDC                     = "pDC",
  cDC1                    = "cDC1",
  cDC2                    = "cDC2",
  Naive_CD4_T_cell        = "Naive CD4+ T cell",
  Activated_CD4_T_cell    = "Activated CD4+ T cell",
  Naive_CD8_T_cell        = "Naive CD8+ T cell",
  Exhausted_CD8_T_cell    = "Exhausted CD8+ T cell",
  Treg                    = "Treg",
  T_gd                    = "gd T cell"
)

read_score <- function(prefix, variant, type = c("polar","cytokine")) {
  type <- match.arg(type)
  suf  <- if (type == "polar") paste0("_scores_", variant, ".csv")
          else                  paste0("_scores_", variant, "_cytokine.csv")
  path <- file.path(DATA_DIR, paste0(prefix, suf))
  if (!file.exists(path)) return(NULL)
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(df)[1] <- "label"
  df$variant  <- variant
  df$celltype <- prefix
  df
}

p13_files   <- list.files(DATA_DIR, pattern = "_scores_P13\\.csv$")
CT_PREFIXES <- sub("_scores_P13\\.csv$", "", p13_files)

DISPLAY_ORDER <- c(
  "Alveolar_Macrophage","Interstitial_Macrophage","Inflammatory_Monocyte",
  "neutrophil","Neutrophil","B_cell","NK_cell",
  "pDC","cDC1","cDC2",
  "Naive_CD4_T_cell","Activated_CD4_T_cell",
  "Naive_CD8_T_cell","Exhausted_CD8_T_cell","Treg","T_gd"
)
CT_PREFIXES <- c(intersect(DISPLAY_ORDER, CT_PREFIXES),
                 setdiff(CT_PREFIXES, DISPLAY_ORDER))

cat("Cell types:", length(CT_PREFIXES), "\n")

## ============================================================
##  PLOT 1 — RADAR GRID
## ============================================================

radar_xy <- function(values, angles)
  data.frame(x = values * cos(angles), y = values * sin(angles))

COL_STOP <- "#B4B4B4"

make_radar <- function(df_p13, df_k337, title, df_stop = NULL) {
  states <- sort(Reduce(union, Filter(Negate(is.null),
                   list(df_p13$label, df_k337$label,
                        if (!is.null(df_stop)) df_stop$label))))
  n      <- length(states)
  angles <- (pi / 2) - 2 * pi * (seq_len(n) - 1) / n

  # Positive-only normalisation + floor so polygon never collapses
  FLOOR <- 0.12
  prep  <- function(df) {
    es     <- setNames(as.numeric(df$`Enrichment Score`), df$label)
    padj_v <- setNames(as.numeric(df$padj), df$label)
    es_v   <- es[states];     es_v[is.na(es_v)]       <- 0
    pv     <- padj_v[states]; pv[is.na(pv)]           <- 1
    pos_v  <- pmax(es_v, 0)
    mx     <- max(pos_v, na.rm = TRUE)
    norm_v <- if (mx > 0) pos_v / mx else pos_v
    norm_v <- pmax(norm_v, FLOOR)
    list(norm = norm_v, es = es_v, padj = pv)
  }

  r13   <- prep(df_p13)
  r337  <- prep(df_k337)
  r_stop <- if (!is.null(df_stop)) prep(df_stop) else NULL

  # Outer polygon (no fill, just border)
  bg_xy   <- rbind(radar_xy(rep(1,   n), angles), radar_xy(rep(1,   1), angles)[1,])
  ring_xy <- rbind(radar_xy(rep(0.5, n), angles), radar_xy(rep(0.5, 1), angles)[1,])

  make_poly <- function(vals, var) {
    xy <- radar_xy(vals, angles); xy <- rbind(xy, xy[1,]); xy$variant <- var; xy
  }
  poly_k337 <- make_poly(r337$norm,  "K337")
  poly_p13  <- make_poly(r13$norm,   "P13")
  poly_stop <- if (!is.null(r_stop)) make_poly(r_stop$norm, "STOP") else NULL

  spoke_df <- data.frame(x0=0, y0=0, x1=cos(angles), y1=sin(angles))

  push <- 1.20
  lab_df <- data.frame(
    label = states,
    x     = push * cos(angles),
    y     = push * sin(angles),
    hjust = ifelse(cos(angles) >  0.15, 0, ifelse(cos(angles) < -0.15, 1, 0.5)),
    vjust = ifelse(sin(angles) >  0.15, 0, ifelse(sin(angles) < -0.15, 1, 0.5))
  )

  make_dots <- function(r) {
    data.frame(
      x       = r$norm * cos(angles),
      y       = r$norm * sin(angles),
      has_dot = r$es > 0,
      stringsAsFactors = FALSE
    )
  }
  dots_k337 <- make_dots(r337)
  dots_p13  <- make_dots(r13)
  dots_stop <- if (!is.null(r_stop)) make_dots(r_stop) else NULL

  # 0.5 label: on the spoke closest to east (smallest |angle|)
  idx_east <- which.min(abs(angles))
  ring_lx  <- 0.5 * cos(angles[idx_east]) + 0.04
  ring_ly  <- 0.5 * sin(angles[idx_east])

  p <- ggplot() +
    # Outer polygon — white background, steelblue border only
    geom_polygon(data = bg_xy, aes(x, y),
                 fill = "white", colour = "black", linewidth = 0.4) +
    # Inner 0.5 ring
    geom_polygon(data = ring_xy, aes(x, y),
                 fill = NA, colour = "black", linewidth = 0.25) +
    # Spokes
    geom_segment(data = spoke_df,
                 aes(x = x0, y = y0, xend = x1, yend = y1),
                  colour = "black", linewidth = 0.25)
    # # 0.5 scale label
    # + annotate("text", x = ring_lx, y = ring_ly,
    #          label = "0.5", size = 2.4, colour = "grey55", hjust = 0)

  # Layer order: K337 bottom, STOP middle, P13 top
  p <- p +
    geom_polygon(data = poly_k337, aes(x, y),
                 fill = COL_K337, colour = COL_K337,
                 alpha = 0.18, linewidth = 0.9)

  if (!is.null(poly_stop)) {
    p <- p +
      geom_polygon(data = poly_stop, aes(x, y),
                   fill = COL_STOP, colour = COL_STOP,
                   alpha = 0.18, linewidth = 0.9)
  }

  p <- p +
    geom_polygon(data = poly_p13, aes(x, y),
                 fill = COL_P13, colour = COL_P13,
                 alpha = 0.18, linewidth = 0.9)

  # Dots: K337, STOP, P13 (top)
  p <- p +
    geom_point(data = filter(dots_k337, has_dot), aes(x, y),
               colour = COL_K337, fill = COL_K337,
               shape = 21, size = 2.0, stroke = 0.7)

  if (!is.null(dots_stop)) {
    p <- p +
      geom_point(data = filter(dots_stop, has_dot), aes(x, y),
                 colour = COL_STOP, fill = COL_STOP,
                 shape = 21, size = 2.0, stroke = 0.7)
  }

  p <- p +
    geom_point(data = filter(dots_p13, has_dot), aes(x, y),
               colour = COL_P13, fill = COL_P13,
               shape = 21, size = 2.0, stroke = 0.7) +
    # State labels
    geom_text(data = lab_df,
              aes(x, y, label = label, hjust = hjust, vjust = vjust),
              size = 2.9, colour = "black") +
    annotate("text", x = 0, y = 1.36, label = title,
             size = 3.1, fontface = "bold", colour = "grey20",
             hjust = 0.5, vjust = 0) +
    coord_fixed(xlim = c(-1.52, 1.52), ylim = c(-1.48, 1.56), clip = "off") +
    labs(title = NULL) +
    theme_void(base_size = 9) +
    theme(plot.margin = margin(1, 1, 1, 1))

  p
}

cat("Building radar plots...\n")
radar_plots <- list()
for (ct in CT_PREFIXES) {
  d13   <- read_score(ct, "P13",  "polar")
  d337  <- read_score(ct, "K337", "polar")
  d_stop <- read_score(ct, "STOP", "polar")   # NULL if file absent
  if (is.null(d13) || is.null(d337)) next
  lbl <- CT_LABELS[ct]; if (is.na(lbl)) lbl <- gsub("_"," ",ct)
  radar_plots[[ct]] <- make_radar(d13, d337, lbl, df_stop = d_stop)
}

ncols <- 4
nrows <- ceiling(length(radar_plots) / ncols)

COL_ALL <- c(P13 = COL_P13, K337 = COL_K337, STOP = COL_STOP)
legend_panel <- ggplot(
  data.frame(variant=c("P13","K337","STOP"), x=c(1,2,3), y=c(1,1,1)),
  aes(x, y, colour=variant, fill=variant)
) +
  geom_point(shape=21, size=5, alpha=0.55, stroke=1.4) +
  geom_text(aes(label=variant), nudge_x=0.30, size=4.5, fontface="bold") +
  scale_colour_manual(values=COL_ALL, guide="none") +
  scale_fill_manual(values=COL_ALL, guide="none") +
  xlim(0.5, 4.2) + ylim(0.5, 1.5) + theme_void()

pdf1_path <- file.path(OUT_DIR, "P13_K337_IREA_cell_polarization_radars.pdf")
pdf(pdf1_path, width=12, height=2.1*nrows+0.5, useDingbats=FALSE)
pw_radar <- wrap_plots(c(radar_plots, list(legend_panel)), ncol=ncols) &
  theme(plot.margin = margin(1, 1, 1, 1))
print(pw_radar)
dev.off()
cat("Saved:", pdf1_path, "\n")

## ============================================================
##  PLOT 2 — PAIRED CYTOKINE BUBBLE (one page, P13+K337 paired per cell type)
## ============================================================
cat("Building cytokine bubble plots...\n")

load_cytokine_all <- function(variant) {
  dfs <- lapply(CT_PREFIXES, function(ct) {
    d <- read_score(ct, variant, "cytokine")
    if (is.null(d)) return(NULL)
    lbl <- CT_LABELS[ct]; if (is.na(lbl)) lbl <- gsub("_"," ",ct)
    d$celltype_label <- lbl; d
  })
  bind_rows(Filter(Negate(is.null), dfs))
}

cyto_p13  <- load_cytokine_all("P13")
cyto_k337 <- load_cytokine_all("K337")

# ── Keep only cell types present in BOTH conditions ──────────────────────────
ct_both <- intersect(unique(cyto_p13$celltype_label),
                     unique(cyto_k337$celltype_label))
cat("Cell types with both P13 and K337:", length(ct_both), "\n")
cyto_p13  <- filter(cyto_p13,  celltype_label %in% ct_both)
cyto_k337 <- filter(cyto_k337, celltype_label %in% ct_both)

# ── Filter: keep cytokines with min padj < 0.01 (-log10 >= 2) in ≥1 condition ──
keep_cyto <- bind_rows(cyto_p13, cyto_k337) %>%
  group_by(label) %>%
  summarise(min_padj = min(padj, na.rm=TRUE), .groups="drop") %>%
  filter(min_padj < 0.001) %>%
  pull(label)
cat("Cytokines kept (min padj < 0.001):", length(keep_cyto), "\n")
cyto_p13  <- filter(cyto_p13,  label %in% keep_cyto)
cyto_k337 <- filter(cyto_k337, label %in% keep_cyto)

# Replace Greek letters with ASCII equivalents so PDF fonts render cleanly.
# Build a lookup from raw bytes in a locale-safe way.
# Replace Greek UTF-8 multi-byte sequences with ASCII letters.
# The file data is UTF-8 but the system locale may not be; match raw bytes.
clean_cyto_name <- function(x) {
  # Each replacement: raw UTF-8 two-byte sequence -> ASCII letter
  subs <- list(
    c("\xce\xb1", "a"),   # α  U+03B1
    c("\xce\xb2", "b"),   # β  U+03B2
    c("\xce\xb3", "g"),   # γ  U+03B3
    c("\xce\xb5", "e"),   # ε  U+03B5
    c("\xce\xba", "k"),   # κ  U+03BA
    c("\xce\xbb", "l"),   # λ  U+03BB
    c("\xce\xb4", "d"),   # δ  U+03B4
    c("\xce\xb1", "a"),   # sometimes stored as ce b1
    c("\xce\xa9", "O")    # Ω  just in case
  )
  xb <- chartr("", "", x)   # no-op copy
  for (s in subs) {
    xb <- gsub(s[1], s[2], xb, fixed=TRUE, useBytes=TRUE)
  }
  xb
}
cyto_p13$label  <- clean_cyto_name(cyto_p13$label)
cyto_k337$label <- clean_cyto_name(cyto_k337$label)
keep_cyto       <- clean_cyto_name(keep_cyto)

# ── Build combined matrix for clustering ──────────────────────────────────────
build_mat <- function(df) {
  df2 <- df %>%
    select(Cytokine=label, CellType=celltype_label, ES=`Enrichment Score`) %>%
    pivot_wider(names_from=Cytokine, values_from=ES)
  mat <- as.matrix(df2[,-1]); rownames(mat) <- df2$CellType; mat
}
mat_p13  <- build_mat(cyto_p13)
mat_k337 <- build_mat(cyto_k337)

expand_mat <- function(mat, rows, cols) {
  out <- matrix(NA_real_, nrow=length(rows), ncol=length(cols),
                dimnames=list(rows,cols))
  out[rownames(mat), colnames(mat)] <- mat; out
}
mat_p13_full  <- expand_mat(mat_p13,  ct_both, keep_cyto)
mat_k337_full <- expand_mat(mat_k337, ct_both, keep_cyto)
mat_comb <- (abs(replace(mat_p13_full,  is.na(mat_p13_full),  0)) +
             abs(replace(mat_k337_full, is.na(mat_k337_full), 0))) / 2

cor_dist <- function(m) {
  cc <- cor(m, use="pairwise.complete.obs", method="spearman"); cc[is.na(cc)] <- 0
  as.dist(1 - cc)
}
col_order <- {hc <- hclust(cor_dist(mat_comb), method="ward.D2"); hc$labels[hc$order]}
row_order  <- {hc <- hclust(dist(mat_comb, method="euclidean"), method="ward.D2");
               rev(hc$labels[hc$order])}   # top = most activated

# ── Prepare long data for each variant with shared column + row order ──────────
prep_df <- function(cyto_df) {
  cyto_df %>%
    rename(Cytokine = label, CellType = celltype_label, ES = `Enrichment Score`) %>%
    filter(!is.na(ES)) %>%
    mutate(
      neg_log10_padj = pmin(-log10(pmax(padj, 1e-10)), 10),
      Cytokine = factor(Cytokine, levels = col_order),
      CellType = factor(CellType, levels = row_order)
    ) %>%
    filter(!is.na(Cytokine), !is.na(CellType))
}
df_p13  <- prep_df(cyto_p13)
df_k337 <- prep_df(cyto_k337)

# ── Shared scale limits across both panels ─────────────────────────────────────
es_lim <- 0.7*max(abs(c(df_p13$ES, df_k337$ES)), na.rm = TRUE)
max_sz <- max(c(df_p13$neg_log10_padj, df_k337$neg_log10_padj), na.rm = TRUE)

n_cyto  <- length(col_order)
n_cells <- length(row_order)

# ── Panel builder ──────────────────────────────────────────────────────────────
# show_y     : show y-axis cell type labels (left panel only)
# show_leg   : show legend (right panel only)
# left_margin: extra left margin in pt for cell type label annotation (left panel)
make_panel <- function(df, title_label, title_col,
                       show_y = TRUE, show_leg = FALSE,
                       left_margin = 8) {

  p <- ggplot(df, aes(x = Cytokine, y = CellType)) +
    geom_point(aes(size = neg_log10_padj, fill = ES, colour = ES),
               shape = 21, stroke = 0.15, alpha = 0.92) +
    scale_size_continuous(
      name   = expression(-log[10](padj)),
      range  = c(0.5, 7),
      limits = c(0, max_sz),
      breaks = c(1, 2, 3, 5),
      guide  = if (show_leg)
        guide_legend(override.aes = list(fill = "grey50", colour = "grey50", alpha = 0.8))
      else "none"
    ) +
    scale_fill_gradientn(
      name    = "Enrichment\nScore",
      colours = c("#2166ac","#4393c3","#92c5de","#f7f7f7",
                  "#f4a582","#d6604d","#b2182b"),
      limits  = c(-es_lim, es_lim), oob = squish,
      guide   = if (show_leg) "colorbar" else "none"
    ) +
    scale_colour_gradientn(
      colours = c("#2166ac","#4393c3","#92c5de","#f7f7f7",
                  "#f4a582","#d6604d","#b2182b"),
      limits  = c(-es_lim, es_lim), oob = squish, guide = "none"
    ) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_y_discrete(expand = expansion(add = 0.45)) +
    # Coloured title above each panel acts as the condition label
    labs(x = NULL, y = NULL,
         title = title_label) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0.5,
                                   colour = title_col, margin = margin(b = 4)),
      axis.text.x   = element_text(angle = 60, hjust = 1, vjust = 1,
                                   size = 8.5, colour = "grey20"),
      axis.text.y   = if (show_y)
                        element_text(size = 10, colour = "grey20")
                      else element_blank(),
      axis.ticks    = element_blank(),
      axis.line.x   = element_line(colour = "grey40", linewidth = 0.4),
      axis.line.y   = element_line(colour = "grey40", linewidth = 0.4),
      panel.grid    = element_blank(),
      legend.position  = if (show_leg) "right" else "none",
      legend.key.size  = unit(0.4, "cm"),
      legend.text      = element_text(size = 9),
      legend.title     = element_text(size = 9.5, face = "bold"),
      plot.margin      = margin(4, 4, 4, left_margin)
    )

  p
}

p_left  <- make_panel(df_p13,  "P13 vs REF",  COL_P13,
                      show_y = TRUE,  show_leg = FALSE, left_margin = 8)
p_right <- make_panel(df_k337, "K337 vs REF", COL_K337,
                      show_y = FALSE, show_leg = TRUE,  left_margin = 4)

# ── Thin vertical divider between the two panels ──────────────────────────────
# # Drawn as a zero-width spacer panel containing a vertical segment
# divider <- ggplot() +
#   annotate("segment", x = 0.5, xend = 0.5, y = 0, yend = 1,
#            colour = "grey60", linewidth = 0.5) +
#   theme_void() +
#   theme(plot.margin = margin(0, 0, 0, 0))

# ── Assemble: left panel | thin divider | right panel ─────────────────────────
combined <- p_left + p_right +
  plot_layout(widths = c(1, 1), guides = "keep")

# ── Save ──────────────────────────────────────────────────────────────────────
pw_pdf <- n_cyto * 0.38 + 7.0
ph_pdf <- n_cells * 0.2 + 2.5

pdf2_path <- file.path(OUT_DIR, "P13_K337_IREA_cytokine_bubbleplot.pdf")
pdf(pdf2_path,
    width  = max(15, pw_pdf),
    height = max(5, ph_pdf),
    useDingbats = FALSE)
print(combined)
dev.off()
cat("Saved:", pdf2_path, "\n")
cat("Done.\n")
