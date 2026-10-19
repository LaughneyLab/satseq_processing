## ============================================================
## iris_4panel.R
##
## Iris (circular bar) plots arranged in a square SVG grid.
## Cell types AND cytokines are selected entirely from raw IREA
## output CSVs using dual padj + |ES| cutoffs — nothing is
## hard-coded; N_PANELS is removed.
##
## Selection logic
##   CYTOKINES
##     Keep a cytokine if, across all cell types and variants,
##     it has at least one observation where
##       padj < CYTO_PADJ_CUTOFF  AND  |ES| > CYTO_ES_CUTOFF
##
##   CELL TYPES
##     Keep a cell type if, restricted to the surviving cytokines,
##     it has at least one observation where
##       padj < CT_PADJ_CUTOFF  AND  |ES| > CT_ES_CUTOFF
##     (must have data for ALL requested variants)
##
##   GRID LAYOUT
##     Number of panels = however many cell types survive.
##     Grid columns = NCOLS (default 2).
##     SVG dimensions are set to be square.
##
## All thresholds are tunable at the top of the script.
## ============================================================
rm(list = ls())
# setwd auto-detect removed — run from script directory

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  # svglite not available — using base svg()
})

# ── User-tunable parameters ───────────────────────────────────────────────────
DATA_DIR  <- "excel"

VARIANTS  <- c("P13", "K337")          # variants overlaid in each panel

# Cytokine selection: both conditions must be met in ≥1 row
CYTO_PADJ_CUTOFF <- 0.01              # min padj threshold for cytokines
CYTO_ES_CUTOFF   <- 0.02             # min |ES|  threshold for cytokines

# Cell type selection: both conditions must be met in ≥1 row
# (evaluated on the already-filtered cytokine set)
CT_PADJ_CUTOFF   <- 0.01             # min padj threshold for cell types
CT_ES_CUTOFF     <- 0.045            # min |ES|  threshold for cell types

NCOLS    <- 2                          # grid columns; rows computed automatically
SVG_SIZE <- 7                         # SVG side in inches (square)

COL_VARIANTS <- c(P13 = "#FBB040", K337 = "#E43B3F")

# ── Helper: replace Greek UTF-8 bytes with ASCII ──────────────────────────────
clean_name <- function(x) {
  subs <- list(
    c("\xce\xb1", "a"), c("\xce\xb2", "b"), c("\xce\xb3", "g"),
    c("\xce\xb5", "e"), c("\xce\xba", "k"), c("\xce\xbb", "l"),
    c("\xce\xb4", "d")
  )
  for (s in subs) x <- gsub(s[1], s[2], x, fixed = TRUE, useBytes = TRUE)
  x
}

# ── Load all cytokine CSV files ───────────────────────────────────────────────
# Pattern: <CellType>_scores_<Variant>_cytokine.csv
cyto_files <- list.files(DATA_DIR,
                         pattern = "_scores_[^_]+_cytokine\\.csv$",
                         full.names = TRUE)

# Parse cell type and variant from filename
parse_meta <- function(path) {
  bn <- sub("_cytokine\\.csv$", "", basename(path))
  # strip known variant suffixes from the end
  for (v in VARIANTS) {
    suf <- paste0("_scores_", v)
    if (endsWith(bn, suf)) {
      return(data.frame(
        path     = path,
        celltype = sub(paste0(suf, "$"), "", bn),
        variant  = v,
        stringsAsFactors = FALSE
      ))
    }
  }
  NULL
}
meta <- bind_rows(Filter(Negate(is.null), lapply(cyto_files, parse_meta)))
# Keep only variants we want
meta <- filter(meta, variant %in% VARIANTS)

cat(sprintf("Found %d CSV files across %d cell types and %d variants.\n",
            nrow(meta),
            length(unique(meta$celltype)),
            length(unique(meta$variant))))

# Read each file
read_cyto <- function(path, ct, vt) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(df)[1] <- "Cytokine"
  df$Cytokine  <- clean_name(df$Cytokine)
  df$celltype  <- ct
  df$variant   <- vt
  df
}
all_raw <- bind_rows(mapply(read_cyto,
                             meta$path, meta$celltype, meta$variant,
                             SIMPLIFY = FALSE))

# ── Keep cell types present in ALL requested variants ─────────────────────────
ct_per_variant <- all_raw %>%
  group_by(celltype, variant) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = variant, values_from = n)

ct_complete <- ct_per_variant %>%
  filter(if_all(all_of(VARIANTS), ~ !is.na(.x))) %>%
  pull(celltype)

cat(sprintf("Cell types with data for all variants (%s): %d\n",
            paste(VARIANTS, collapse = " & "), length(ct_complete)))

all_data <- filter(all_raw, celltype %in% ct_complete)

# ── Select cytokines: padj < CYTO_PADJ_CUTOFF AND |ES| > CYTO_ES_CUTOFF ──────
# A cytokine passes if at least one (celltype × variant) row satisfies BOTH.
cyto_pass <- all_data %>%
  rename(ES = `Enrichment Score`) %>%
  group_by(Cytokine) %>%
  summarise(
    min_padj = min(padj,    na.rm = TRUE),
    max_absES = max(abs(ES), na.rm = TRUE),
    # flag: at least one row with padj AND |ES| both passing
    any_pass = any(padj < CYTO_PADJ_CUTOFF & abs(ES) > CYTO_ES_CUTOFF,
                   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(any_pass) %>%
  arrange(min_padj, desc(max_absES))

keep_cyto <- cyto_pass$Cytokine
cat(sprintf(
  "Cytokines with padj < %.4f AND |ES| > %.4f in ≥1 row: %d\n",
  CYTO_PADJ_CUTOFF, CYTO_ES_CUTOFF, length(keep_cyto)
))
if (length(keep_cyto) == 0)
  stop("No cytokines pass. Lower CYTO_PADJ_CUTOFF or CYTO_ES_CUTOFF.")

# ── Select cell types: padj < CT_PADJ_CUTOFF AND |ES| > CT_ES_CUTOFF ─────────
# Evaluated only on the kept cytokines.
ct_pass <- all_data %>%
  filter(Cytokine %in% keep_cyto) %>%
  rename(ES = `Enrichment Score`) %>%
  group_by(celltype) %>%
  summarise(
    n_sig     = sum(padj < CT_PADJ_CUTOFF & abs(ES) > CT_ES_CUTOFF,
                    na.rm = TRUE),
    min_padj  = min(padj,    na.rm = TRUE),
    max_absES = max(abs(ES), na.rm = TRUE),
    any_pass  = any(padj < CT_PADJ_CUTOFF & abs(ES) > CT_ES_CUTOFF,
                    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(any_pass) %>%
  arrange(desc(n_sig), min_padj, desc(max_absES))

cat("\nCell types passing padj < ", CT_PADJ_CUTOFF,
    " AND |ES| > ", CT_ES_CUTOFF, ":\n", sep = "")
print(ct_pass[, c("celltype","n_sig","min_padj","max_absES")], n = Inf)

selected_ct <- ct_pass$celltype
cat(sprintf("\n%d cell type(s) selected: %s\n",
            length(selected_ct), paste(selected_ct, collapse = ", ")))
if (length(selected_ct) == 0)
  stop("No cell types pass. Lower CT_PADJ_CUTOFF or CT_ES_CUTOFF.")

# Friendly display labels: underscores → spaces, keep "+" and "+" notation
ct_label <- function(ct) gsub("_", " ", ct)

# ── Filter data to selected cell types and cytokines ─────────────────────────
plot_data <- filter(all_data,
                    celltype %in% selected_ct,
                    Cytokine %in% keep_cyto)

# ── Cytokine order: Ward-D2 clustering on mean |ES| per (cytokine × celltype) ─
mat_for_clust <- plot_data %>%
  rename(ES = `Enrichment Score`) %>%
  group_by(Cytokine, celltype) %>%
  summarise(ES = mean(abs(ES)), .groups = "drop") %>%
  pivot_wider(names_from = celltype, values_from = ES, values_fill = 0) %>%
  tibble::column_to_rownames("Cytokine")

cyto_order <- rownames(mat_for_clust)[
  hclust(dist(mat_for_clust), method = "ward.D2")$order
]
n_cyto <- length(cyto_order)
cat(sprintf("Cytokines in plot: %d\n", n_cyto))

# ── Global max |ES| for consistent radial scaling ────────────────────────────
global_max <- 0.07#max(abs(plot_data$`Enrichment Score`), na.rm = TRUE)

# ── Polar geometry ────────────────────────────────────────────────────────────
BAR_GAP  <- 2.0    # degrees gap between cytokine pairs
BAR_PAD  <- 0.07   # fraction of bar arc kept as intra-pair padding

total_deg <- 360 - n_cyto * BAR_GAP
pair_w    <- total_deg / n_cyto
bar_w     <- pair_w / length(VARIANTS)  # one slot per variant

# Centre angle of each variant bar within a cytokine pair
# Variants are spread evenly across pair_w
var_offsets <- seq(bar_w/2, pair_w - bar_w/2, length.out = length(VARIANTS)) -
               pair_w/2   # relative to pair centre; pair centre is at pair_mid

pair_mid <- BAR_GAP/2 +
            (seq_len(n_cyto) - 1) * (pair_w + BAR_GAP) +
            pair_w/2

# For each variant, absolute bar centre angles across all cytokines
var_mids <- setNames(lapply(seq_along(VARIANTS), function(vi)
  pair_mid + var_offsets[vi]
), VARIANTS)

R_BASE <- 3.5   # zero-line radius
R_OUT  <- 3.5   # max outward extent (positive ES)
R_IN   <- 3.5/2   # max inward extent  (negative ES)

polar_rect <- function(theta_c, half_w, r_bot, r_top, n = 14) {
  t0     <- theta_c - half_w * (1 - BAR_PAD)
  t1     <- theta_c + half_w * (1 - BAR_PAD)
  th     <- c(seq(t0, t1, length.out = n), seq(t1, t0, length.out = n))
  rs     <- c(rep(r_bot, n), rep(r_top, n))
  th_rad <- (90 - th) * pi / 180
  data.frame(x = rs * cos(th_rad), y = rs * sin(th_rad))
}

# ── Build one iris panel ──────────────────────────────────────────────────────
make_iris <- function(ct) {
  lbl  <- ct_label(ct)
  data_ct <- filter(plot_data, celltype == ct)

  # Build bar polygons for every cytokine × variant
  polys <- bind_rows(lapply(seq_along(cyto_order), function(ci) {
    cyto <- cyto_order[ci]
    bind_rows(lapply(seq_along(VARIANTS), function(vi) {
      vt   <- VARIANTS[vi]
      col  <- COL_VARIANTS[vt]
      row  <- filter(data_ct, Cytokine == cyto, variant == vt)
      es   <- if (nrow(row) > 0) row$`Enrichment Score`[1] else 0
      h    <- abs(es) / global_max * R_OUT
      if (h < 1e-5) return(NULL)
      poly <- if (es >= 0)
        polar_rect(var_mids[[vt]][ci], bar_w / 2, R_BASE,     R_BASE + h)
      else
        polar_rect(var_mids[[vt]][ci], bar_w / 2, R_BASE - h, R_BASE)
      poly$id      <- paste(cyto, vt)
      poly$fill    <- unname(col)
      poly$variant <- vt
      poly
    }))
  }))

  # Reference rings
  ring_th <- seq(0, 2 * pi, length.out = 400)
  mk_ring <- function(r, tp)
    data.frame(x = r * cos(ring_th), y = r * sin(ring_th), type = tp)
  ring_df <- bind_rows(
    mk_ring(R_BASE,              "base"),
    mk_ring(R_BASE + R_OUT,      "outer"),
    mk_ring(R_BASE - R_IN,       "inner"),
    mk_ring(R_BASE + R_OUT * .5, "mid_out"),
    mk_ring(R_BASE - R_IN  * .5, "mid_in")
  )

  # Spoke lines
  spoke_th <- (90 - pair_mid) * pi / 180
  spoke_df <- data.frame(
    x0 = (R_BASE - R_IN)  * cos(spoke_th),
    y0 = (R_BASE - R_IN)  * sin(spoke_th),
    x1 = (R_BASE + R_OUT) * cos(spoke_th),
    y1 = (R_BASE + R_OUT) * sin(spoke_th)
  )

  # Cytokine labels — always radially outward
  lab_r  <- R_BASE + R_OUT + 0.42
  lab_th <- (90 - pair_mid) * pi / 180
  lab_df <- data.frame(
    cyto  = cyto_order,
    x     = lab_r * cos(lab_th),
    y     = lab_r * sin(lab_th),
    angle = ifelse(pair_mid <= 180, -(pair_mid - 90), -(pair_mid - 270)),
    hjust = ifelse(pair_mid <= 180, 0, 1)
  )

  # Viewport: outer label radius + margin
  vp <- lab_r + 1.8
  vp <- max(vp, R_BASE + R_OUT + 4.0)

  ggplot() +
    # faint dashed mid rings
    geom_path(data = filter(ring_df, type %in% c("mid_out","mid_in")),
              aes(x, y, group = type),
              colour = "grey88", linewidth = 0.22, linetype = "dashed") +
    # outer/inner boundary rings
    geom_path(data = filter(ring_df, type %in% c("outer","inner")),
              aes(x, y, group = type),
              colour = "grey72", linewidth = 0.42) +
    # spoke guides
    geom_segment(data = spoke_df,
                 aes(x = x0, y = y0, xend = x1, yend = y1),
                 colour = "grey92", linewidth = 0.20) +
    # bars
    geom_polygon(data = polys,
                 aes(x, y, group = id, fill = fill),
                 colour = NA) +
    scale_fill_identity() +
    # zero ring on top of bars
    geom_path(data = filter(ring_df, type == "base"),
              aes(x, y, group = type),
              colour = "grey35", linewidth = 0.75) +
    # hub dot
    annotate("point", x = 0, y = 0, size = 3, colour = "grey55", shape = 16) +
    # cytokine labels
    geom_text(data = lab_df,
              aes(x, y, label = cyto, angle = angle, hjust = hjust),
              size = 3.8, colour = "grey18", vjust = 0.5) +
    coord_fixed(xlim = c(-vp, vp), ylim = c(-vp, vp), clip = "off") +
    labs(title = lbl) +
    theme_void() +
    theme(
      plot.title  = element_text(size = 14, face = "bold", colour = "grey10",
                                 hjust = 0.5, vjust = 1,
                                 margin = margin(b = 2)),
      plot.margin = margin(18, 8, 8, 8)
    )
}

# ── Build patchwork grid (NCOLS wide, rows computed from n panels) ─────────────
N_PANELS   <- length(selected_ct)
nrows_grid <- ceiling(N_PANELS / NCOLS)
cat(sprintf("Building %d panel(s) in a %d × %d grid...\n",
            N_PANELS, nrows_grid, NCOLS))

panels <- lapply(selected_ct, make_iris)

rows <- lapply(seq_len(nrows_grid), function(r) {
  idx <- ((r - 1) * NCOLS + 1) : min(r * NCOLS, N_PANELS)
  wrap_plots(panels[idx], ncol = length(idx))
})
grid_out <- Reduce(`/`, rows) &
  theme(plot.margin = margin(1, 1, 1, 1))

# ── Standalone legend panel centred over the full grid ────────────────────────
# A small ggplot with coloured dots + labels, no axes, transparent background.
n_var  <- length(VARIANTS)
leg_df <- data.frame(
  label   = VARIANTS,
  col     = unname(COL_VARIANTS[VARIANTS]),
  y       = seq(n_var - 1, 0, by = -1),   # top variant at highest y
  stringsAsFactors = FALSE
)

legend_panel <- ggplot(leg_df, aes(x = 0, y = y)) +
  geom_point(aes(colour = I(col), fill = I(col)),
             shape = 21, size = 5, stroke = 0) +
  geom_text(aes(x = 0.3, label = label, colour = I(col)),
            hjust = 0, vjust = 0.5, size = 5, fontface = "bold") +
  scale_x_continuous(limits = c(-0.2, 1.5)) +
  scale_y_continuous(limits = c(-0.5, n_var - 0.5)) +
  theme_void() +
  theme(plot.background = element_rect(fill = NA, colour = NA))

# inset_element coords are fractions of the full assembled grid [0,1].
# Each panel has a title margin of 18 pt at the top and 8 pt elsewhere.
# With SVG_SIZE inches and nrows_grid rows, estimate the title fraction
# so we can shift the legend down to the true visual centre of the irises.
# title_pt / (SVG_SIZE * 72 / nrows_grid) gives the fraction one title bar
# takes per row; the legend needs to shift down by half that per row.
# Legend will be drawn at the absolute centre of the SVG device using grid,
# after the patchwork is printed. No inset needed.

# ── Save as square SVG ────────────────────────────────────────────────────────
out_path <- "IREA_cytokine_cell_type_iris.svg"
svg(out_path, width = SVG_SIZE, height = SVG_SIZE)
print(grid_out)

# Draw legend at the absolute centre of the device using grid primitives.
# grid npc units: (0,0) = bottom-left, (1,1) = top-right of the device.
# x = 0.5, y = 0.5 is always the exact pixel centre regardless of margins.
library(grid)
n_var   <- length(VARIANTS)
leg_dy  <- 0.04          # vertical spacing between items (npc)
leg_y0  <- 0.5 + (n_var - 1) / 2 * leg_dy   # top item y, centred at 0.5
dot_r   <- 0.006          # dot radius (npc)
txt_off <- 0.012          # text x offset from dot centre (npc)
fsize   <- SVG_SIZE * 2.2 # font size in pt, scales with SVG size

for (vi in seq_along(VARIANTS)) {
  vt  <- VARIANTS[vi]
  col <- unname(COL_VARIANTS[vt])
  y_i <- leg_y0 - (vi - 1) * leg_dy

  grid.circle(x = 0.5, y = y_i, r = dot_r,
              gp = gpar(fill = col, col = col))
  grid.text(label = vt,
            x = 0.5 + txt_off, y = y_i,
            just = c("left", "centre"),
            gp = gpar(col = col, fontface = "bold", fontsize = fsize))
}

dev.off()

# Print limits of the iris
ES_outer <-  global_max                          # +ve ring
ES_inner <- -global_max * (R_IN / R_OUT)         # -ve ring
cat(sprintf("Outer ring ES: +%.4f | Inner ring ES: %.4f\n", ES_outer, ES_inner))
