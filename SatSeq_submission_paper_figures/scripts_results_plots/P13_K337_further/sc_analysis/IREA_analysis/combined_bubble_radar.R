suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)
})
DATA_DIR <- "excel"
COL_P13  <- "#FBB040"; COL_K337 <- "#E43B3F"

clean_name <- function(x) {
  subs <- list(c("\xce\xb1","a"),c("\xce\xb2","b"),c("\xce\xb3","g"),
               c("\xce\xb5","e"),c("\xce\xba","k"),c("\xce\xbb","l"),c("\xce\xb4","d"))
  for (s in subs) x <- gsub(s[1],s[2],x,fixed=TRUE,useBytes=TRUE); x
}

FOCUS <- c("Inflammatory_Monocyte","T_gd","Naive_CD8_T_cell","cDC2")
CT_LABELS <- c(
  Inflammatory_Monocyte="Inflammatory Monocyte", T_gd="gd T cell",
  Naive_CD8_T_cell="Naive CD8+ T cell", cDC2="cDC2"
)

POLAR_INFO <- list(
  Inflammatory_Monocyte = list(
    "Mono-a"=list(cy="IFNa1, IFNb",           mk="Ifit3, Rsad2, Ifit2\nIrf7, Isg15"),
    "Mono-b"=list(cy="IFNg, IL-2, IL-12\nIL-15, IL-18", mk="Cxcl10, Socs1, Gbp5\nFam26f, Ifi47"),
    "Mono-c"=list(cy="IL-1a, IL-1b, IL-36a",  mk="Sdc4, Srgn, Dab2\nC5ar1, Adam8"),
    "Mono-d"=list(cy="TNF",                   mk="Mmp14, Gpr84\nBcl2a1b, Nr1h3"),
    "Mono-e"=list(cy="GM-CSF, IL-3, IL-33",   mk="Fabp5, Naaa, Eif5a\nCish, Dok2")
  ),
  T_gd = list(
    "Tgd-a" =list(cy="IFNa1, IFNb\nIFNe, IFNk", mk="Isg15, Ifit1, Bst2\nIfi203, Rnf213"),
    "Tgd-b" =list(cy="IFNg, IL-2, IL-12\nIL-15, IL-18", mk="Igtp, Stat1, Gbp4\nIfi47, Gbp2"),
    "Tgd-c" =list(cy="IL-1a, IL-1b",          mk="Npm1, Ncl, Eif5a\nGzmc, Odc1"),
    "Tgd-d" =list(cy="TNF, TL1A, IL-17E",     mk="Odc1, Adam8, Ltb\nCd82, Cyba"),
    "Tgd-e" =list(cy="IL-2, IL-7, IL-15\nIL-18, IL-36a",mk="Eif5a, Ppa1, Hspa5\nIgtp, Ifi47"),
    "Tgd-f" =list(cy="IL-23",                 mk="Cpd, Il22, Gzmb\nSerpinb1a, Ffar2")
  ),
  Naive_CD8_T_cell = list(
    "T8-a"  =list(cy="IFNa1, IFNb\nIFNe, IFNk", mk="Isg15, Ifit3, Ifit1\nBst2, Usp18"),
    "T8-b"  =list(cy="IFNg, IL-2, IL-12\nIL-15, IL-18", mk="Stat1, Ifi27l2a, Igtp\nPsmb9, Gbp4"),
    "T8-c"  =list(cy="IL-1a, IL-1b",          mk="Igfbp4, Gadd45g\nBcl3, Ly6a, Tcf7"),
    "T8-e"  =list(cy="IL-2, IL-15, IL-36a",   mk="Ncl, Ptma, Eif5a\nNpm1, Ran"),
    "T8-f"  =list(cy="IL-4, IL-7",            mk="Eif5a, Ptma, Eif4a1\nNcl, Nme1")
  ),
  cDC2 = list(
    "cDC2-a"=list(cy="IFNa1, IFNb\nIFNe, IFNk", mk="Ifitm3, Isg15, Irf7\nPhf11b, Rtp4"),
    "cDC2-b"=list(cy="IFNg, IL-2, IL-12\nIL-15, IL-18", mk="Cxcl9, Gbp2, Ifi47\nGbp4, Samhd1"),
    "cDC2-c"=list(cy="IL-1a, IL-1b, TNF",     mk="Fabp5, Ccl17, Cish\nVdr, Socs2"),
    "cDC2-e"=list(cy="IL-3, LIF, OSM, CT-1",  mk="Cyp7b1, Cd209d, Tcf4\nApod, Cd209e")
  )
)

# ── Load data ──────────────────────────────────────────────────────────────────
p13_files   <- list.files(DATA_DIR, pattern="_scores_P13\\.csv$")
CT_PREFIXES <- sub("_scores_P13\\.csv$","",p13_files)

read_score <- function(prefix, variant, type="cytokine") {
  suf  <- if (type=="cytokine") paste0("_scores_",variant,"_cytokine.csv")
          else paste0("_scores_",variant,".csv")
  path <- file.path(DATA_DIR, paste0(prefix, suf))
  if (!file.exists(path)) return(NULL)
  df <- read.csv(path, stringsAsFactors=FALSE, check.names=FALSE)
  names(df)[1] <- "label"
  if (type=="cytokine") df$label <- clean_name(df$label)
  df$variant <- variant; df$celltype <- prefix; df
}

cyto_p13  <- bind_rows(Filter(Negate(is.null), lapply(FOCUS, read_score, "P13",  "cytokine")))
cyto_k337 <- bind_rows(Filter(Negate(is.null), lapply(FOCUS, read_score, "K337", "cytokine")))

keep_cyto <- bind_rows(cyto_p13, cyto_k337) %>%
  filter(celltype %in% FOCUS) %>%
  group_by(label) %>%
  summarise(min_padj=min(padj,na.rm=TRUE),.groups="drop") %>%
  filter(min_padj < 0.001) %>% pull(label)
cat("Cytokines:", length(keep_cyto), "\n")

cyto_p13  <- filter(cyto_p13,  label %in% keep_cyto)
cyto_k337 <- filter(cyto_k337, label %in% keep_cyto)

mat_wide <- bind_rows(cyto_p13, cyto_k337) %>%
  rename(ES=`Enrichment Score`) %>%
  group_by(label, celltype) %>%
  summarise(ES=mean(abs(ES)),.groups="drop") %>%
  pivot_wider(names_from=celltype,values_from=ES,values_fill=0) %>%
  tibble::column_to_rownames("label")
col_order <- rownames(mat_wide)[hclust(dist(mat_wide),method="ward.D2")$order]
mat_ct    <- t(mat_wide) %>% as.data.frame()
row_order <- intersect(FOCUS,
  rev(rownames(mat_ct)[hclust(dist(mat_ct),method="ward.D2")$order]))

es_lim <- max(abs(c(cyto_p13$`Enrichment Score`, cyto_k337$`Enrichment Score`)),na.rm=TRUE)
max_sz  <- max(-log10(pmax(c(cyto_p13$padj, cyto_k337$padj),1e-10)),na.rm=TRUE)

# ── BUBBLE PANELS ─────────────────────────────────────────────────────────────
make_panel <- function(df, title_col, show_y, show_leg) {
  df_plot <- df %>%
    rename(Cytokine=label, CellType=celltype, ES=`Enrichment Score`) %>%
    mutate(
      nlp      = pmin(-log10(pmax(padj,1e-10)),10),
      Cytokine = factor(Cytokine, levels=col_order),
      CellType = factor(CT_LABELS[CellType], levels=CT_LABELS[row_order])
    ) %>% filter(!is.na(Cytokine),!is.na(CellType))

  ggplot(df_plot, aes(Cytokine, CellType)) +
    geom_point(aes(size=nlp, fill=ES, colour=ES), shape=21, stroke=0.15, alpha=0.92) +
    scale_size_continuous(name=expression(-log[10](padj)),
      range=c(0.8,10), limits=c(0,max_sz), breaks=c(1,2,3,5),
      guide=if(show_leg) guide_legend(override.aes=list(fill="grey50",colour="grey50",alpha=0.8)) else "none") +
    scale_fill_gradientn(name="ES",
      colours=c("#2166ac","#4393c3","#92c5de","#f7f7f7","#f4a582","#d6604d","#b2182b"),
      limits=c(-es_lim,es_lim),oob=squish,
      guide=if(show_leg) "colorbar" else "none") +
    scale_colour_gradientn(
      colours=c("#2166ac","#4393c3","#92c5de","#f7f7f7","#f4a582","#d6604d","#b2182b"),
      limits=c(-es_lim,es_lim),oob=squish,guide="none") +
    scale_x_discrete(expand=expansion(add=0.5)) +
    scale_y_discrete(expand=expansion(add=0.42)) +
    labs(x=NULL, y=NULL,
         title=if(show_y) "P13 vs REF" else "K337 vs REF") +
    coord_cartesian(clip="off") +
    theme_minimal(base_size=13) +
    theme(
      plot.title      = element_text(size=13,face="bold",hjust=0.5,colour=title_col,margin=margin(b=2)),
      axis.text.x     = element_text(angle=60,hjust=1,vjust=1,size=11,colour="grey20"),
      axis.text.y     = if(show_y) element_text(size=12,colour="grey15") else element_blank(),
      axis.ticks      = element_blank(),
      axis.line.x     = element_line(colour="grey40",linewidth=0.4),
      axis.line.y     = element_line(colour="grey40",linewidth=0.4),
      panel.grid      = element_blank(),
      legend.position = if(show_leg) "right" else "none",
      legend.key.size = unit(0.45,"cm"),
      legend.text     = element_text(size=11),
      legend.title    = element_text(size=11,face="bold"),
      plot.margin     = margin(2,2,0,2)
    )
}

p_p13  <- make_panel(cyto_p13,  COL_P13,  show_y=TRUE,  show_leg=FALSE)
p_k337 <- make_panel(cyto_k337, COL_K337, show_y=FALSE, show_leg=TRUE)

bubble_row <- p_p13 + p_k337 +
  plot_layout(widths=c(1,1), guides="keep")

# ── RADAR PLOTS ────────────────────────────────────────────────────────────────
FLOOR <- 0.12
radar_xy <- function(values, angles)
  data.frame(x=values*cos(angles), y=values*sin(angles))

# Fixed-size annotation boxes in data coordinates.
# Box width is determined per-cell-type by the arc chord at the push radius,
# so boxes never overlap regardless of polygon vertex count.
# Both cy and mk boxes have the same width and height — text centred inside.
BH      <- 0.34   # half-height of each box (same for cy and mk)
GAP     <- 0.4   # radial gap between the two stacked boxes
PUSH_CY <- 4   # radial distance to cy box centre (clear of vertex labels at ~1.22)
PUSH_MK <- PUSH_CY + 2*BH + GAP  # mk box centre directly above cy box

make_radar <- function(ct) {
  d13  <- read_score(ct,"P13", "polar")
  d337 <- read_score(ct,"K337","polar")
  if (is.null(d13)||is.null(d337)) return(ggplot()+theme_void())

  # STOP removed as requested
  states <- sort(union(d13$label, d337$label))
  n      <- length(states)
  angles <- (pi/2) - 2*pi*(seq_len(n)-1)/n

  prep <- function(df) {
    es <- setNames(as.numeric(df$`Enrichment Score`),df$label)
    pv <- setNames(as.numeric(df$padj),df$label)
    es_v <- es[states]; es_v[is.na(es_v)] <- 0
    pv_v <- pv[states]; pv_v[is.na(pv_v)] <- 1
    pos  <- pmax(es_v,0); mx <- max(pos,na.rm=TRUE)
    norm <- if(mx>0) pos/mx else pos
    list(norm=pmax(norm,FLOOR), es=es_v, padj=pv_v)
  }
  r13  <- prep(d13); r337 <- prep(d337)

  mk_poly <- function(vals,var) {
    xy <- radar_xy(vals,angles); xy <- rbind(xy,xy[1,]); xy$variant <- var; xy
  }
  poly_k337 <- mk_poly(r337$norm,"K337")
  poly_p13  <- mk_poly(r13$norm, "P13")
  bg_xy     <- rbind(radar_xy(rep(1,n),angles),   radar_xy(rep(1,1),angles)[1,])
  ring_xy   <- rbind(radar_xy(rep(0.5,n),angles), radar_xy(rep(0.5,1),angles)[1,])
  spoke_df  <- data.frame(x0=0,y0=0,x1=cos(angles),y1=sin(angles))

  push <- 1.22
  lab_df <- data.frame(
    label = states,
    x     = push*cos(angles), y=push*sin(angles),
    hjust = ifelse(cos(angles)>0.15,0,ifelse(cos(angles)<(-0.15),1,0.5)),
    vjust = ifelse(sin(angles)>0.15,0,ifelse(sin(angles)<(-0.15),1,0.5))
  )

  idx_east <- which.min(abs(angles))
  ring_lx  <- 0.5*cos(angles[idx_east])+0.04
  ring_ly  <- 0.5*sin(angles[idx_east])

  mk_dots <- function(r) {
    df2 <- data.frame(x=r$norm*cos(angles),y=r$norm*sin(angles),has_dot=r$es>0)
    filter(df2,has_dot)
  }

  # Per-cell-type safe box half-width: 85% of half the arc chord at PUSH_CY
  BW <- min(
    0.85 * PUSH_CY * sin(pi/n),   # chord-based width at cy push
    0.85 * PUSH_MK * sin(pi/n)    # (also ensure mk boxes don't overlap)
  )

  info  <- POLAR_INFO[[ct]]

  # Find the topmost box position for cell type title placement
  top_mk_y <- max(sapply(seq_along(states), function(i)
    PUSH_MK * sin(angles[i]) + BH
  ))

  p <- ggplot() +
    geom_polygon(data=bg_xy,   aes(x,y), fill="white",colour="steelblue3",linewidth=0.4) +
    geom_polygon(data=ring_xy, aes(x,y), fill=NA, colour="steelblue3",linewidth=0.25,linetype="dashed") +
    geom_segment(data=spoke_df,aes(x=x0,y=y0,xend=x1,yend=y1),colour="steelblue3",linewidth=0.25) +
    annotate("text",x=ring_lx,y=ring_ly,label="0.5",size=3.0,colour="grey55",hjust=0) +
    geom_polygon(data=poly_k337,aes(x,y),fill=COL_K337,colour=COL_K337,alpha=0.18,linewidth=0.9) +
    geom_polygon(data=poly_p13, aes(x,y),fill=COL_P13, colour=COL_P13, alpha=0.18,linewidth=0.9) +
    geom_point(data=mk_dots(r337),aes(x,y),colour=COL_K337,fill=COL_K337,shape=21,size=2.0,stroke=0.7) +
    geom_point(data=mk_dots(r13), aes(x,y),colour=COL_P13, fill=COL_P13, shape=21,size=2.0,stroke=0.7) +
    geom_text(data=lab_df,aes(x,y,label=label,hjust=hjust,vjust=vjust),
              size=3.2,colour="steelblue4",fontface="bold")

  # Draw fixed-size boxes per vertex
  for (i in seq_along(states)) {
    st  <- states[i]; inf <- info[[st]]
    if (is.null(inf)) next
    ax <- cos(angles[i]); ay <- sin(angles[i])
    cx_cy <- PUSH_CY * ax; cy_cy <- PUSH_CY * ay
    cx_mk <- PUSH_MK * ax; cy_mk <- PUSH_MK * ay

    p <- p +
      annotate("rect",
               xmin=cx_cy-BW, xmax=cx_cy+BW,
               ymin=cy_cy-BH,  ymax=cy_cy+BH,
               fill="#E8F5E9", colour="#388E3C", linewidth=0.3) +
      annotate("text", x=cx_cy, y=cy_cy, label=inf$cy,
               size=2.3, colour="#1B5E20", hjust=0.5, vjust=0.5, lineheight=0.88) +
      annotate("rect",
               xmin=cx_mk-BW, xmax=cx_mk+BW,
               ymin=cy_mk-BH,  ymax=cy_mk+BH,
               fill="#E8EAF6", colour="#3949AB", linewidth=0.3) +
      annotate("text", x=cx_mk, y=cy_mk, label=inf$mk,
               size=2.3, colour="#1A237E", hjust=0.5, vjust=0.5, lineheight=0.88)
  }

  # Cell type label above the topmost box
  p +
    annotate("text", x=0, y=top_mk_y+0.22, label=CT_LABELS[ct],
             size=4.0, fontface="bold", colour="grey15", hjust=0.5, vjust=0) +
    coord_fixed(xlim=c(-5.0,5.0), ylim=c(-5.0,5.0), clip="off") +
    labs(title=NULL) +
    theme_void() +
    theme(plot.margin=margin(0,0,0,0))
}

cat("Building radar panels...\n")
radar_panels <- lapply(row_order, make_radar)
radar_row    <- wrap_plots(radar_panels, ncol=length(FOCUS)) &
  theme(plot.margin=margin(0,0,0,0))

# ── ASSEMBLE ───────────────────────────────────────────────────────────────────
n_cyto   <- length(col_order)
n_cells  <- length(row_order)
pw_pdf   <- n_cyto * 0.40 + 7.0
bubble_ht <- n_cells * 0.50 + 1.4
radar_ht  <- 6.2

combined <- bubble_row / radar_row +
  plot_layout(heights=c(bubble_ht, radar_ht))

pdf("IREA_bubble_radar_combined.pdf",
    width=max(20,pw_pdf),
    height=bubble_ht + radar_ht + 0.2,
    useDingbats=FALSE)
print(combined)
dev.off()
cat("Saved Plot4_bubble_radar_combined.pdf\n")
