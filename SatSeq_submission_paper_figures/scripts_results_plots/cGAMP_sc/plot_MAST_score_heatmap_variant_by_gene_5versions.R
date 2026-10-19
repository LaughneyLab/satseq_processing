rm(list=ls()) # Clear environment

library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(patchwork)
library(ggplot2)
library(PCAtools)
library(tidyverse)
library(svglite)

# set working directory to where this script is
setwd(getSrcDirectory(function(){})[1])

# Where the SVGs/PDFs land
out_dir <- "heatmap_svg_versions"
dir.create(out_dir, showWarnings = FALSE)
out_dir_pdf <- "heatmap_pdf_versions"
dir.create(out_dir_pdf, showWarnings = FALSE)

# ============================================================
# Sizing helpers (only affects canvas size, not coloring/aesthetics)
# ============================================================
svg_width  <- function(n_cols) max(8, 2 + n_cols * 0.12)
svg_height <- function(n_rows) max(6, 2 + n_rows * 0.15)

# For the LOF-only plots: keep enough width/height to fit the content
# (all columns labeled, all rows shown), but bend the canvas toward a
# 3:2 (wide:height) aspect ratio rather than the free-form sizing above.
svg_dims_ratio <- function(n_cols, n_rows, ratio = 3/2) {
  content_w <- svg_width(n_cols)
  content_h <- svg_height(n_rows)
  # Whichever dimension the content needs more of "wins"; the other is
  # derived from the ratio so nothing gets clipped.
  if (content_w / content_h >= ratio) {
    list(width = content_w, height = content_w / ratio)
  } else {
    list(width = content_h * ratio, height = content_h)
  }
}

save_heatmap_svg <- function(ht, lgd, file, width, height) {
  svglite::svglite(file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  ComplexHeatmap::draw(ht, heatmap_legend_list = list(lgd), merge_legend = TRUE)
}

# Same plot, saved as a PDF alongside the SVG. `file` is the SVG's file
# path (as passed to save_heatmap_svg); the PDF is written under
# out_dir_pdf with the same base name and a .pdf extension.
save_heatmap_pdf <- function(ht, lgd, file, width, height) {
  pdf_file <- file.path(out_dir_pdf, sub("\\.svg$", ".pdf", basename(file)))
  pdf(pdf_file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  ComplexHeatmap::draw(ht, heatmap_legend_list = list(lgd), merge_legend = TRUE)
}

# ============================================================
# Clustering-order helpers
#
# These reproduce "as if you had clustered X annotation only" by running
# hclust() directly on the small binary annotation matrix (not on the
# MAST-score heatmap matrix), using ComplexHeatmap's own defaults
# (euclidean distance, complete linkage), then reading off the leaf order.
# Used for the row pathway-annotation clustering (plots 1/3/4): "cluster
# the rows of this small binary matrix" where those rows are genes
# (UPR/IFN/NFkB pathway membership).
# ============================================================
get_annotation_cluster_order <- function(bin_mat) {
  bin_mat <- as.matrix(bin_mat)
  storage.mode(bin_mat) <- "numeric"
  hc <- hclust(dist(bin_mat, method = "euclidean"), method = "complete")
  rownames(bin_mat)[hc$order]
}

# Deterministic "how many pathways" grouping for the LOF columns: items
# hit in all 3 pathways come first, then 2-of-3, then 1, then 0. Within
# each count tier, variants are further kept contiguous by their EXACT
# combination of pathways (e.g. all "UPR-only" together, all "IFN+NFkB"
# together, etc.) — so a single-LOF or pairwise-LOF block never gets
# split apart by an unrelated sort key. `tie_order` (residue position) is
# only used as the last, innermost tiebreak among variants that share
# the exact same combination — it can never break up a combination block,
# it just orders the (otherwise arbitrary) rows within one.
get_lof_group_order <- function(bin_df, tie_order) {
  bin_df <- as.data.frame(bin_df)
  mat <- as.matrix(bin_df)
  pathway_count <- rowSums(mat)
  combo_key <- apply(mat, 1, function(r) paste(colnames(mat)[r == 1], collapse = "_"))
  tie_rank <- match(rownames(bin_df), tie_order)
  ord <- order(-pathway_count, combo_key, tie_rank)
  rownames(bin_df)[ord]
}

# ============================================================
# Reusable annotation builders (identical styling to the original script;
# only the data/column order fed in changes between plot versions)
# ============================================================

# Column (variant/position) side: top LOF-pathway strip + expression
# barplots, bottom cell-count barplot + position labels.
build_col_annotations <- function(variant_order, anno_data, lof_df, xgb_data,
                                   col_fun_bar, count_colors,
                                   pathway_binary_cols, cell_line,
                                   show_all_labels = FALSE) {

  anno_sub <- anno_data[match(variant_order, anno_data$variant), , drop = FALSE]

  xgb_sub <- tibble(variant = variant_order) %>%
    left_join(xgb_data, by = "variant") %>%
    mutate(
      count_label_1 = coalesce(count_label_1, 0),
      count_label_0 = coalesce(count_label_0, 0)
    )
  count_matrix_sub <- as.matrix(xgb_sub[, c("count_label_0", "count_label_1")])

  col_anno_sub <- data.frame(variant = variant_order) %>%
    left_join(lof_df, by = "variant")

  top_anno <- if (cell_line == "EO771") {
    NULL
  } else {
    HeatmapAnnotation(
      `UPR LOF` = anno_simple(
        as.character(grepl("UPR", col_anno_sub$group)),
        col = pathway_binary_cols$UPR,
        border = F, height = unit(4, "mm")
      ),
      `IFN LOF` = anno_simple(
        as.character(grepl("IFNa", col_anno_sub$group)),
        col = pathway_binary_cols$IFN,
        border = F, height = unit(4, "mm")
      ),
      `NFkB LOF` = anno_simple(
        as.character(grepl("NFkB", col_anno_sub$group)),
        col = pathway_binary_cols$NFkB,
        border = F, height = unit(4, "mm")
      ),

      STING = anno_barplot(anno_sub$STING_expression, border = FALSE,
                           gp = gpar(col = NA, fill = col_fun_bar(anno_sub$STING_expression))),
      Proliferation = anno_barplot(anno_sub$prolif_expression, border = FALSE,
                                   gp = gpar(col = NA, fill = col_fun_bar(anno_sub$prolif_expression))),
      ATF3 = anno_barplot(anno_sub$ATF3_expression, border = FALSE,
                          gp = gpar(col = NA, fill = col_fun_bar(anno_sub$ATF3_expression))),
      IFN = anno_barplot(anno_sub$Dual_IFN_expression, border = FALSE,
                         gp = gpar(col = NA, fill = col_fun_bar(anno_sub$Dual_IFN_expression))),
      NFkB = anno_barplot(anno_sub$NFkB_expression, border = FALSE,
                          gp = gpar(col = NA, fill = col_fun_bar(anno_sub$NFkB_expression))),

      gap = unit(c(1, 1, 2, 1, 1, 1, 1), "mm"),
      show_legend = c(UPR_LOF = FALSE, IFN_LOF = FALSE, NFkB_LOF = FALSE,
                      STING = TRUE, Proliferation = TRUE, ATF3 = TRUE, NFkB = TRUE, IFN = TRUE)
    )
  }

  # nums is used only to decide WHICH columns get a label (every 5th, or
  # all of them); the label text itself is always the full variant string
  # as it appears in FACS_weights$variant (e.g. "S365", WT residue letter
  # included), never just the bare position number.
  nums <- as.numeric(gsub("[^0-9]", "", variant_order))
  specials <- c("STOP", "dCTT", "RRDD")
  clean_labels <- if (show_all_labels) {
    # LOF-only plots: label every position, not just every 5th
    variant_order
  } else {
    ifelse(
      variant_order %in% specials | (nums %% 5 == 0 & !is.na(nums)),
      variant_order,
      ""
    )
  }
  idx_to_show <- which(clean_labels != "")

  bottom_anno <- HeatmapAnnotation(
    `Cell\ncount` = anno_barplot(
      count_matrix_sub,
      gp = gpar(fill = count_colors, col = NA),
      border = FALSE,
      bar_width = 0.8,
      height = unit(1.2, "cm")
    ),

    `Variant` = anno_mark(
      at = idx_to_show,
      labels = clean_labels[idx_to_show],
      side = "bottom",
      labels_gp = gpar(fontsize = 9, rot = 45)
    ),

    show_annotation_name = TRUE,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10)
  )

  list(top_anno = top_anno, bottom_anno = bottom_anno,
       clean_labels = clean_labels, col_anno_sub = col_anno_sub)
}

# Row (gene) side: pathway indicator strip
build_row_annotation <- function(gene_order, gene_pathway_map, pathway_binary_cols) {
  row_anno_data <- gene_pathway_map %>%
    filter(Gene %in% gene_order) %>%
    column_to_rownames("Gene")
  row_anno_data <- row_anno_data[gene_order, , drop = FALSE]

  row_anno <- rowAnnotation(
    UPR = anno_simple(as.character(row_anno_data$`ER Stress`), col = pathway_binary_cols$UPR, border = F),
    IFN = anno_simple(as.character(row_anno_data$`IFN Response`), col = pathway_binary_cols$IFN, border = F),
    NFkB = anno_simple(as.character(row_anno_data$`NFkB signaling`), col = pathway_binary_cols$NFkB, border = F),
    annotation_name_gp = gpar(fontsize = 10),
    annotation_width = unit(c(4, 4, 4), "mm"),
    gap = unit(1, "mm")
  )

  list(row_anno = row_anno, row_anno_data = row_anno_data)
}

# Heatmap builder: same coloring/legend/fontsize params as the original
# script throughout — only the matrix, annotations, column labels, and
# cluster_rows flag change between versions.
build_variant_heatmap <- function(mat, top_anno, bottom_anno, row_anno,
                                   clean_labels, col_fun_heat, cluster_rows, cell_line,
                                   cluster_columns = FALSE) {
  Heatmap(
    mat,
    column_labels = clean_labels,
    col = col_fun_heat,
    na_col = "gray80",
    top_annotation = top_anno,
    bottom_annotation = bottom_anno,
    right_annotation = row_anno,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 11),
    heatmap_legend_param = list(
      title = paste0(cell_line, " MAST score"),
      title_position = "leftcenter-rot",
      title_gp = gpar(fontsize = 12, fontface = "bold")
    )
  )
}

# ============================================================
# sc-based LOF variant calling (plot 6 only)
#
# Replicates the LOF-calling logic from
# plot_variant_sc_score_FACS_weight_dotplots.R (is_hit_score computed from
# the single-cell MAST pathway scores in BOTH cell lines, combined with
# the matching FACS weight), but keeps ONLY the "is_hit_score & weight <
# -1" branch of that script's `final_hit` formula — the separate "OR
# weight < -2" branch is intentionally dropped per request, so a variant
# only counts as an sc-only LOF hit when both the sc score and the FACS
# weight agree.
#
# This needs both cell lines' sc scores at once (is_hit_score compares
# sign/magnitude across 4T1 AND EO771), so it's computed once, up front,
# rather than inside the per-cell_line loop below.
# ============================================================
STING_df_sc <- rbind(read.csv('4T1_interaction/fcHurdle_STOP.csv'),
                     read.csv('EO771_interaction/fcHurdle_STOP.csv'))
STING_genes_sc <- unique(STING_df_sc$primerid[STING_df_sc$fdr < 0.5])
df_cgamp_stop_sc <- read.csv('4T1_variant_vs_REF_cGAMP/fcHurdle_4T1_STOP.csv')
STING_genes_sc <- intersect(STING_genes_sc, df_cgamp_stop_sc$primerid[df_cgamp_stop_sc$coef < 0])

hallmark_df_sc <- read.csv('mouse_hallmarks_cGAMP_extended.csv')
er_genes_sc <- intersect(STING_genes_sc,
  unique(hallmark_df_sc$Gene[hallmark_df_sc$Pathway %in%
    c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE_CORE")]))
ifna_genes_sc <- intersect(STING_genes_sc,
  unique(hallmark_df_sc$Gene[hallmark_df_sc$Pathway %in% c("HALLMARK_INTERFERON_ALPHA_RESPONSE")]))
nfkb_genes_sc <- intersect(STING_genes_sc,
  unique(hallmark_df_sc$Gene[hallmark_df_sc$Pathway %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB")]))
all_hallmark_genes_sc <- unique(c(er_genes_sc, ifna_genes_sc, nfkb_genes_sc))

load_scores_sc <- function(directory, cell_line) {
  files <- list.files(directory,
                      pattern = paste0("^fcHurdle_", cell_line, "_.*\\.csv$"),
                      full.names = TRUE)
  bind_rows(lapply(files, function(f) {
    variant_name <- gsub(paste0("fcHurdle_", cell_line, "_|.csv"), "", basename(f))
    df <- read.csv(f)
    df$score <- -log10(df$fdr + 1e-300) * sign(df$coef)
    df$variant <- variant_name
    df[, c("primerid", "score", "variant")]
  }))
}

sc_4t1_mat <- load_scores_sc("4T1_variant_vs_REF_cGAMP", "4T1") %>%
  filter(primerid %in% all_hallmark_genes_sc) %>%
  pivot_wider(names_from = variant, values_from = score) %>%
  column_to_rownames("primerid") %>%
  as.matrix()

sc_EO771_mat <- load_scores_sc("EO771_variant_vs_REF_cGAMP", "EO771") %>%
  filter(primerid %in% all_hallmark_genes_sc) %>%
  pivot_wider(names_from = variant, values_from = score) %>%
  column_to_rownames("primerid") %>%
  as.matrix()

get_pathway_avg_sc <- function(genes, label, mat) {
  valid_genes <- intersect(genes, rownames(mat))
  colMeans(mat[valid_genes, , drop = FALSE], na.rm = TRUE) %>%
    enframe(name = "variant", value = label)
}

sc_4T1_avg <- get_pathway_avg_sc(er_genes_sc, "er_stress", sc_4t1_mat) %>%
  left_join(get_pathway_avg_sc(ifna_genes_sc, "ifna_response", sc_4t1_mat), by = "variant") %>%
  left_join(get_pathway_avg_sc(nfkb_genes_sc, "nfkb_response", sc_4t1_mat), by = "variant")

sc_EO771_avg <- get_pathway_avg_sc(er_genes_sc, "er_stress", sc_EO771_mat) %>%
  left_join(get_pathway_avg_sc(ifna_genes_sc, "ifna_response", sc_EO771_mat), by = "variant") %>%
  left_join(get_pathway_avg_sc(nfkb_genes_sc, "nfkb_response", sc_EO771_mat), by = "variant")

sc_4T1_long_sc <- sc_4T1_avg %>%
  pivot_longer(cols = c(er_stress, ifna_response, nfkb_response), names_to = "pathway", values_to = "score_4T1")
sc_EO771_long_sc <- sc_EO771_avg %>%
  pivot_longer(cols = c(er_stress, ifna_response, nfkb_response), names_to = "pathway", values_to = "score_EO771")

facs_weights_sc <- read.csv('FACS_weights_all.csv')

merged_sc_lof <- inner_join(sc_4T1_long_sc, sc_EO771_long_sc, by = c("variant", "pathway")) %>%
  left_join(facs_weights_sc, by = "variant") %>%
  mutate(
    weight = case_when(
      pathway == "nfkb_response" ~ NFkB_expression,
      pathway == "er_stress"     ~ ATF3_expression,
      pathway == "ifna_response" ~ Dual_IFN_expression
    )
  ) %>%
  group_by(pathway) %>%
  mutate(
    mean_4T1   = mean(score_4T1, na.rm = TRUE),
    sd_4T1     = sd(score_4T1, na.rm = TRUE),
    mean_EO771 = mean(score_EO771, na.rm = TRUE),
    sd_EO771   = sd(score_EO771, na.rm = TRUE),
    is_hit_score = (score_4T1 < (mean_4T1 - 1 * sd_4T1) | score_EO771 < (mean_EO771 - 1 * sd_EO771)) &
      (sign(score_4T1) == sign(score_EO771)),
    # Per request: sc-only LOF calling keeps ONLY the "is_hit_score &
    # weight < -1" branch of the uploaded script's final_hit formula —
    # the "OR weight < -2" branch is dropped.
    sc_lof_hit = is_hit_score & (weight < -1)
  ) %>%
  ungroup()

sc_lof_variants <- unique(merged_sc_lof$variant[merged_sc_lof$sc_lof_hit])

# -----------------------------
# USER INPUT
# -----------------------------
for (cell_line in c('4T1','EO771')){

# List of STING-dependent genes which are down in STOP vs REF in 4T1 cGAMP
STING_df <- rbind(read.csv('4T1_interaction/fcHurdle_STOP.csv'),
                  read.csv('EO771_interaction/fcHurdle_STOP.csv'))

STING_genes <- unique(STING_df$primerid[STING_df$fdr < 0.5])

df_cgamp_stop <- read.csv('4T1_variant_vs_REF_cGAMP/fcHurdle_4T1_STOP.csv')
STING_genes <- intersect(STING_genes,df_cgamp_stop$primerid[df_cgamp_stop$coef<0])


# Read Hallmark genes
hallmark_df <- read.csv('mouse_hallmarks_cGAMP_extended.csv')

# ER Stress
er_genes <- intersect(STING_genes,
                      unique(hallmark_df$Gene[hallmark_df$Pathway %in%
                                                c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                                                  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE_CORE"
                                                )]))
# IFNa response
ifna_genes <- intersect(STING_genes,
                        unique(hallmark_df$Gene[hallmark_df$Pathway %in%
                                                  c("HALLMARK_INTERFERON_ALPHA_RESPONSE")]))
# NFkB
nfkb_genes <- intersect(STING_genes,
                        unique(hallmark_df$Gene[hallmark_df$Pathway %in%
                                                  c("HALLMARK_TNFA_SIGNALING_VIA_NFKB")]))

# Combine to get all genes
all_hallmark_genes <- unique(c(er_genes, ifna_genes, nfkb_genes))

# Intersect with STING genes
gene_list <- intersect(STING_genes, all_hallmark_genes)

# Create a dataframe mapping genes to pathways
gene_pathway_map <- tibble(Gene = gene_list) %>%
  mutate(
    `ER Stress` = Gene %in% er_genes,
    `IFN Response` = Gene %in% ifna_genes,
    `NFkB signaling` = Gene %in% nfkb_genes
  ) %>%
  # Create a combined label for genes in multiple pathways, or just one
  mutate(Pathway = case_when(
    `ER Stress` & `IFN Response` & `NFkB signaling` ~ "IFNa & NFkB & UPR",
    `ER Stress` & `IFN Response` ~ "IFNa & UPR",
    `IFN Response` & `NFkB signaling` ~ "IFNa & NFkB",
    `ER Stress` & `NFkB signaling` ~ "ER stress + NFkB",
    `ER Stress` ~ "Only UPR",
    `IFN Response` ~ "Only IFNa",
    `NFkB signaling` ~ "Only NFkB",
    TRUE ~ "Other"
  ))

# -----------------------------
# Load data
# -----------------------------
load_scores <- function(directory, cell_line) {
  files <- list.files(directory,
                      pattern = paste0("^fcHurdle_", cell_line, "_.*\\.csv$"),
                      full.names = TRUE)

  bind_rows(lapply(files, function(f) {
    variant_name <- gsub(paste0("fcHurdle_", cell_line, "_|.csv"), "",
                         basename(f))

    df <- read.csv(f)
    df$score <- -log10(df$fdr + 1e-300) * sign(df$coef)
    df$variant <- variant_name

    df[, c("primerid", "score", "variant")]
  }))
}

df_scores <- load_scores(paste0(cell_line, "_variant_vs_REF_cGAMP"), cell_line)
facs_weights <- read.csv('FACS_weights_all.csv')

# -----------------------------
# Heatmap matrix
# -----------------------------
heat_mat <- df_scores %>%
  filter(primerid %in% gene_list) %>%
  pivot_wider(names_from = variant, values_from = score) %>%
  column_to_rownames("primerid") %>%
  as.matrix()

lof_df <- read.csv('LOF_variants.csv')

# controls are total loss
lof_df <- rbind(lof_df,
                data.frame(
  variant = c("STOP", "dCTT", "RRDD"),
  group = "IFNa & NFkB & UPR",
  residue = NA,
  stringsAsFactors = FALSE
))


# -----------------------------
# Annotation
# -----------------------------
anno_data <- facs_weights %>%
  select(variant, ATF3_expression, Dual_IFN_expression,NFkB_expression,
         STING_expression, prolif_expression) %>%
  distinct() %>%
  # extract numeric position from variant (e.g. K337 -> 337)
  mutate(position = as.numeric(gsub("^[A-Za-z]+", "", variant)))


# Rearrange columns to match anno_data and fill missing variants with NA
heat_mat <- heat_mat[, match(anno_data$variant, colnames(heat_mat))]
colnames(heat_mat) <- anno_data$variant

# 1. Read the XGBoost file
xgb_data <- read_csv(paste0("XGBoost_output_",cell_line,"_cGAMP.csv"))

# -----------------------------
# Colors
# -----------------------------
rd_bu <- rev(brewer.pal(7, "RdBu"))

col_fun_heat <- colorRamp2(seq(-1, 1, length.out = 7), rd_bu)
col_fun_bar  <- colorRamp2(seq(-3, 3, length.out = 7), rd_bu)

count_colors <- c("count_label_0" = "#2166AC", "count_label_1" = "gray80")

variant_col_map <- setNames(lof_df$group, lof_df$variant)

pathway_binary_cols <- list(
  UPR  = c("TRUE" = "#35b6e3", "FALSE" = "white"),
  IFN  = c("TRUE" = "#e3312a", "FALSE" = "white"),
  NFkB = c("TRUE" = "#7ac943", "FALSE" = "white")
)

lgd_counts <- Legend(
  labels = c("STOP-like","REF-like"),
  title = "XGBoost predicted class",
  legend_gp = gpar(fill = count_colors),
  title_position = "leftcenter",
  direction = "vertical",
  nrow = 1
)

# Row annotation for the full/natural gene order (used by plots 1 & 2,
# which keep value-based row clustering — ComplexHeatmap reorders this
# alongside the matrix internally when cluster_rows = TRUE)
row_anno_full <- build_row_annotation(rownames(heat_mat), gene_pathway_map, pathway_binary_cols)
row_anno <- row_anno_full$row_anno
row_anno_data <- row_anno_full$row_anno_data

# ============================================================
# Orders shared across the LOF-only, pathway/LOF-clustered versions
# ============================================================

# Pathway-only row order (plots 3, 4, 5): cluster JUST the row
# pathway-indicator annotation (UPR/IFN/NFkB), not the MAST-score matrix
pathway_bin_mat <- row_anno_data %>%
  transmute(
    UPR  = as.numeric(as.logical(`ER Stress`)),
    IFN  = as.numeric(as.logical(`IFN Response`)),
    NFkB = as.numeric(as.logical(`NFkB signaling`))
  )
rownames(pathway_bin_mat) <- rownames(row_anno_data)
pathway_row_order <- get_annotation_cluster_order(pathway_bin_mat)

# Row annotation reordered to match — reused by plots 1, 3 & 4, all of
# which now order rows by pathway-annotation clustering rather than by
# MAST-score value
row_anno_pathway <- build_row_annotation(pathway_row_order, gene_pathway_map, pathway_binary_cols)$row_anno

# LOF-only column subset (plots 2, 3, 4, 5): variants present in
# LOF_variants.csv (including the STOP/dCTT/RRDD total-LOF controls),
# kept in their original left-to-right column order for plots 2 & 3
lof_variant_set <- unique(lof_df$variant)
lof_columns <- colnames(heat_mat)[colnames(heat_mat) %in% lof_variant_set]

# Residue-position column order (plot 3): ascending residue number,
# STOP/dCTT/RRDD (no residue) placed last
lof_residue_lookup <- lof_df[match(lof_columns, lof_df$variant), c("variant", "residue")]
residue_col_order <- lof_residue_lookup$variant[order(lof_residue_lookup$residue, na.last = TRUE)]

# LOF-pathway-count column order (plot 4): variants hit in all 3 pathways
# grouped to one side, then 2-of-3, then 1, then 0 — mirroring how the
# rows are grouped by pathway membership, so the row and column groupings
# visually correspond.
lof_col_anno <- data.frame(variant = lof_columns) %>% left_join(lof_df, by = "variant")
lof_anno_bin_mat <- data.frame(
  UPR  = as.numeric(grepl("UPR", lof_col_anno$group)),
  IFN  = as.numeric(grepl("IFNa", lof_col_anno$group)),
  NFkB = as.numeric(grepl("NFkB", lof_col_anno$group))
)
rownames(lof_anno_bin_mat) <- lof_col_anno$variant
lof_group_col_order <- get_lof_group_order(lof_anno_bin_mat, tie_order = residue_col_order)

# Common height (gene count never changes across the 5 versions)
h <- svg_height(nrow(heat_mat))

# ============================================================
# Plot 1: existing heatmap — all positions, rows ordered by pathway-
# annotation-only clustering (same row order as plots 3 & 4, for
# consistency across all the heatmap versions)
# ============================================================
heat_mat_pathwayrows <- heat_mat[pathway_row_order, , drop = FALSE]
anno1 <- build_col_annotations(colnames(heat_mat_pathwayrows), anno_data, lof_df, xgb_data,
                               col_fun_bar, count_colors, pathway_binary_cols, cell_line)
ht1 <- build_variant_heatmap(heat_mat_pathwayrows, anno1$top_anno, anno1$bottom_anno, row_anno_pathway,
                             anno1$clean_labels, col_fun_heat, cluster_rows = FALSE, cell_line)
file1 <- file.path(out_dir, sprintf("%s_1_all-positions_rows-pathway-clustered.svg", cell_line))
save_heatmap_svg(ht1, lgd_counts, file1, svg_width(ncol(heat_mat_pathwayrows)), h)
save_heatmap_pdf(ht1, lgd_counts, file1, svg_width(ncol(heat_mat_pathwayrows)), h)

# ============================================================
# Plot 2: same as 1, but LOF positions only (rows still value-clustered)
#
# LOF-only plots: every LOF position gets an x-axis label (show_all_labels
# = TRUE), and the canvas is sized toward a 3:2 wide:height aspect ratio
# rather than the free-form sizing plot 1 uses.
# ============================================================
mat_lof <- heat_mat[, lof_columns, drop = FALSE]
anno2 <- build_col_annotations(lof_columns, anno_data, lof_df, xgb_data,
                               col_fun_bar, count_colors, pathway_binary_cols, cell_line,
                               show_all_labels = TRUE)
ht2 <- build_variant_heatmap(mat_lof, anno2$top_anno, anno2$bottom_anno, row_anno,
                             anno2$clean_labels, col_fun_heat, cluster_rows = TRUE, cell_line)
dims2 <- svg_dims_ratio(ncol(mat_lof), nrow(mat_lof))
file2 <- file.path(out_dir, sprintf("%s_2_LOF-positions-only_rows-value-clustered.svg", cell_line))
save_heatmap_svg(ht2, lgd_counts, file2, dims2$width, dims2$height)
save_heatmap_pdf(ht2, lgd_counts, file2, dims2$width, dims2$height)

# ============================================================
# LOF-only matrix with rows in the same pathway-clustered order as plot 1
# (row_anno_pathway was already built above, shared by plots 1, 3 & 4)
# ============================================================
mat_lof_pathwayrows <- mat_lof[pathway_row_order, , drop = FALSE]

# ============================================================
# Plot 3: LOF positions only, rows pathway-clustered, columns ordered by
# residue position. (This used to also exist as a separate "natural LOF
# column order" version, but for this dataset the natural order already
# is residue order, so that was a duplicate — only one is kept now.)
# ============================================================
mat_pathwayrows_residuecols <- mat_lof_pathwayrows[, residue_col_order, drop = FALSE]
anno3 <- build_col_annotations(residue_col_order, anno_data, lof_df, xgb_data,
                               col_fun_bar, count_colors, pathway_binary_cols, cell_line,
                               show_all_labels = TRUE)
ht3 <- build_variant_heatmap(mat_pathwayrows_residuecols, anno3$top_anno, anno3$bottom_anno, row_anno_pathway,
                             anno3$clean_labels, col_fun_heat, cluster_rows = FALSE, cell_line)
dims3 <- svg_dims_ratio(ncol(mat_pathwayrows_residuecols), nrow(mat_pathwayrows_residuecols))
file3 <- file.path(out_dir, sprintf("%s_3_LOF-positions-only_rows-pathway-clustered_cols-residue-order.svg", cell_line))
save_heatmap_svg(ht3, lgd_counts, file3, dims3$width, dims3$height)
save_heatmap_pdf(ht3, lgd_counts, file3, dims3$width, dims3$height)

# ============================================================
# Plot 4: same rows as 3, but columns grouped by LOF-pathway count — all
# 3 pathways to one side, then 2-of-3, then 1, then 0 — so the column
# grouping visually echoes the pathway-based row grouping.
# ============================================================
mat_pathwayrows_lofgroupcols <- mat_lof_pathwayrows[, lof_group_col_order, drop = FALSE]
anno4 <- build_col_annotations(lof_group_col_order, anno_data, lof_df, xgb_data,
                               col_fun_bar, count_colors, pathway_binary_cols, cell_line,
                               show_all_labels = TRUE)
ht4 <- build_variant_heatmap(mat_pathwayrows_lofgroupcols, anno4$top_anno, anno4$bottom_anno, row_anno_pathway,
                             anno4$clean_labels, col_fun_heat, cluster_rows = FALSE, cell_line)
dims4 <- svg_dims_ratio(ncol(mat_pathwayrows_lofgroupcols), nrow(mat_pathwayrows_lofgroupcols))
file4 <- file.path(out_dir, sprintf("%s_4_LOF-positions-only_rows-pathway-clustered_cols-LOF-pathway-count-grouped.svg", cell_line))
save_heatmap_svg(ht4, lgd_counts, file4, dims4$width, dims4$height)
save_heatmap_pdf(ht4, lgd_counts, file4, dims4$width, dims4$height)

# ============================================================
# Plot 5 (new): LOF positions only - same column order as plot 4
# (lof_group_col_order: LOF-pathway-count grouping), but rows clustered by
# the MAST-score heatmap values themselves (like plot 2, NOT the
# pathway-annotation row order used by plots 1/3/4).
# ============================================================
mat_lof_colclustered <- mat_lof[, lof_group_col_order, drop = FALSE]
anno5 <- build_col_annotations(lof_group_col_order, anno_data, lof_df, xgb_data,
                               col_fun_bar, count_colors, pathway_binary_cols, cell_line,
                               show_all_labels = TRUE)
ht5 <- build_variant_heatmap(mat_lof_colclustered, anno5$top_anno, anno5$bottom_anno, row_anno,
                             anno5$clean_labels, col_fun_heat, cluster_rows = TRUE, cell_line)
dims5 <- svg_dims_ratio(ncol(mat_lof_colclustered), nrow(mat_lof_colclustered))
file5 <- file.path(out_dir, sprintf("%s_5_LOF-positions-only_cols-LOF-pathway-count-grouped_rows-value-clustered.svg", cell_line))
save_heatmap_svg(ht5, lgd_counts, file5, dims5$width, dims5$height)
save_heatmap_pdf(ht5, lgd_counts, file5, dims5$width, dims5$height)

# ============================================================
# Plot 6 (new): sc-only LOF positions - variants called as LOF purely from
# the single-cell MAST-score + FACS-weight criterion computed above
# (sc_lof_variants), NOT the LOF_variants.csv set used by plots 2-5. Both
# rows AND columns are clustered directly on the MAST-score heatmap values
# (cluster_rows = TRUE, cluster_columns = TRUE) - no pathway-annotation or
# manual grouping order for either axis.
# ============================================================
sc_lof_columns <- colnames(heat_mat)[colnames(heat_mat) %in% sc_lof_variants]
mat_sc_lof <- heat_mat[, sc_lof_columns, drop = FALSE]
anno6 <- build_col_annotations(sc_lof_columns, anno_data, lof_df, xgb_data,
                               col_fun_bar, count_colors, pathway_binary_cols, cell_line,
                               show_all_labels = TRUE)
ht6 <- build_variant_heatmap(mat_sc_lof, anno6$top_anno, anno6$bottom_anno, row_anno,
                             anno6$clean_labels, col_fun_heat, cluster_rows = TRUE, cell_line,
                             cluster_columns = TRUE)
dims6 <- svg_dims_ratio(ncol(mat_sc_lof), nrow(mat_sc_lof))
file6 <- file.path(out_dir, sprintf("%s_6_sc-only-LOF-positions_rows-cols-value-clustered.svg", cell_line))
save_heatmap_svg(ht6, lgd_counts, file6, dims6$width, dims6$height)
save_heatmap_pdf(ht6, lgd_counts, file6, dims6$width, dims6$height)

cat("Wrote 6 SVGs + 6 PDFs for", cell_line, "to", out_dir, "and", out_dir_pdf, "\n")

}
