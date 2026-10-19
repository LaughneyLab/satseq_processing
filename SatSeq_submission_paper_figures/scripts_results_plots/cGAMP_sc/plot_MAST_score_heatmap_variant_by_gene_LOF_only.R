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


heat_mat <- heat_mat[,intersect(lof_df$variant,colnames(heat_mat))]
# -----------------------------
# Annotation
# -----------------------------
anno_data <- facs_weights %>%
  select(variant, ATF3_expression, Dual_IFN_expression,NFkB_expression,
         STING_expression, prolif_expression) %>%
  distinct() %>%
  filter(variant %in% colnames(heat_mat)) %>%
  
  # extract numeric position from variant (e.g. K337 -> 337)
  mutate(position = as.numeric(gsub("^[A-Za-z]+", "", variant))) %>%
  
  # define domains
  mutate(domain = case_when(
    position >= 1   & position <= 150 ~ "Transmembrane",
    position >= 151 & position <= 342 ~ "Ligand-Binding",
    position >= 343 & position <= 378 ~ "C-Terminal Tail",
    TRUE                             ~ "Controls"
  ))

# reorder heatmap columns
heat_mat <- heat_mat[, anno_data$variant]

# -----------------------------
# Colors
# -----------------------------
rd_bu <- rev(brewer.pal(7, "RdBu"))

col_fun_heat <- colorRamp2(seq(-1, 1, length.out = 7), rd_bu)
col_fun_bar  <- colorRamp2(seq(-3, 3, length.out = 7), rd_bu)

domain_cols <- c(
  "Transmembrane"   = "#d95f02",
  "Ligand-Binding"  = "#7570b3",
  "C-Terminal Tail" = "#1b9e77",
  "Controls"           = "grey80"
)

# -----------------------------
# Annotations
# -----------------------------
variant_col_map <- setNames(lof_df$group, lof_df$variant)

# Colors for pathways
pathway_binary_cols <- list(
  UPR  = c("TRUE" = "#35b6e3", "FALSE" = "white"),
  IFN  = c("TRUE" = "#e3312a", "FALSE" = "white"),
  NFkB = c("TRUE" = "#7ac943", "FALSE" = "white")
)

# Ensure lof_df is ordered exactly like your heatmap columns
col_anno_data <- data.frame(variant = colnames(heat_mat)) %>%
  left_join(lof_df, by = "variant")

top_anno <- HeatmapAnnotation(
  # --- New Pathway Binary Indicators ---
  `UPR LOF` = anno_simple(
    as.character(grepl("UPR", col_anno_data$group)), 
    col = pathway_binary_cols$UPR,
    border = F, height = unit(4, "mm")
  ),
  `IFN LOF` = anno_simple(
    as.character(grepl("IFNa", col_anno_data$group)), 
    col = pathway_binary_cols$IFN,
    border = F, height = unit(4, "mm")
  ),
  `NFkB LOF` = anno_simple(
    as.character(grepl("NFkB", col_anno_data$group)), 
    col = pathway_binary_cols$NFkB,
    border = F, height = unit(4, "mm")
  ),
  
  # --- Original Barplots ---
  STING = anno_barplot(anno_data$STING_expression, border = FALSE, 
                       gp = gpar(col=NA, fill = col_fun_bar(anno_data$STING_expression))),
  Proliferation = anno_barplot(anno_data$prolif_expression, border = FALSE, 
                               gp = gpar(col=NA, fill = col_fun_bar(anno_data$prolif_expression))),
  ATF3 = anno_barplot(anno_data$ATF3_expression, border = FALSE, 
                      gp = gpar(col=NA, fill = col_fun_bar(anno_data$ATF3_expression))),
  IFN = anno_barplot(anno_data$Dual_IFN_expression, border = FALSE, 
                     gp = gpar(col=NA, fill = col_fun_bar(anno_data$Dual_IFN_expression))),
  NFkB = anno_barplot(anno_data$NFkB_expression, border = FALSE, 
                      gp = gpar(col=NA, fill = col_fun_bar(anno_data$NFkB_expression))),
  
  
  # Global settings
  gap = unit(c(1, 1, 2, 1, 1, 1, 1), "mm"), # Adjust gaps between specific rows
  show_legend = c(UPR_LOF = FALSE, IFN_LOF = FALSE, NFkB_LOF = FALSE, 
                  STING = TRUE, Proliferation = TRUE, ATF3 = TRUE, NFkB = TRUE, IFN = TRUE)
)


bottom_anno <- HeatmapAnnotation(
  Domain = anno_simple(anno_data$domain, col = domain_cols),
  show_annotation_name = FALSE
)

# Ensure row annotation aligns with heatmap rows
row_anno_data <- gene_pathway_map %>%
  filter(Gene %in% rownames(heat_mat)) %>%
  column_to_rownames("Gene")

# Match row order to heatmap
row_anno_data <- row_anno_data[rownames(heat_mat), , drop = FALSE]

# ==============================================================================
# 4. Colors & Row Annotation
# ==============================================================================

row_anno <- rowAnnotation(
  
  UPR = anno_simple(
    as.character(row_anno_data$`ER Stress`),
    col = pathway_binary_cols$UPR,
    border = F
  ),
  
  IFN = anno_simple(
    as.character(row_anno_data$`IFN Response`),
    col = pathway_binary_cols$IFN,
    border = F
  ),
  
  NFkB = anno_simple(
    as.character(row_anno_data$`NFkB signaling`),
    col = pathway_binary_cols$NFkB,
    border = F
  ),
  
  annotation_name_gp = gpar(fontsize = 10),
  annotation_width = unit(c(4,4,4), "mm"),
  gap = unit(1, "mm")
)

# ==============================================================================
# 5. Updated Heatmap
# ==============================================================================
ht <- Heatmap(
  heat_mat,
  name = "MAST score",
  col = col_fun_heat,
  top_annotation = top_anno,
  bottom_annotation = bottom_anno,
  right_annotation = row_anno,
  cluster_rows = TRUE,
  cluster_columns = T,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8), # Reduced font size for better readability
  column_names_gp = gpar(fontsize = 11)
)

lgd_domain <- Legend(
  labels = names(domain_cols),
  title = "Domain",
  legend_gp = gpar(fill = domain_cols),
  title_position = "leftcenter",
  direction = "horizontal", # Forces items to spread left-to-right
  nrow = 1                  # Explicitly limits to one row
)

saveRDS(ht, file = paste0("cGAMP_gene_variant_MAST_heatmap_", cell_line, ".rds"))

# Open svglite device
svglite(paste0("cGAMP_gene_variant_MAST_heatmap_", cell_line, ".svg"),
        width = 15, height = 11)

# Draw the heatmap
draw(ht, 
     annotation_legend_list = list(lgd_domain), 
     heatmap_legend_side = "right", 
     annotation_legend_side = "bottom")

# Close the device
dev.off()

}
