rm(list=ls()) # Clear environment
setwd(getSrcDirectory(function(){})) # Set working directory to script location

library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(svglite)
library(ggplot2)
library(cowplot)
library(patchwork)

# ==============================================================================
# 1. PARAMETERS & INPUT CONFIGURATION
# ==============================================================================
target_genes <- read.csv('STING_pathway_cGAMP_genes.csv')[,2]
hallmark_df <- read.csv('../../cGAMP_sc/mouse_hallmarks_cGAMP_extended.csv')
ifn_genes <- hallmark_df$Gene[hallmark_df$Pathway=='HALLMARK_INTERFERON_ALPHA_RESPONSE']
#target_genes <- intersect(target_genes,ifn_genes)

rep1_tsv <- "7J9GL7_results/7J9GL7-expression-matrix.tsv"
rep1_mapping_dir <- "7J9GL7_results/7J9GL7_mapping-stats"

rep2_tsvs <- c("3F5XJK_results/3F5XJK-expression-matrix.tsv", "Q57ZYK_results/Q57ZYK-expression-matrix.tsv")
rep2_mapping_dirs <- c("3F5XJK_results/3F5XJK_mapping-stats", "Q57ZYK_results/Q57ZYK_mapping-stats")

# ==============================================================================
# 2. EXTRACT AND PARSE REPLICATE 1 VALUES
# ==============================================================================
files_rep1 <- list.files(rep1_mapping_dir, pattern = "\\.tsv$", full.names = FALSE)

metadata_rep1 <- data.frame(filename = files_rep1, stringsAsFactors = FALSE) %>%
  mutate(
    sample_id = str_extract(filename, "^[A-Z0-9]+_\\d+"),
    Variant = str_split_i(filename, "_", 4),
    Time_raw = str_split_i(str_remove(filename, "\\.tsv$"), "_", 5)
  ) %>%
  filter(Time_raw != "chronic") %>% 
  mutate(
    Time = paste0(Time_raw, "h"),
    raw_colname = paste0("X", sample_id, "_cpm")
  ) %>%
  select(raw_colname, Variant, Time)

df_rep1 <- read.delim(rep1_tsv, sep = "\t",
                      stringsAsFactors = FALSE, check.names = TRUE)

rep1_processed <- df_rep1 %>% 
  filter(tolower(gene_name) %in% tolower(target_genes)) %>%
  group_by(gene_name) %>%
  summarise(across(ends_with("_cpm"), \(x) sum(x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(cols = ends_with("_cpm"), names_to = "raw_colname", values_to = "CPM") %>%
  inner_join(metadata_rep1, by = "raw_colname") 
# ==============================================================================
# 3. EXTRACT AND PARSE REPLICATE 2 VALUES
# ==============================================================================
files_rep2 <- unlist(lapply(rep2_mapping_dirs, function(d) list.files(d, pattern = "\\.tsv$", full.names = FALSE)))

metadata_rep2 <- data.frame(filename = files_rep2, stringsAsFactors = FALSE) %>%
  mutate(
    # Captures the single-index pattern: e.g., "3F5XJK_1" or "Q57ZYK_10"
    sample_id = str_extract(filename, "^[A-Z0-9]+_\\d+"),
    Variant = str_split_i(filename, "_", 3),
    # Captures trailing numeric digits or explicit text strings safely
    Time_raw = str_extract(str_remove(filename, "\\.tsv$"), "(chronic|\\d+)$")
  ) %>%
  filter(Time_raw != "chronic") %>% 
  mutate(
    Time = paste0(Time_raw, "h"),
    prefix = ifelse(str_detect(sample_id, "^[0-9]"), "X", ""),
    raw_colname = paste0(prefix, sample_id, "_cpm")
  ) %>%
  select(raw_colname, Variant, Time)

rep2_processed_list <- lapply(rep2_tsvs, function(p) {
  df_block <- read.delim(p, sep = "\t", stringsAsFactors = FALSE, check.names = TRUE)
  cpm_cols <- grep("_cpm$", colnames(df_block), value = TRUE)
  if(length(cpm_cols) == 0) return(NULL)
  
  df_block %>%
    filter(tolower(gene_name) %in% tolower(target_genes)) %>%
    group_by(gene_name) %>%
    summarise(across(all_of(cpm_cols), \(x) sum(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(cols = all_of(cpm_cols), names_to = "raw_colname", values_to = "CPM") %>%
    inner_join(metadata_rep2, by = "raw_colname") 
})

rep2_processed <- bind_rows(rep2_processed_list[!sapply(rep2_processed_list, is.null)])


# 1. Combine replicates and average CPM
combined_reps <- bind_rows(rep1_processed, rep2_processed) %>%
  group_by(gene_name, Variant, Time) %>%
  summarize(mean_CPM = mean(CPM, na.rm = TRUE), .groups = "drop") %>%
  mutate(log2_CPM = log2(mean_CPM + 1))

# 2. Extract REF values to calculate log2 fold change relative to REF at each Time
ref_data <- combined_reps %>%
  filter(Variant == "REF") %>%
  select(gene_name, Time, ref_log2_CPM = log2_CPM)

combined_data <- combined_reps %>%
  inner_join(ref_data, by = c("gene_name", "Time")) %>%
  mutate(
    final_val = log2_CPM - ref_log2_CPM,
    sample_group_id = paste(Variant, Time, sep = "_")
  )

# 3. Apply your custom setup and aesthetics
time_order <- c("0h", "2h", "6h", "12h")
custom_colors <- c("REF" = "#4E4E4E", "STOP" = "#B4B4B4",
                   "P13A" = "#FBB040", "K337R" = "#E43B3F")

metadata_final <- combined_data %>% 
  select(sample_group_id, Variant, Time) %>% 
  distinct() %>% 
  mutate(
    Variant = factor(Variant, levels = names(custom_colors)),
    Time = factor(Time, levels = time_order)
  ) %>% 
  arrange(Variant, Time) %>%
  filter(Variant != "REF") # Exclude REF from the final plot metadata


wide_df <- combined_data %>% 
  select(gene_name, sample_group_id, final_val) %>% 
  pivot_wider(names_from = sample_group_id, values_from = final_val)

gene_names_vector <- wide_df$gene_name
wide_df$gene_name <- NULL
mat_raw <- as.matrix(wide_df)
rownames(mat_raw) <- gene_names_vector
class(mat_raw) <- "numeric"

ordered_groups <- metadata_final$sample_group_id
mat_ordered <- mat_raw[, ordered_groups, drop = FALSE]

# 4. Transpose matrix and flip annotations to make the heatmap horizontal 
# (Genes on X-axis, Time-Variant groups on Y-axis)
mat_horizontal <- t(mat_ordered)

annotation_base <- as.data.frame(metadata_final)
rownames(annotation_base) <- annotation_base$sample_group_id
annotation_base$sample_group_id <- NULL

variant_vector <- as.character(annotation_base$Variant)
gap_indices <- which(variant_vector[-1] != variant_vector[-length(variant_vector)])

time_colors <- c("0h" = "#E4D00A", "2h" = "#A1D99B", "6h" = "#74C476", "12h" = "#238443")
annotation_colors <- list(Variant = custom_colors, Time = time_colors)
heatmap_colors <- colorRampPalette(rev(brewer.pal(7, "RdBu")))(50)

palette_breaks <- seq(-0.5, 0.5, length.out = 51)

# 5. Generate horizontal heatmap (swapped clustering, gaps, and label visibility)
p_heatmap <- pheatmap(
  mat_horizontal,
  annotation_row = annotation_base,      # Changed from annotation_col
  annotation_colors = annotation_colors,
  color = heatmap_colors,
  breaks = palette_breaks,
  cluster_rows = FALSE,                  # Do not cluster the ordered sample groups
  cluster_cols = TRUE,                   # Cluster the genes on the X-axis
  gaps_row = gap_indices,                # Changed from gaps_col
  scale = "none",
  show_colnames = TRUE,                  # Show gene names on X-axis
  show_rownames = FALSE,                  # Show sample groups on Y-axis
  border_color = "white"
)

# Inverted width and height for horizontal aspect ratio
svglite("rnaseq_heatmap_logFC_from_REF_full.svg", width = 10, height = 8)
print(p_heatmap)
dev.off()