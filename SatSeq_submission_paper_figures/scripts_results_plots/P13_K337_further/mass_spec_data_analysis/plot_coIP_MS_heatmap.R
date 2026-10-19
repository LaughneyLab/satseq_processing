rm(list=ls()) # Clear environment
setwd(getSrcDirectory(function(){})[1]) # Set working directory to script location

# Load required libraries
library(readxl)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(circlize)
library(svglite)
library(scales) # Needed for squish

# 1. Load the data
df <- read_excel("SK-1-71_For_Amelia_1.xlsx", sheet = "Sheet1")

# 2. Define specific genes from the screenshot
target_genes <- c("Anxa2", "Atp2a2", "Copa", "Copb2", "Copg1", "Ifit1", "Sec23ip", 
  "Tns3", "Traf3", "Tuba1b", "Tubb2a", "Tubb4b", "Ubap2", 
  "Ubc", "Vdac1", "Vdac2", "Vdac3", "Vim", "Vps16")

# 3. Define time point mapping (excluding chronic)
tp_mapping <- c(
  "0"   = "0h",
  "0.5" = "30m",
  "4"   = "4h",
  "10"  = "10h",
  "24"  = "24h"
)

# 4. Clean and reshape data (keeping separate values for REF, P13A, K337R)
processed_data <- df %>%
  select(Gene_Symbol = `Gene Symbol`, everything(), -Description) %>%
  filter(Gene_Symbol %in% target_genes) %>%
  mutate(across(-Gene_Symbol, ~ as.numeric(as.character(.)))) %>%
  pivot_longer(
    cols = -Gene_Symbol,
    names_to = c("Condition", "Time"),
    names_pattern = "([A-Za-z0-9]+)_(.*)"
  ) %>%
  filter(Time %in% names(tp_mapping)) %>% # Excludes chronic time point
  group_by(Gene_Symbol, Condition, Time) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Time_Label = factor(tp_mapping[Time], levels = c("0h", "30m", "4h", "10h", "24h")),
    Condition = factor(Condition, levels = c("REF", "P13A", "K337R"))
  ) %>%
  arrange(Condition, Time_Label)

# 5. Prepare Matrix (Rows = Genes, Columns = Variant_Time)
# Create a unique column key combining Condition and Time
wide_data <- processed_data %>%
  mutate(Col_Name = paste(Condition, Time_Label, sep = "_")) %>%
  select(Gene_Symbol, Col_Name, value) %>%
  pivot_wider(names_from = Col_Name, values_from = value)

# Set Row Names and convert to Matrix
heatmap_matrix <- wide_data %>%
  column_to_rownames("Gene_Symbol") %>%
  as.matrix()

# Ensure the rows follow the exact order of target_genes list
heatmap_matrix <- heatmap_matrix[match(target_genes,
                                       rownames(heatmap_matrix)), ]

# 6. Create Column Top Annotation (Condition/Variant blocks)
col_conditions <- factor(
  rep(c("REF", "P13A", "K337R"), each = 5), 
  levels = c("REF", "P13A", "K337R"))

col_time_labels <- rep(c("0h", "30m", "4h", "10h", "24h"), 3)

top_annotation <- HeatmapAnnotation(
  Condition = col_conditions,
  col = list(Condition = c("REF" = "#4E4E4E",
                           "P13A" = "#FBB040",
                           "K337R" = "#E43B3F")),
  show_annotation_name = FALSE
)

# 7. Establish exact visual scale mapping constraints (5 to 15)
heatmap_colors <- colorRampPalette(rev(brewer.pal(7, "RdBu")))(50) # Adjusted palette for absolute scale intensity
color_mapping <- colorRamp2(seq(5, 15, length.out = 50), heatmap_colors)

# 8. Generate the final heatmap
p <- Heatmap(
  heatmap_matrix,
  col = color_mapping,
  border = TRUE,
  top_annotation = top_annotation,
  
  # Dendrogram and Label placement configurations
  cluster_rows = TRUE,      
  cluster_columns = FALSE,   # Keeps timecourse chronological within variants
  
  column_split = col_conditions,
  column_gap = unit(3, "mm"), 
  column_title = NULL,
  column_labels = col_time_labels,
  
  column_names_gp = gpar(fontsize = 18),
  row_names_gp = gpar(fontsize = 18),
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = NULL,
    title_gp = gpar(fontsize = 18),
    labels_gp = gpar(fontsize = 18) 
  )
)

# Save output plot
svglite("coIP-MS_heatmap.svg", width = 12, height = 8)
draw(p)
dev.off()


# 5. Calculate log2 fold change relative to REF
diff_data <- processed_data %>%
  pivot_wider(names_from = Condition, values_from = value) %>%
  mutate(
    P13A_vs_REF = P13A - REF,
    K337R_vs_REF = K337R - REF
  ) %>%
  select(Gene_Symbol, Time_Label, P13A_vs_REF, K337R_vs_REF)

# 6. Prepare Matrix Blocks
time_levels <- c("0h", "30m", "4h", "10h", "24h")

matrix_p13a <- diff_data %>%
  select(Gene_Symbol, Time_Label, P13A_vs_REF) %>%
  pivot_wider(names_from = Time_Label, values_from = P13A_vs_REF) %>%
  column_to_rownames("Gene_Symbol") %>%
  as.matrix()
matrix_p13a <- matrix_p13a[match(target_genes, rownames(matrix_p13a)), time_levels]

matrix_k337r <- diff_data %>%
  select(Gene_Symbol, Time_Label, K337R_vs_REF) %>%
  pivot_wider(names_from = Time_Label, values_from = K337R_vs_REF) %>%
  column_to_rownames("Gene_Symbol") %>%
  as.matrix()
matrix_k337r <- matrix_k337r[match(target_genes, rownames(matrix_k337r)), time_levels]

# 7. Generate a UNIFIED clustering dendrogram across both datasets
# Concatenating the matrices side-by-side calculates distance using all data points simultaneously
combined_matrix_for_clustering <- cbind(matrix_p13a, matrix_k337r)
shared_dist <- dist(combined_matrix_for_clustering)
shared_dendrogram <- as.dendrogram(hclust(shared_dist))

# 8. Establish Shrunk Colorbar Range with OOB Squish Mapping
# Define your narrowed limits (e.g., -1.5 to 1.5 log2 Fold Change)
color_limit <- 5 

# Squish the matrices so values outside limits gracefully take on the maximum scale color 
matrix_p13a_squished  <- scales::squish(matrix_p13a,  range = c(-color_limit, color_limit))
matrix_k337r_squished <- scales::squish(matrix_k337r, range = c(-color_limit, color_limit))

# Set color mapping anchored precisely inside our narrowed limit bounds
heatmap_colors <- colorRampPalette(rev(brewer.pal(7, "RdBu")))(50)
color_mapping <- colorRamp2(seq(-color_limit, color_limit, length.out = 50), heatmap_colors)

# 9. Define block top annotations
top_ann_p13a <- HeatmapAnnotation(
  Condition = factor(rep("P13A", 5), levels = c("P13A")),
  col = list(Condition = c("P13A" = "#FBB040")),
  show_annotation_name = FALSE
)

top_ann_k337r <- HeatmapAnnotation(
  Condition = factor(rep("K337R", 5), levels = c("K337R")),
  col = list(Condition = c("K337R" = "#E43B3F")),
  show_annotation_name = FALSE
)

# 10. Generate individual heatmap blocks passing the single shared dendrogram
p_p13a <- Heatmap(
  matrix_p13a_squished,
  col = color_mapping,
  border = TRUE,
  top_annotation = top_ann_p13a,
  cluster_rows = shared_dendrogram, # Applies identical gene hierarchy tree
  cluster_columns = FALSE,
  row_dend_width = unit(15, "mm"),   # Shows the shared tree on the far left side
  column_title = "P13A vs REF",
  column_title_gp = gpar(fontsize = 18),
  column_labels = c("0h", "30m", "4h", "10h", "24h"),
  column_names_gp = gpar(fontsize = 18),
  row_names_gp = gpar(fontsize = 18),
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = "logFC",
    title_gp = gpar(fontsize = 18),
    labels_gp = gpar(fontsize = 18)
  )
)

p_k337r <- Heatmap(
  matrix_k337r_squished,
  col = color_mapping,
  border = TRUE,
  top_annotation = top_ann_k337r,
  cluster_rows = shared_dendrogram, # Applies identical gene hierarchy tree
  cluster_columns = FALSE,
  show_row_dend = FALSE,             # Hide duplicate tree on second block
  column_title = "K337R vs REF",
  column_title_gp = gpar(fontsize = 18),
  column_labels = c("0h", "30m", "4h", "10h", "24h"),
  column_names_gp = gpar(fontsize = 18),
  row_names_gp = gpar(fontsize = 18),
  show_row_names = TRUE,
  row_names_side = "right", 
  show_column_names = TRUE,
  column_names_rot = 45,
  show_heatmap_legend = FALSE 
)

# 11. Combine plots side-by-side
combined_heatmap <- p_p13a + p_k337r

# Save output plot
svglite("coIP-MS_diff_heatmap.svg", width = 8, height = 8)
draw(combined_heatmap, ht_gap = unit(10, "mm"))
dev.off()

