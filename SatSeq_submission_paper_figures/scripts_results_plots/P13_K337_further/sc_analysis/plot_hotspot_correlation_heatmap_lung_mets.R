rm(list=ls())#Clear environment
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# 1. Load Required Libraries Safely
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(cowplot)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(magick)
library(grid)
library(svglite)

# Define Paths relative to your project workspace
dist_matrix_path <- "hotspot_gene_distances_lung_mets.csv"
modules_path     <- "hotspot_gene_modules_lung_mets.csv"
output_pdf_path  <- "hotspot_heatmap.pdf"

# 2. Load the Dataset
dist_matrix <- as.matrix(read.csv(dist_matrix_path, row.names = 1))
modules_df  <- read.csv(modules_path, row.names = 1)

# Ensure the modules dataframe matches the matrix rows initially
all_matrix_genes <- rownames(dist_matrix)
modules_df       <- modules_df[all_matrix_genes, , drop = FALSE]

# Filter out module "-1" while preserving original matrix row/col order
keep_logical <- as.character(modules_df$Module) != "-1" & 
  all_matrix_genes %in% colnames(dist_matrix)

dist_matrix <- dist_matrix[keep_logical, all_matrix_genes[keep_logical], drop = FALSE]
modules_df  <- modules_df[keep_logical, , drop = FALSE]
filtered_genes <- rownames(dist_matrix)

# ==============================================================================
# NEW STEP: Remap numeric modules to biological group names
# ==============================================================================
modules_df <- modules_df %>%
  mutate(Module_Str = as.character(Module)) %>%
  mutate(Functional_Module = paste0("Module ",Module_Str))

# Extract the named vector for mapping
module_vector <- modules_df$Functional_Module
# ==============================================================================

# 3. Clean and Explicit Gene Selection List
target_hallmarks <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_UV_RESPONSE_DN",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")

# Read the GMT file line by line
gmt_lines <- readLines("mouse_hallmarks.gmt")

# Split each line by tabs to turn it into a list of elements
parsed_gmt <- strsplit(gmt_lines, "\t")

# Filter the list to keep only elements where the first item (pathway name) matches your targets
matched_lines <- parsed_gmt[sapply(parsed_gmt, function(x) x[1] %in% target_hallmarks)]

# Extract genes (skipping columns 1 and 2) and flatten into a single vector
all_extracted_genes <- unlist(lapply(matched_lines, function(x) x[3:length(x)]))

# Clean up empty strings and get unique gene names
gene_list <- unique(all_extracted_genes[all_extracted_genes != ""])

# Spike in some proliferation genes not included
gene_list <- c(gene_list,
               'Mki67',"Myc","Itga2", "Diaph3","Pcna","Mre11a")

# Identify pointer indices on the updated row sequence
genes_to_label <- intersect(gene_list, filtered_genes)
label_indices  <- which(filtered_genes %in% genes_to_label)

# 4. Color Scale Configurations
# Heatmap interior matrix scale (-100 to 100)
brewer_base_colors <- rev(brewer.pal(11, "RdBu"))
vibrant_rdbu_ramp  <- colorRampPalette(brewer_base_colors)(101)

heatmap_col_fun <- colorRamp2(
  seq(-100, 100, length.out = 101), 
  vibrant_rdbu_ramp
)

# Shared color specifications for the named biological categories
module_colors <- c(
  "Module 1" = "#e41a1c", 
  "Module 2" = "#377eb8", 
  "Module 3" = "#4daf4a", 
  "Module 4" = "#984ea3", 
  "Module 5" = "#ff7f00")

# 5. Generate Symmetric Axis Annotations
# Row annotation track (Left side)
row_ann <- rowAnnotation(
  Module = module_vector,
  col = list(Module = module_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Module = list(
      title_gp  = gpar(fontsize = 20, fontface = "bold",
                       title_position = "lefttop-rot"), # Legend title size ("Module")
      labels_gp = gpar(fontsize = 20)  # Module names/categories font size
    )
  )
)

# Column annotation track (Top side)
col_ann <- columnAnnotation(
  Module = module_vector,
  col = list(Module = module_colors),
  show_legend = FALSE,         # Hides duplicate legend for clean layout
  show_annotation_name = FALSE
)

# Text annotation track (Right side pointer links)
text_ann <- rowAnnotation(
  Genes = anno_mark(
    at = label_indices, 
    labels = filtered_genes[label_indices],
    labels_gp = gpar(fontface = "italic", fontSize = 6.5)))

# =========================================================================
# CRITICAL FIX: Use a temporary graphics device to lock anno_mark layouts
# =========================================================================

# Open a temporary device matching your target canvas width/height aspect ratio
pdf(tf <- tempfile(), width = 15, height = 14)

grid::grid.newpage()
ht <- Heatmap(
  dist_matrix,
  name = "Z-score",
  col = heatmap_col_fun,
  cluster_rows = FALSE,      
  cluster_columns = FALSE,   
  show_row_names = FALSE,    
  show_column_names = FALSE, 
  left_annotation = row_ann,   # Color blocks matching rows
  top_annotation = col_ann,    # Identical color blocks matching columns
  right_annotation = text_ann,  # Pointer lines for targeted labels
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 16), # Scale text size
    title_gp  = gpar(fontsize = 20, fontface = "bold",
                     title_position = "lefttop-rot" )  # Legend title (Z-score) size
  )
)

# Draw the heatmap inside the active device so coordinates evaluate correctly
draw(ht)

# Grab the fully compiled graphical object tree (grob)
grob_ht <- grid::grid.grab()

# Close the temporary device and clean up the file immediately
dev.off()
if(file.exists(tf)) file.remove(tf)

# =========================================================================
# Wrap into Patchwork and Assemble Full Figure
# =========================================================================

# Wrap the frozen grob and eliminate outer margins so it stacks tightly
ht_grob <- wrap_elements(grob_ht) + 
  theme(
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"), # Add 15pt cushion at bottom
    plot.background = element_blank(),
    panel.background = element_blank()
  )

# Save the final consolidated multi-panel PDF
pdf('hotspot_module_gene_correlation_lung_mets.pdf', width = 23, height = 20)
grid::grid.draw(ht_grob)
dev.off()



