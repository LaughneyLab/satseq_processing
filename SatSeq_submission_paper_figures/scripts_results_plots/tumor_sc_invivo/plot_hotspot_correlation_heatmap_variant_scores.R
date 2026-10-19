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

# 1. Load the dataset
df <- read.csv("hotspot_module_scores_tumor_4T1.csv", check.names = FALSE)

# Fix for empty/NA column names
names(df)[names(df) == "" | is.na(names(df))] <- paste0("v", seq_along(names(df)[names(df) == "" | is.na(names(df))]))

# 2. Pivot the dataset into long format
df_long <- df %>%
  pivot_longer(
    cols = c("1", "2", "3"),
    names_to = "pathway_id",
    values_to = "score"
  ) %>%
  mutate(
    pathway = case_when(
      pathway_id == "1" ~ "Proliferation",
      pathway_id == "2" ~ "IFN",
      pathway_id == "3" ~ "NFkB"
    ),
    pathway = factor(pathway, levels = c("Proliferation", "IFN", "NFkB"))
  )

# 3. Process WT data (Assigned to position 0)
df_wt <- df_long %>%
  filter(str_detect(variant_substitution, "mSTING_WT")) %>%
  mutate(position = 0)

# Extract precise baseline mean scalars per pathway
mean_prolif <- mean(df_wt$score[df_wt$pathway == "Proliferation"], na.rm = TRUE)
mean_ifn    <- mean(df_wt$score[df_wt$pathway == "IFN"], na.rm = TRUE)
mean_nfkb   <- mean(df_wt$score[df_wt$pathway == "NFkB"], na.rm = TRUE)

# 4. Process Variant data (Positions 1:378)
df_variants <- df_long %>%
  filter(!str_detect(variant_substitution, "mSTING_WT")) %>%
  mutate(position = as.numeric(str_extract(variant_substitution, "\\d+"))) %>%
  filter(!is.na(position) & position >= 1 & position <= 378)

# 5. Compute statistics per position (Wilcoxon, mean, max)
position_stats <- df_variants %>%
  group_by(pathway, position) %>%
  summarize(
    max_score = max(score, na.rm = TRUE),
    local_mean = mean(score, na.rm = TRUE),
    p_value = tryCatch({
      current_scores <- score
      wt_reference <- df_wt$score[df_wt$pathway == first(pathway)]
      wilcox.test(current_scores, wt_reference)$p.value
    }, error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(
    stars = ifelse(!is.na(p_value) & p_value < 0.01, "*", ""),
    position_factor = factor(position, levels = c(0, 1:378))
  )

# Add position 0 (WT) mean statistics to the tracking dataset for drawing mean bars
wt_stats <- df_wt %>%
  group_by(pathway) %>%
  summarize(
    max_score = max(score, na.rm = TRUE),
    local_mean = mean(score, na.rm = TRUE),
    stars = "",
    .groups = "drop"
  ) %>%
  mutate(position = 0, position_factor = factor(0, levels = c(0, 1:378)))

all_stats <- bind_rows(wt_stats, position_stats)
significance_data <- all_stats %>% filter(stars == "*")

# 6. Combine WT and Variants into a single unified dataframe for plotting
df_plot <- bind_rows(df_wt, df_variants) %>%
  mutate(position_factor = factor(position, levels = c(0, 1:378)))

# Fix Cropping: Find maximum values and expand coordinate upper ceiling
star_heights <- df_plot %>%
  group_by(pathway) %>%
  summarize(
    global_max = max(score, na.rm = TRUE),
    global_min = min(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Set star height at 12% above max
    constant_y = global_max + (abs(global_max) * 0.12),
    # Expand upper plot limit to 25% above max to prevent cropping
    y_max_limit = global_max + (abs(global_max) * 0.25)
  )

# Join calculations back to plotting tables
significance_data <- significance_data %>% left_join(star_heights, by = "pathway")
df_plot <- df_plot %>% left_join(star_heights, by = "pathway")

df_plot <- df_plot %>%
  mutate(
    wt_mean_line = case_when(
      pathway == "Proliferation" ~ mean_prolif,
      pathway == "IFN" ~ mean_ifn,
      pathway == "NFkB" ~ mean_nfkb
    )
  )

# 7. Define color palette
pathway_colors <- c(
  "Proliferation" = "black",
  "IFN" = "#e3312a",
  "NFkB" = "#7ac943"
)

# 8. Build Figure 1: Main Pathway Grid
p1 <- ggplot(df_plot, aes(x = position_factor, y = score, fill = pathway)) +
  # Thin violins in background
  geom_violin(
    linewidth = 0.1, 
    alpha = 0.8, 
    scale = "width"
  ) +
  # Little black crossbars to denote the mean score of each individual violin
  geom_errorbar(
    data = all_stats,
    aes(x = position_factor, y = local_mean, ymin = local_mean, ymax = local_mean),
    inherit.aes = FALSE,
    color = "gray50",
    linewidth = 0.2,
    width = 0.6  # Adjust width of the mean bar inside the violin
  ) +
  # Solid black horizontal WT global mean baseline
  geom_line(
    aes(y = wt_mean_line, group = pathway), 
    color = "black", 
    linetype = "solid", 
    linewidth = 0.2
  ) +
  # Vertical thin line segments connecting violin tips to the stars
  geom_segment(
    data = significance_data,
    aes(x = position_factor, xend = position_factor,
        y = max_score+1, yend = constant_y),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.2
  ) +
  # Stars annotation layer
  geom_text(
    data = significance_data,
    aes(x = position_factor, y = constant_y, label = stars),
    inherit.aes = FALSE,
    color = "black",
    size = 3,
    vjust = 0.3
  ) +
  facet_grid(pathway ~ ., scales = "free_y") + # free_y works with custom limits
  geom_blank(aes(y = y_max_limit)) +          # Dynamically forces upper padding safety window
  geom_blank(aes(y = global_min)) +          # Anchors lower range bounds
  scale_fill_manual(values = pathway_colors) +
  scale_x_discrete(drop = FALSE, breaks = NULL) +
  theme_cowplot(font_size = 20) +
  scale_x_discrete(
    drop = FALSE, 
    breaks = c(0, seq(5, 378, by = 5)),
    labels = function(x) ifelse(x == "0", "REF", x)
  ) +
  theme(
    strip.text.y = element_text(angle = -90, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 30),
    strip.background = element_blank(),
    axis.title.y = element_text(size=30),
    
    # Make facet ribbons thinner by minimizing inner element padding margin bounds
    strip.text = element_text(margin = margin(t = 0, b = 0, l = 0, r = 0)), 
    legend.position = "none") +
  labs(y = "Hotspot module score", x=NULL)

# 9. Build Figure 2: Bottom Gray Barplot for Cell Counts (Excluding REF)
df_counts <- df_plot %>%
  filter(position != 0) %>%
  group_by(position_factor, pathway) %>%
  summarize(cell_count = n(), .groups = "drop") %>%
  filter(pathway == "Proliferation") 

p2 <- ggplot(df_counts, aes(x = position_factor, y = cell_count)) +
  geom_bar(stat = "identity", fill = "gray40", width = 0.8) +
  scale_x_discrete(
    drop = FALSE, 
    breaks = c(0, seq(5, 378, by = 5)),
    labels = NULL) +
  theme_cowplot(font_size = 30) +
  theme(
    strip.text.y = element_text(color = "white", angle = 0, hjust = 0),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  scale_y_continuous(
    breaks = c(250,750)) +
  labs(
    x = "Position",
    y = "Cell\ncount")

# Get LOF track
## Make track of LOF marks
# 1. Load your LOF data
lof_data <- read.csv("../cGAMP_sc/LOF_variants.csv")

# Create a master dataframe for all 378 residues
# 1. Prepare data (ensure all 378 residues are present)
df_tracks <- data.frame(residue = 1:378) %>%
  left_join(lof_data, by = c("residue" = "residue")) %>%
  mutate(
    `IFN LOF`  = ifelse(grepl("IFNa", group), "TRUE", "FALSE"),
    `NFkB LOF` = ifelse(grepl("NFkB", group), "TRUE", "FALSE"),
    `UPR LOF`  = ifelse(grepl("UPR", group), "TRUE", "FALSE")
  ) %>%
  mutate(across(c(`IFN LOF`, `NFkB LOF`, `UPR LOF`), ~replace_na(., "FALSE"))) %>%
  pivot_longer(cols = c(`IFN LOF`, `NFkB LOF`, `UPR LOF`), names_to = "Pathway", values_to = "Status")

# Create the top panel containing just the 3 pixel rows
p_tracks <- ggplot(df_tracks, aes(x = residue, y = Pathway,
                                  fill = interaction(Pathway, Status))) +
  geom_tile(color = "white", linewidth = 0) + # Light grid border for crisp pixels
  scale_fill_manual(
    values = c(
      "UPR LOF.TRUE"  = "#35b6e3", "UPR LOF.FALSE"  = "white",
      "IFN LOF.TRUE"  = "#e3312a", "IFN LOF.FALSE"  = "white",
      "NFkB LOF.TRUE" = "#7ac943", "NFkB LOF.FALSE" = "white"
    ),
    guide = "none"
  ) +
  scale_x_continuous(limits = c(0, 378),
                     breaks = seq(5, 375, 5),
                     expand = c(0.0008, 0.0008)) +
  theme_cowplot(font_size = 30) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(), # Hide X text since the bottom plot covers it
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )

# Load mets barplots
p_mets <- readRDS('mets_bars_arcs.rds')

# Define Paths relative to your project workspace
dist_matrix_path <- "hotspot_gene_distances.csv"
modules_path     <- "hotspot_gene_modules.csv"
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
  mutate(Functional_Module = case_when(
    Module_Str == "1" ~ "Proliferation",
    Module_Str == "2"          ~ "IFN",
    Module_Str == "3"          ~ "NFkB"))

# Extract the named vector for mapping
module_vector <- modules_df$Functional_Module
# ==============================================================================

# 3. Clean and Explicit Gene Selection List
target_hallmarks <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
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
               'Mki67','Ccnd1','Ccnf','Ccnt1','Cdc20','Cdc25a','Cdc25b',
               'Cdc27','Cdc45','Cdc6','Cdc7','Cdk1','Cdk4','Cdkn1b','Cdkn2c',
               'Cdkn3',"Myc")

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
  "Proliferation" = "black", 
  "IFN"           = "#e3312a",
  "NFkB"          = "#7ac943" 
)

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
    plot.margin = margin(t = 200, r = 0, b = 0, l = 200, unit = "pt"), # Add 15pt cushion at bottom
    plot.background = element_blank(),
    panel.background = element_blank()
  )

# Load UMAP png
p_umap   <- image_read("hotspot_module_gene_umaps_tumor_4T1.png") 
p_umap   <- wrap_elements(rasterGrob(p_umap)) +
  theme(plot.margin = margin(0,0,0,0,"pt"))

# Load densty plot
p_density <- readRDS('Hotspot_score_2d_density_mets_bonds.rds')


# Combine all sequential rows vertically
# Define the nested layout with a spacer to balance the heatmap legend width
left_column <- ht_grob / (p_umap + plot_spacer() + plot_layout(widths = c(8, 1))) +
  plot_layout(heights = c(3, 1))

p_final <- (left_column | (p_tracks / p1 / p2 / p_mets / p_density + 
                             plot_layout(heights = c(1, 7, 0.75, 4, 4.25))))  &
  theme(
    plot.margin = margin(10, 2, 2, 10, "pt"))

# Coordinates for the left-side Heatmap
xa <- 0.0; ya <- 1  # Panel "a"

# Coordinates for the right-side vertical tracks stack
xb <- 0.0; yb <- 0.29
xc <- 0.415; yc <- 1
xd <- 0.415; yd <- 0.5
xe <- 0.415; ye <- 0.29  

labels_grob_layer <- grobTree(
  textGrob("a", x = unit(xa, "npc"), y = unit(ya, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 40)),
  textGrob("b", x = unit(xb, "npc"), y = unit(yb, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 40)),
  textGrob("c", x = unit(xc, "npc"), y = unit(yc, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 40)),
  textGrob("d", x = unit(xd, "npc"), y = unit(yd, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 40)),
  textGrob("e", x = unit(xe, "npc"), y = unit(ye, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 40)))

p_final_grob <- grobTree(
  patchwork::patchworkGrob(p_final),
  labels_grob_layer
)

# Save the final consolidated multi-panel PDF
pdf('hotspot_module_gene_correlation_sc_score.pdf', width = 40, height = 20)
grid::grid.draw(p_final_grob)
dev.off()

svglite("hotspot_module_gene_correlation_umap.svg", width=18, height=20)
print(left_column)
dev.off()

svglite("hotspot_sc_score_per_position.svg", width=20, height=10)
print(p1 / p2 + plot_layout(heights = c(7, 0.75)))
dev.off()
