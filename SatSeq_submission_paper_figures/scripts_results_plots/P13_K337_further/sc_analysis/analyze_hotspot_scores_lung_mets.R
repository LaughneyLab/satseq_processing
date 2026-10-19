rm(list=ls())#Clear environment
# Set working directory to where this script is
setwd(getSrcDirectory(function(){})[1]) 

# 1. Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggsignif)
library(ggpubr)
library(scales)
library(patchwork)
library(tidyr)
library(readxl)
library(stringr)
library(svglite)
library(stringr)
library(magick)
library(grid)
library(RColorBrewer)

# Do module gsea dotplot

# --- 1. Define Variables ---
file_path <- "hotspot_modules_gsea_lung_mets.xlsx"

# --- 2. Load Data from Excel Sheets ---
df_pval <- read_excel(file_path, sheet = "Adjusted P-value")
df_over <- read_excel(file_path, sheet = "Overlap")

# --- 3. Reshape and Filter P-values ---
df_pval_long <- df_pval %>%
  pivot_longer(cols = -Term, names_to = "Module", values_to = "P_val") %>%
  mutate(Module = factor(Module, levels = unique(Module)))

# Identify terms where at least one module meets the significance threshold
keep_terms <- df_pval_long %>%
  filter(P_val < 0.05) %>%
  pull(Term) %>%
  unique()

# --- 4. Reshape Overlap Matrix and Calculate Fractions ---
df_over_long <- df_over %>%
  pivot_longer(cols = -Term, names_to = "Module", values_to = "Overlap_Str") %>%
  separate(Overlap_Str, into = c("Overlap_Genes", "Total_Genes"), sep = "/", convert = TRUE) %>%
  mutate(
    Fraction = ifelse(Total_Genes > 0, Overlap_Genes / Total_Genes, 0)
  )

# --- 5. Merge Datasets and Apply Filter ---
df_combined <- df_pval_long %>%
  inner_join(df_over_long, by = c("Term", "Module")) %>%
  filter(Term %in% keep_terms)

# --- 6. Perform Hierarchical Clustering for Term Order ---
mat_clustering <- df_combined %>%
  mutate(NegLogP = -log10(P_val)) %>% 
  select(Term, Module, NegLogP) %>%
  pivot_wider(names_from = Module, values_from = NegLogP, values_fill = 0) %>%
  column_to_rownames("Term")

hr <- hclust(dist(mat_clustering), method = "complete")
pathway_order <- rownames(mat_clustering)[hr$order]

# --- 7. Clean Text Names and Apply Factor Order ---
df_plot <- df_combined %>%
  mutate(
    Term_Clean = Term %>% 
      str_remove_all("^HALLMARK_|^GOBP_") %>% 
      str_replace_all("_", " ") %>% 
      str_to_title() %>% 
      str_replace("Tnfa Signaling Via NfkB", "TNFa Signaling via NFkB") %>% 
      str_replace("Il2 Stat5", "IL2 STAT5") %>% 
      str_replace("Dna", "DNA") %>% 
      str_replace("Tgf Beta", "TGFb") %>% 
      str_replace("E2f", "E2F") %>% 
      str_replace("G2m", "G2M") %>% 
      str_replace("Uv", "UV") %>% 
      str_replace("NfKb", "NFkB") %>% 
      str_replace("Il6 Jak Stat3", "IL6 JAK STAT3") %>% 
      str_replace("Reactive Oxygen Species", "ROS"),
    Term = factor(Term, levels = pathway_order)
  )


df_plot$Module <- paste0("Module ", df_plot$Module)

p_gsea <- ggplot(df_plot, aes(x = -log10(P_val), y = Term)) +
  # Create solid black bars with uniform heights
  geom_col(fill = "black", color = "white", linewidth = 0.2, width = 0.8) + 
  # Apply the cleaned up labels directly to the shared y-axis
  scale_y_discrete(labels = setNames(df_plot$Term_Clean, df_plot$Term)) +
  # Force layout into a single row with a globally shared y-axis
  facet_wrap(~Module, nrow = 1, scales = "fixed") + 
  theme_cowplot() +
  geom_vline(xintercept = -log10(0.05),linetype="dashed",color='gray') +
  
  labs(
    x = "-log10(P value)", 
    y = "Hallmark Pathway"
  ) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    # Style the module facet banner labels
    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "gray85", colour = "black", linewidth = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10))
  )

svglite("Hotspot_module_enrichR_barplot.svg", width = 6, height = 4)
print(p_gsea)
dev.off()



# 3. Load hotspot data & assign structural groups
hotspot_df <- read.csv('hotspot_module_scores_4T1_lung_mets.csv',row.names = 1) %>%
  mutate(variant = as.character(sample) %>% str_split_i("_", 1)) %>%
  pivot_longer(cols = starts_with("X"), names_to = "module", values_to = "score") %>%
  mutate(
    # Extract the number from 'X1', 'X2', etc.
    module_num = as.numeric(str_extract(module, "\\d+")),
    # Create the new clean name
    module_clean = paste0("Module ", module_num),
    # Reorder the factor levels based on the numeric order
    module_clean = factor(module_clean, levels = paste0("Module ", 1:5))
  )  %>% 
  filter(variant %in% c("REF","STOP", "P13","K337")) %>%
  mutate(
    # Explicitly set the x-axis display order for your plot
    variant = fct_relevel(variant, "REF","STOP", "P13","K337"))

# 4. Define specific Wilcoxon comparisons
my_comparisons <- list(
  c("REF","STOP"),c("REF","P13"),c("REF", "K337"), c("P13","K337"))

# 5. Plot the Violin/Boxplot with Significance
p_score <- ggplot(hotspot_df, aes(x = variant, y = score, color = variant)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) + 
  #stat_summary(fun = mean, geom = "crossbar", width = 0.3,
  #color = "black", linewidth = 0.3, show.legend = FALSE) +
  geom_hline(yintercept = 0, color="black", linetype="dashed", linewidth=0.4) +
  facet_wrap(~module_clean, scales = "free_y", nrow = 1) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format", 
    bracket.size = 0.5,
    size = 4) +
  scale_color_manual(
    values = c("P13" = "#FBB040",
               "K337" = "#E43B3F",
               "STOP" = "#B4B4B4",
               "REF" = "#4E4E4E")
  ) +
  theme_cowplot() +
  labs(x = NULL, y = "Module Score", fill = "Variant") +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "none",
    strip.text = element_text(margin = margin(t = 1, b = 1))
  )

# Save plot
svglite("Hotspot_score_lung_mets.svg", width = 8, height = 4)
print(p_score)
dev.off()

# Load umaps
# Import 3d plots of weights
image_add <- image_read('hotspot_module_score_umaps_lung_mets.png')
image_add <- wrap_elements(rasterGrob(image_add)) + 
  theme(plot.margin = margin(0, 0, 0, 0, "pt"))

p_all <- ((p_gsea | p_score) + plot_layout(widths = c(1, 4))) / image_add
svglite("Hotspot_lung_mets.svg", width = 20, height = 12)
print(p_all)
dev.off()
# ==========================================
# Plot sc 2d density contours by Variant
# ==========================================

# Step A: Pair IFN and NFkB scores for each individual cell using variant groupings
plot_wide <- hotspot_df %>% 
  filter(module_clean %in% c("NFkB", "IFN")) %>% 
  # Group by variables present in the data frame to index cells uniquely per variant
  group_by(sample, variant, module_clean) %>% 
  mutate(cell_index = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(
    id_cols = c(sample, variant, cell_index), 
    names_from = module_clean, 
    values_from = score
  ) %>% 
  drop_na(NFkB, IFN)

# Define color mapping exactly matching your original boxplot colors
variant_colors <- c(
  "REF"  = "#4E4E4E",
  "STOP" = "#B4B4B4",
  "P13"  = "#FBB040", 
  "K337" = "#E43B3F"
)

p4 <- ggplot(plot_wide, aes(x = IFN, y = NFkB)) +
  # 1. Filled 2D density contour bands (Gradients from white to variant color)
  geom_density_2d_filled(
    aes(fill = variant, alpha = after_stat(level)), # Map fill to variant, alpha to density depth
    bins = 8, 
    show.legend = FALSE
  ) +
  # Control the alpha range so the lowest density level is white/transparent
  scale_alpha_manual(values = seq(0, 0.85, length.out = 8)) +
  # 2. Extreme outer contour lines matching the variant color
  geom_density_2d(
    aes(color = variant), 
    bins = 8, 
    linewidth = 0.5, 
    show.legend = FALSE
  ) +
  # Intersecting baseline grids at 0
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
  # Separate the plot into 4 distinct panels based on your variant factor levels
  facet_wrap(~variant, ncol = 4) +
  # Apply your specific variant colors to fill and lines
  scale_fill_manual(values = variant_colors) +
  scale_color_manual(values = variant_colors) +
  # Labels and themes
  labs(
    x = "IFN Module score", 
    y = "NFkB Module score"
  ) +
  theme_cowplot(font_size = 20) +
  # Use coord_cartesian to prevent data clipping near borders
  coord_cartesian(xlim = c(-20, 15), ylim = c(-15, 20)) +
  theme(
    legend.position = "none", 
    strip.background = element_blank(), 
    strip.text = element_text(size = 20), 
    axis.line = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75)
  )

# Save plot
svglite("Hotspot_score_2d_density_lung_mets.svg", width = 10, height = 6)
print(p4)
dev.off()
