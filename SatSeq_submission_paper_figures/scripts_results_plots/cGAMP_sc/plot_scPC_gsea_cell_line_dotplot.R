rm(list=ls())#Clear environment

#################
## libraries ####
#################
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(readxl)
library(stringr)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggpmisc)
library(tidyverse)
library(reshape2)
library(svglite)



setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is



## Plot PC hallmark NES------
# Load Variance Data
var_df <- read.csv("PC_variance_explained.csv")
colnames(var_df)[1] <- "PC"

# Common plot function
process_gsea_results <- function(file_path, cell_line) {
  df <- read.csv(file_path)
  
  # 1. Filter: Identify terms that meet the threshold
  keep_terms <- df %>%
    filter(NOM.p.val < 0.05 & FDR.q.val < 0.15) %>%
    pull(Term) %>%
    unique()
  
  df_filtered <- df %>% filter(Term %in% keep_terms)
  
  # 2. Reshape to Matrix for Clustering
  mat_nes <- df_filtered %>%
    select(Term, Name, NES) %>%
    pivot_wider(names_from = Name, values_from = NES, values_fill = 0) %>%
    column_to_rownames("Term")
  
  # 3. Perform Hierarchical Clustering
  hr <- hclust(dist(mat_nes), method = "complete")
  pathway_order <- rownames(mat_nes)[hr$order]
  
  # 4. Clean Names and Factor Levels
  df_plot <- df_filtered %>%
    mutate(
      Term_Clean = Term %>%
        str_remove_all("^HALLMARK_|^GOBP_") %>%
        str_replace_all("_", " ") %>%
        str_to_title() %>%
        str_replace("Tnf", "TNF") %>%
        str_replace("Il2 Stat5", "IL2 STAT5") %>%
        str_replace("Dna", "DNA") %>%
        str_replace("Tgf Beta", "TGFb") %>%
        str_replace("E2f", "E2F") %>%
        str_replace("G2m", "G2M") %>%
        str_replace("Uv", "UV") %>%
        str_replace("Nfkb", "NFkB") %>%
        #str_replace("Epithelial Mesenchymal Transition", "EMT") %>%
        #str_replace("Macroautophagy", "Autophagy") %>%
        str_replace("Il6 Jak Stat3", "IL6 JAK STAT3") %>%
        str_replace("Reactive Oxygen Species", "ROS"),
      Term = factor(Term, levels = pathway_order),
      Name = factor(Name, levels = paste0("PC", 1:10))
    )
  
  # --- 5. Create the Variance Barplot (p_top) ---
  # Extract relevant variance and convert to percentage
  line_var <- var_df %>%
    select(PC, !!sym(cell_line)) %>%
    filter(PC %in% paste0("PC", 1:10)) %>%
    mutate(
      PC = factor(PC, levels = paste0("PC", 1:10)),
      VarPercent = !!sym(cell_line) * 100
    )
  
  # 5. Generate the Plot
  p_main <- ggplot(df_plot, aes(x = Name, y = Term)) +
    geom_point(aes(size = -log10(NOM.p.val + 1e-10), color = NES),
               show.legend = !(cell_line=="EO771")) + # Added small epsilon to avoid log(0)
    scale_color_gradientn(colors = rev(brewer.pal(7, "RdBu")), limits = c(-2, 2)) +
    scale_y_discrete(labels = setNames(df_plot$Term_Clean, df_plot$Term)) +
    scale_size_continuous(range = c(0.5, 4), breaks = c(2, 4, 6, 8)) + 
    theme_minimal() +
    labs(
      x = NULL,
      y = paste0(sub("_.*", "", basename(file_path)), " pathway"),
      size = expression(-log[10](P)),
      color = "Enrichment"
    ) +
    theme(
      axis.text.x = element_text(size = 10,angle = 45, hjust = 1),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      strip.background = element_rect(fill = "grey90"),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    )
  
  p_top <- ggplot(line_var, aes(x = PC, y = VarPercent)) +
    geom_bar(stat = "identity", fill = "black", width = 0.7) +
    theme_minimal() +
    scale_y_continuous(expand = c(0, 0)) +
    labs(y =NULL, x = NULL) +
    labs(title = "% Variance explained", x = NULL, y = NULL) +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5, face = "plain",
                                margin = margin(b = 5)),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(), # Removed all grid lines
      # Moves title very close to the axis line
      axis.text.y = element_text(size = 10),
      plot.margin = margin(b = -10, l = -5, r = 5, t = 5) # Tighten bottom margin to connect to heatmap
    )
  
  # --- 6. Combine using Patchwork ---
  # Align the plots vertically
  p_combined <- p_top / p_main + plot_layout(heights = c(1, 6))
  return(p_combined)
}

# Save plot
p8 <- process_gsea_results("4T1_PC_cGAMP_gsea.csv", "X4T1")
p9 <- process_gsea_results("EO771_PC_cGAMP_gsea.csv", "EO771")

# Save plot
 # ggsave(plot=p10, 'scPC_pathways_EO771_4T1.pdf',
 #        width=5, height=9)



svglite("PC1-10_pathways_heatmap_4T1.svg", width = 6, height = 4.5)
print(p8)
dev.off()

saveRDS(p8, file = "PC1-10_pathways_heatmap_4T1.rds")

svglite("PC1-10_pathways_heatmap_EO771.svg", width = 6, height = 4.5)
print(p9)
dev.off()

saveRDS(p9, file = "PC1-10_pathways_heatmap_EO771.rds")
