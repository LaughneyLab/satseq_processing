rm(list=ls())#Clear environment
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(patchwork)
library(svglite)


plot_gene_cpm <- function(tsv_path, target_gene,
                          mapping_dir = "7J9GL7_mapping-stats") {
  
  # 1. Dynamically read and parse filenames from the mapping directory
  if (!dir.exists(mapping_dir)) {
    stop(paste("Directory", mapping_dir, "does not exist."))
  }
  
  files <- list.files(mapping_dir, pattern = "\\.tsv$", full.names = FALSE)
  
  if (length(files) == 0) {
    stop(paste("No .tsv files found in", mapping_dir))
  }
  
  # Build metadata from file patterns like "7J9GL7_4_4_K337R_0.tsv"
  sample_metadata <- data.frame(filename = files, stringsAsFactors = FALSE) %>%
    mutate(
      sample_id = str_extract(filename, "^7J9GL7_\\d+"),
      Condition = str_split_i(filename, "_", 4),
      Time_raw = str_split_i(str_remove(filename, "\\.tsv$"), "_", 5)
    ) %>%
    mutate(
      Time = ifelse(Time_raw == "chronic", "Chronic", paste0(Time_raw, "h")),
      sample_id_R = paste0("X", sample_id)
    ) %>%
    select(sample_id_R, Condition, Time)
  
  # Enforce chronological ordering on X-axis
  time_order <- c("0h", "2h", "6h", "12h", "Chronic")
  actual_times <- intersect(time_order, unique(sample_metadata$Time))
  sample_metadata$Time <- factor(sample_metadata$Time, levels = actual_times)
  
  # Enforce strict group ordering and link to custom colors
  custom_colors <- c(
    "REF"   = "#4D4D4D", 
    "STOP"  = "#B3B3B3", 
    "P13A"  = "#FF5210", 
    "K337R" = "#FF1D25"
  )
  
  sample_metadata$Condition <- factor(
    sample_metadata$Condition, 
    levels = names(custom_colors)
  )
  
  # 2. Read and filter the RNA-seq data matrix
  df <- read.delim(tsv_path, sep = "\t", stringsAsFactors = FALSE, check.names = TRUE)
  
  # Filter for the target gene case-insensitively
  gene_df <- df %>% 
    filter(tolower(gene_name) == tolower(target_gene))
  
  if (nrow(gene_df) == 0) {
    stop(paste("Gene", target_gene, "not found in the dataset."))
  }
  
  # 3. Shape data long and map data structures together
  plot_data <- gene_df %>%
    select(gene_name, ends_with("_cpm")) %>%
    pivot_longer(
      cols = ends_with("_cpm"),
      names_to = "raw_colname",
      values_to = "CPM"
    ) %>%
    mutate(CPM = log10(1+CPM)) %>%
    mutate(sample_id_R = str_remove(raw_colname, "_cpm")) %>%
    inner_join(sample_metadata, by = "sample_id_R")
  
  if (nrow(plot_data) == 0) {
    stop("Could not map any active columns in TSV to the file profiles found in mapping-stats.")
  }
  
  # 4. Generate the bar plot with custom styling
  p <- ggplot(plot_data, aes(x = Time, y = CPM, fill = Condition)) +
    geom_bar(
      stat = "identity", 
      position = position_dodge(width = 0.8), 
      width = 0.7
    ) +
    # Force custom explicit colors
    scale_fill_manual(values = custom_colors) +
    theme_cowplot() +
    labs(
      title = target_gene,
      x = NULL,
      y = "Log CPM",
      fill = NULL
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 60, hjust = 1))
  
  return(p)
}


# Plot per gene
p1 <- plot_gene_cpm("7J9GL7_results/7J9GL7-expression-matrix.tsv",
              "Cxcl10", "7J9GL7_results/7J9GL7_mapping-stats")

p2 <- plot_gene_cpm("7J9GL7_results/7J9GL7-expression-matrix.tsv",
                    "Isg15", "7J9GL7_results/7J9GL7_mapping-stats")

p3 <- plot_gene_cpm("7J9GL7_results/7J9GL7-expression-matrix.tsv",
                    "Ifit3", "7J9GL7_results/7J9GL7_mapping-stats")

p4 <- plot_gene_cpm("7J9GL7_results/7J9GL7-expression-matrix.tsv",
                    "Ifnb1", "7J9GL7_results/7J9GL7_mapping-stats")

p <- (p1 | p4) / (p2 | p3)

svglite("bulk_RNAseq_P13_K337.svg", width=8, height=5)
print(p)
dev.off()
