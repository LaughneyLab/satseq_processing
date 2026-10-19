rm(list=ls())#Clear environment
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# 1. Load required libraries
library(readxl)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(svglite)
library(dplyr)
library(readr)
library(scales)
library(ggplot2)
library(ggrepel)

# 2. Define parameters and file paths
excel_file   <- "Dual3_LOF_flow_forAmitava.xlsx"
csv_file     <- "LOF_variants.csv"
sheet_ifn    <- "IFN_log2(FC)"
sheet_nfkb   <- "NFkB_log2(FC)"

# ==========================================
# 3. Read and Process the LOF Variant List
# ==========================================
# Read CSV file and handle potential header spaces/case-sensitivity
lof_meta <- read_csv(csv_file) %>%
  rename_all(tolower) %>%
  mutate(
    variant = str_trim(as.character(variant)),
    group   = str_trim(as.character(group))
  )

# Build logic lists: Match text strings dynamically & inject the 'STOP' controls
ifn_lof_list <- lof_meta %>% 
  filter(str_detect(group, "IFNa")) %>% 
  pull(variant) %>% 
  c("STOP","dCTT","RRDD") %>% 
  unique()

nfkb_lof_list <- lof_meta %>% 
  filter(str_detect(group, "NFkB")) %>% 
  pull(variant) %>% 
  c("STOP","dCTT","RRDD") %>% 
  unique()

# ==========================================
# 4. Helper Function to Process Each Sheet
# ==========================================
process_screen_sheet <- function(file_path, sheet_name, target_lof_list, pathway_label) {
  
  # Read Excel sheet
  raw_df <- read_excel(file_path, sheet = sheet_name)
  
  # Ensure the first column has a standard tracking label
  colnames(raw_df)[1] <- "residue"
  replicate_cols <- colnames(raw_df)[-1]
  
  # Clean residue string elements
  processed_df <- raw_df %>%
    filter(str_ends(residue, "X")|residue %in% c("REF", "STOP", "dCTT","RRDD")) %>%
    mutate(residue = str_trim(as.character(residue)),
           residue = case_when(
             residue %in% c("REF", "STOP", "dCTT","RRDD") ~ residue,
             str_ends(residue, "X")        ~ str_sub(residue, 1, -2),
             TRUE                          ~ residue
           ))
  
  # Extract the reference baseline adjustments
  ref_row <- processed_df %>% filter(residue == "REF")
  if (nrow(ref_row) == 0) {
    stop(paste("Error: 'REF' identifier missing in sheet:", sheet_name))
  }
  ref_values <- as.numeric(ref_row[, replicate_cols])
  
  # Center and baseline the raw values via vector broadcasting
  normalized_df <- processed_df %>%
    mutate(across(all_of(replicate_cols), ~ . - ref_values[match(cur_column(), replicate_cols)])) %>%
    filter(residue != "REF")
  
  # Condense downstream metrics down to single-row residue averages
  final_df <- normalized_df %>%
    mutate(
      mean_normalized_log2FC = rowMeans(select(., all_of(replicate_cols)), na.rm = TRUE),
      Pathway = pathway_label,
      Group = if_else(residue %in% target_lof_list, "LOF\nvariants",
                      "Other\nvariants"),
      Group = factor(Group, levels = c("LOF\nvariants","Other\nvariants"))
    ) %>%
    filter(!is.na(mean_normalized_log2FC)) %>%
    select(residue, mean_normalized_log2FC, Pathway, Group) %>%
    filter(!residue %in% c("P13A","K337R"))
  
  # Remove individually tested residues P13A and K337R
  return(final_df)
}

# ==========================================
# 5. Parse Data Sheets and Merge Results
# ==========================================
ifn_clean  <- process_screen_sheet(excel_file, sheet_ifn, ifn_lof_list, "IFN")
nfkb_clean <- process_screen_sheet(excel_file, sheet_nfkb, nfkb_lof_list, "NFkB")

# Read sc scores
scores <- read_csv("sc_MAST_pathway_score_4T1.csv") 

# Combine both parsed screens into a single long-format pipeline data frame
combined_data <- bind_rows(ifn_clean, nfkb_clean) %>%
  mutate(Pathway = factor(Pathway, levels = c("IFN", "NFkB"))) %>%
# Create mapping helper variables to cleanly control fills vs borders
  mutate(
    # LOF gets the pathway color code; LOF is forced to "Hollow" (White)
    fill_color = if_else(Group == "Other\nvariants", as.character(Pathway),
                         "Hollow"),
    # Both groups pass their literal Pathway string name down to map line colors
    border_color = as.character(Pathway)
  ) %>%
  # Add sc scores
  left_join(scores, by = c("residue" = "variant")) %>%
  mutate(
    sc_score = case_when(
      Pathway == "IFN"  ~ ifna_response,
      Pathway == "NFkB" ~ nfkb_response,
      TRUE              ~ NA_real_
    )
  ) %>%
  # 3. Clean up the intermediate merge columns
  select(-er_stress, -ifna_response, -nfkb_response) %>%
  mutate(
    #Replace NA with 0 so it takes the smallest "base" size on the plot
    sc_score = if_else(is.na(sc_score), 0, sc_score)
  ) %>%
  mutate(label = if_else(sc_score< -5, residue, NA))

# ==========================================
# 6. Generate the Faceted Violin Plot via ggpubr
# ==========================================
custom_fill_palette  <- c("IFN" = "gray", "NFkB" = "gray",
                          "Hollow" = "#FFFFFF")
custom_color_palette <- c("IFN" = "gray", "NFkB" = "gray")


plot_output <- ggplot(combined_data, aes(x = Group, y = mean_normalized_log2FC)) +
  # Violin Layer: We explicitly control the border width (linewidth) using a conditional vector.
  # LOF gets a crisp 0.8 border line; LOF gets a 0 line width (no border).
  geom_violin(
    aes(fill = fill_color, color = border_color, linewidth = Group),
    trim = FALSE, alpha = 0.85
  ) +
  # Map line widths explicitly based on the Group factor
  scale_linewidth_manual(values = c("Other\nvariants" = 0,
                                    "LOF\nvariants" = 0.8),
                         guide = "none") +
  
  # Jitter Data Points Layer
  geom_jitter(
    aes(size = scales::squish(sc_score, c(-5, 0))),
    shape = 16,
    color = "black", 
    width = 0.06, alpha = 0.5, show.legend = T
  ) +
  stat_summary(
    fun = mean, 
    geom = "crossbar", 
    width = 0.25, 
    color = "black", 
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  
  # Apply Custom Color Palettes
  scale_fill_manual(values = custom_fill_palette, guide = "none") +
  scale_color_manual(values = custom_color_palette, guide = "none") +
  scale_size_continuous(
    name = "Single cell score",
    limits = c(-5, 0),
    range = c(3, 1)) +
  
  # Apply facet separation grid using the Pathway column factor
  facet_wrap(~ Pathway, scales = "free_y") +
  # Calculate and print independent Wilcoxon significance tags per window box
  stat_compare_means(comparisons = list( c("LOF\nvariants",
                                           "Other\nvariants") ),
                     method = "wilcox.test",
                     method.args = list(alternative = "less"),
                     label = "p.format") +
  
  # Format layout styles to match standard journal publication defaults
  theme_cowplot() +
  labs(y="logFC change from REF") +
  theme(
    strip.text = element_text(size = 12), # Styling facet titles
    strip.background = element_rect(fill = "grey85"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(),
    axis.title.x = element_blank()
  )

svglite("FACS_LOF_validation_Dual_Reporter.svg", width = 7, height = 3.5)
print(plot_output)
dev.off()
