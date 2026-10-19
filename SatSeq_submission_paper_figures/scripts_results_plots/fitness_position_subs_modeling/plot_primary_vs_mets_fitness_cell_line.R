rm(list=ls()) # Clear environment

#################
## libraries ####
#################
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(RColorBrewer)
library(ggpmisc)
library(stringr)
library(patchwork)
library(ggrepel)
library(svglite)

# Set working directory to where this script is located
setwd(getSrcDirectory(function(){})) 

### Per position weights---------------------

#################
## load data ####
#################
## 4T1 weights
primary_weights_4T1 <- read.csv('inferred_weights_Initial_vs_Primary_4T1.csv')
colnames(primary_weights_4T1) <- c('position', 'primary_additive_4T1','primary_multiplicative_4T1')

mets_weights_4T1 <- read.csv('inferred_weights_Initial_vs_Met_4T1.csv')
colnames(mets_weights_4T1) <- c('position', 'mets_additive_4T1','mets_multiplicative_4T1')

primary2mets_weights_4T1 <- read.csv('inferred_weights_Primary_vs_Met_4T1.csv')
colnames(primary2mets_weights_4T1) <- c('position', 'primary_mets_additive_4T1','primary_mets_multiplicative_4T1')

tailvein2mets_weights_4T1 <-read.csv('inferred_weights_Tail_Vein_vs_Met_4T1.csv')
colnames(tailvein2mets_weights_4T1) <- c('position', 'tailvein_mets_additive_4T1','tailvein_mets_multiplicative_4T1')

  
## EO771 weight
primary_weights_EO771 <- read.csv('inferred_weights_Initial_vs_Primary_EO771.csv')
colnames(primary_weights_EO771) <- c('position', 'primary_additive_EO771', 'primary_multiplicative_EO771')


#####################
## data processing ##
#####################

# 1. Top Panel Data: Initial vs Mets 4T1 (Mapping as requested)
panel_top <- mets_weights_4T1 %>%
  inner_join(primary_weights_4T1, by = "position") %>%
  select(position, x = primary_additive_4T1, y = mets_additive_4T1) %>%
  mutate(point_label = case_when(
    position == 12  ~ "P13",
    position == 336 ~ "K337",
    TRUE            ~ NA_character_
  ))

# 1. Top Panel Data: Initial vs Mets 4T1 (Mapping as requested)
panel_middle <- primary2mets_weights_4T1 %>%
  inner_join(tailvein2mets_weights_4T1, by = "position") %>%
  select(position, x = tailvein_mets_additive_4T1, y = primary_mets_additive_4T1) %>%
  mutate(point_label = case_when(
    position == 12  ~ "P13",
    position == 336 ~ "K337",
    position == 208  ~ "D209",
    position == 180  ~ "F181",
    position == 197  ~ "Y198",
    position == 106  ~ "F107",
    position == 77  ~ "Q78",
    position == 174  ~ "Q175",
    position == 204  ~ "C205",
    position == 212  ~ "V213",
    TRUE            ~ NA_character_
  ))


# 2. Bottom Panel Data: Primary vs Mets 4T1 (Mapping as requested)
panel_bottom <- primary_weights_4T1 %>%
  inner_join(primary_weights_EO771, by = "position") %>%
  select(position, x = primary_additive_4T1, y = primary_additive_EO771) %>%
  mutate(point_label = case_when(
    position == 12  ~ "P13",
    position == 336 ~ "K337",
    TRUE            ~ NA_character_
  ))


#################
## plotting #####
#################

# Base plot template to keep configurations identical
base_plot <- function(df) {
  ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.75, size = 1.25, color='gray50') +
    geom_smooth(method = "lm", linewidth = 1, color = 'black') +
    geom_text_repel(
      aes(label = point_label),
      color = 'black',
      size = 5,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.5,
      min.segment.length = 0
    ) +
    stat_correlation(
      method = "pearson",
      label.y = "top",
      label.x = 'left',
      mapping = aes(label = paste("atop(",
                                  after_stat(r.label),
                                  ",", after_stat(p.value.label), ")")),
      size = 5,
      parse = TRUE # Evaluates atop() dynamically as stacked text
    )   +
    theme_cowplot() +
    theme(
      # --- Increased Font Sizes ---
      axis.text.x = element_text(size = 16), # Numeric X-tick labels
      axis.text.y = element_text(size = 16), # Numeric Y-tick labels
      axis.title.x = element_text(size = 16), # X-axis main title
      axis.title.y = element_text(size = 16)  # Y-axis main title
    )
}

# Create individual plots with customized Y axis labels
p1 <- base_plot(panel_top) + 
  labs(x = NULL, y = "4T1 Metastatic fitness")

p2 <- base_plot(panel_bottom) + 
  labs(x = "4T1 primary tumor fitness", y = "EO771 primary tumor fitness")

p3 <- base_plot(panel_middle) +
  labs(x = "Tail vein metastatc fitness", y = "Primary metastatic fitness")

# Combine into a single column using patchwork syntax
# p1 goes on top, p2 on the bottom
final_plot <- p1 / p2

svglite("Fitness_primary_mets_cell_lines.svg", width=4, height=6)
print(final_plot)
dev.off()

svglite("Fitness_mets_primary_vs_tail_Vein_4T1.svg", width=6, height=5)
print(p3)
dev.off()


# Per variant weights ---------------
##########################
## processing function ###
##########################

# Standardized function to load, clean, and average per-variant data
load_and_average_fitness <- function(file_path) {
  df <- read_csv(file_path) %>%
    filter(is_control == FALSE) %>%
    mutate(
      wt_residue = str_extract(variant, "^[A-Za-z]+"),
      position = as.numeric(str_extract(variant, "[0-9.]+")),
      mutant_residue = str_extract(variant, "[A-Za-z]+$")
    ) %>%
    filter(
      !is.na(wt_residue), !is.na(position), !is.na(mutant_residue),
      wt_residue %in% c('R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W'),
      mutant_residue %in% c('R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W'),
      wt_residue != mutant_residue
    ) %>%
    # Group by variant to collapse across experiment replicates
    group_by(variant, position) %>%
    summarize(mean_fitness = mean(observed_fitness, na.rm = TRUE),
              .groups = "drop")
  
  return(df)
}

#################
## load data ####
#################

# Load and process the three files per-variant
fit_initial_vs_met_4T1 <- load_and_average_fitness("../../datasets/bulk_fitness/Mammary_Fat_Pad_Initial_vs_Met_boot_4T1.csv")
fit_initial_vs_primary_4T1 <- load_and_average_fitness("../../datasets/bulk_fitness/Mammary_Fat_Pad_Initial_vs_Primary_boot_4T1.csv")
fit_initial_vs_primary_EO771 <- load_and_average_fitness("../../datasets/bulk_fitness/Mammary_Fat_Pad_Initial_vs_Primary_boot_EO771.csv")

#####################
## data processing ##
#####################

# 1. Top Panel Data: Comparison of 4T1 variants between Primary and Met selections
panel_top <- fit_initial_vs_primary_4T1 %>%
  inner_join(fit_initial_vs_met_4T1, by = c("variant", "position"),
             suffix = c("_primary", "_met")) %>%
  select(position, variant, x = mean_fitness_primary, y = mean_fitness_met) %>%
  mutate(point_label = case_when(
    position == 13  ~ "P13",
    position == 337 ~ "K337",
    TRUE            ~ NA_character_
  ))

# 2. Bottom Panel Data: Comparison of Initial vs Primary variants between 4T1 and EO771
panel_bottom <- fit_initial_vs_primary_4T1 %>%
  inner_join(fit_initial_vs_primary_EO771, by = c("variant", "position"), suffix = c("_4T1", "_EO771")) %>%
  select(position, variant, x = mean_fitness_4T1, y = mean_fitness_EO771) %>%
  mutate(point_label = case_when(
    position == 13  ~ "P13",
    position == 337 ~ "K337",
    TRUE            ~ NA_character_
  ))


#################
## plotting #####
#################

# Base plot template to keep 2D density and directional statistics layouts identical
base_plot_2d <- function(df, low_color, high_color) {
  
  # Filter out data specifically for P13 (position 12) and K337 (position 336)
  highlight_points <- df %>%
    filter(position %in% c(13, 337)) %>%
    mutate(mutant_id = case_when(
      position == 13  ~ "P13",
      position == 337 ~ "K337"
    ))
  
  ggplot(df, aes(x = x, y = y)) +
    # 1. Filled 2D Density Contours (Background)
    geom_density_2d_filled(bins = 50, linewidth = 0, alpha=1,
                           show.legend = FALSE) +
    scale_fill_manual(values = colorRampPalette(c(low_color,
                                                  high_color))(50)) +
    # 2. Linear fit line
    geom_smooth(method = "lm", linewidth = 0.5,
                color = 'black',
                se = FALSE) +
    
    # 3. Layer individual highlighted points explicitly over the density contours
    geom_point(
      data = highlight_points,
      aes(x = x, y = y, color = mutant_id),
      size = 1,
      shape = 16
    ) +
    
    # 4. Explicit manual hex styling for the target locations
    scale_color_manual(
      values = c("P13" = "#d01c8b", "K337" = "#5e3c99")
    ) +
    
    # 5. Multi-line Correlation showing BOTH signed R and P value
    stat_correlation(
      method = "pearson",
      label.y = "top",
      label.x = 'left',
      mapping = aes(label = paste("atop(",
                                  after_stat(r.label),
                                  ",", after_stat(p.value.label), ")")),
      size = 4,
      parse = TRUE # Evaluates atop() dynamically as stacked text
    )  +
    theme_cowplot() +
    guides(fill = "none") + # Turn off density fill legend
    theme(legend.title = element_blank()) # Clean title off the point legend
}

# Create individual plots with customized color palettes and Y axis labels
p1 <- base_plot_2d(panel_top,"white", "gray10") + 
  labs(x = NULL, y = "4T1 Metastatic fitness")

p2 <- base_plot_2d(panel_bottom, "white", "gray10") + 
  labs(x = "4T1 primary tumor fitness", y = "EO771 primary tumor fitness")

# Combine into a single column using patchwork syntax
final_plot <- p1 / p2

# Save plot
svglite("Fitness_primary_mets_cell_lines_per_variant.svg", width=5, height=6)
print(final_plot)
dev.off()




