rm(list=ls()) # Clear environment

library(ggplot2)
library(scales)
library(dplyr)
library(tidyverse)
library(cowplot)
library(patchwork)
library(readr)
library(readxl)
library(ggrepel)

# Set working directory to where this script is
setwd(getSrcDirectory(function(){})) 

# 1. Load fitness and sd data
variants <- read.csv("FACS_mets_weights_all.csv")[,2]
fitness_mean <- read.csv("../fitness_position_subs_modeling/inferred_weights_Primary_vs_Met_4T1.csv")[,2]
fitness_sd <- read.csv("../fitness_position_subs_modeling/inferred_weights_sd_Primary_vs_Met_4T1.csv")[,2]

# Full variant df with both position and variant info
variants_full <- read.csv("FACS_mets_weights_all.csv")[,c(1,2)]

# Load standard error data, convert type, and average duplicates
sderr_data <- read_xlsx("4T1_MFP_Primary_vs_Met_Barplot.xlsx", sheet = "Fully Aggregated Data") %>% 
  select(position, standard_error) %>%
  mutate(position = as.integer(position)) %>%
  group_by(position) %>%
  summarise(standard_error = mean(standard_error, na.rm = TRUE), .groups = "drop")

# Merge standard_error by matching the position column
fitness_df <- data.frame(variants, fitness_mean, fitness_sd) %>%
  left_join(variants_full, by = c("variants" = "variant")) %>% 
  left_join(sderr_data, by = "position")

# 2. Calculate the standard deviation OF the fitness_mean values
mets_cutoff <- sd(fitness_df$fitness_mean)

# 3. Subset data for hi mets to add labels
labeled_variants <- subset(fitness_df, fitness_mean > mets_cutoff)

# 4. Create the plot
p <- ggplot(fitness_df, aes(x = fitness_mean, y = fitness_sd)) +
  geom_point(aes(size = standard_error, color = standard_error), alpha = 0.6) +
  scale_color_viridis_c() + 
  # Larger size for smaller error (reverses the scale)
  scale_size_continuous(trans = "reverse", range = c(0.5, 3)) + 
  geom_vline(xintercept = mets_cutoff, linetype = "dashed",
             color = "gray40", linewidth = 1) +
  geom_text_repel(data = labeled_variants, aes(label = variants), max.overlaps = 35) +
  theme_cowplot() +
  labs(
    x = "Fitness posterior mean", 
    y = "Fitness posterior SD",
    color = "Standard Error",
    size = "Standard Error"
  )

ggsave(plot=p, 'Mets_mean_sd_sderr.pdf',
       width=8, height=6)