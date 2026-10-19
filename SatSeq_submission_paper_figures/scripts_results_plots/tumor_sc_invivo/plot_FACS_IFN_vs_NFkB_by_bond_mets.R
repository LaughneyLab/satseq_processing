rm(list=ls())#Clear environment

library(stringr)
library(ggplot2)
library(scales)
library(dplyr)
library(tidyverse)
library(cowplot)
library(patchwork)
library(readr)
library(ggpubr)
library(rstatix)
library(readxl)
library(ggrepel)
library(svglite)
library(ggnewscale)


setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# 1. Load contact matrices
bonds_AA  <- read.csv('bonds_monomer.csv', row.names = 1)
bonds_AB  <- read.csv('bonds_dimer.csv', row.names = 1)
bonds_ACD <- read.csv('bonds_tetramer.csv', row.names = 1)

# 2. Load fitness and mets data
fitness_data <- read.csv("../cGAMP_sc/FACS_weights_all.csv")[1:378,]
mets_data <- read.csv('../fitness_position_subs_modeling/inferred_weights_Primary_vs_Met_4T1.csv')
seq_conservation <- read.csv('conservation.csv')

fitness_data$mets_add <- mets_data$additive
fitness_data$position <- 1:378

# 3. Build the wide dataframe using fast column-summing
mets_bonds_wide <- data.frame(
  position = 1:378,
  monomer  = colSums(bonds_AA[, 1:378, drop = FALSE], na.rm = TRUE),
  dimer    = colSums(bonds_AB[, 1:378, drop = FALSE], na.rm = TRUE),
  tetramer = colSums(bonds_ACD[, 1:378, drop = FALSE], na.rm = TRUE)
) %>% 
  # 4. Bring in the fitness column directly
  left_join(select(fitness_data, position, variant, mets_add,
                   STING_expression, Dual_IFN_expression, NFkB_expression,
                   ATF3_expression),
            by = "position")


# 1. Calculate fractions and categorize groups
plot_data <- mets_bonds_wide %>%
  mutate(
    total_bonds = monomer + dimer + tetramer,
    monomer_fraction = ifelse(total_bonds > 0, monomer / total_bonds, 0),
    dimer_fraction   = ifelse(total_bonds > 0, dimer / total_bonds, 0),
    tetramer_fraction   = ifelse(total_bonds > 0, tetramer / total_bonds, 0),
    
    # Classify based on your fractional rules
    bond_group = case_when(
      monomer_fraction > dimer_fraction ~ "Intra",
      dimer_fraction >= monomer_fraction ~ "Inter",
      TRUE                   ~ NA_character_ 
    )) %>%
      filter(total_bonds > 0, !is.na(bond_group)) %>%
    mutate(bond_group = factor(bond_group,
                               levels = c("Intra", "Inter"))) 

plot_data$label <- ifelse(plot_data$mets_add > sd(plot_data$mets_add),
                          plot_data$variant, NA)


p <- ggplot(plot_data, aes(x = Dual_IFN_expression, y = NFkB_expression)) +
  # 1. INTRA CONTOURS (White to #5e3c99)
  geom_density_2d_filled(
    data = ~ filter(.x, bond_group == "Intra"),
    bins = 5, linewidth = 0, alpha=1, show.legend = FALSE) +
  
  scale_fill_manual(values = colorRampPalette(c("white", "#5e3c99"))(5)) +
  
  # Clear the fill scale to start a new one for Inter
  new_scale_fill() +
  
  # 2. INTER CONTOURS (White to #d01c8b)
  geom_density_2d_filled(
    data = ~ filter(.x, bond_group == "Inter"),
    bins = 5, alpha = 0.25, linewidth = 0, show.legend = FALSE) +
  
  scale_fill_manual(values = colorRampPalette(c("white", "#d01c8b"))(5)) + 
  
  new_scale_fill() +
  
  geom_point(
    aes(fill = bond_group, size = scales::squish(mets_add, c(0, 0.5))), 
    shape = 21,       # Filled circle shape that allows black borders
    color = "black",  # Border color
    stroke = 0.5,     # Border line thickness
    alpha = 0.8) +
  geom_hline(yintercept = 0, color = "black") + 
  geom_vline(xintercept = 0, color = "black") + 
  
  scale_x_continuous(trans = pseudo_log_trans(sigma = 0.25, base = 10),
                     breaks = scales::breaks_pretty(n = 8)) +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 0.25, base = 10),
                     breaks = scales::breaks_pretty(n = 8)) +
  # Repel text labels safely away from points
  geom_text_repel(
    aes(label = label),
    size = 3.5,                 # Text font size
    color = "black",          # Text color (keeps it readable)
    max.overlaps = 12,        # Prevents clutter by hiding overlapping labels
    show.legend = FALSE       # Prevents 'a' from appearing in the legend
  ) +
  # Visual formatting
  scale_fill_manual(values = c("Intra" = "#5e3c99",
                                "Inter" = "#d01c8b")) +
  scale_size_continuous(
    range = c(0.5, 5),                # Set the minimum and maximum point sizes
  ) + 
  labs(
    x = "IFN FACS score", 
    y = "NFkB FACS score",
    color = "Residue type",           # Custom title for the color legend
    size = "Metastatic\nfitness"             # Custom title for the size legend
  )  +
  theme_cowplot() +
  theme(panel.border = element_blank())

svglite("FACS_mets_bonds_IFN_vs_NFkB.svg", width = 8, height = 6)
print(p)
dev.off()

## Plot LOF groups by monomer/dimer counts
lof_df <- read.csv('../cGAMP_sc/LOF_variants.csv')
lof_bonds <- data.frame(
  residue = 1:378,
  monomer  = colSums(bonds_AA[, 1:378, drop = FALSE], na.rm = TRUE),
  dimer    = colSums(bonds_AB[, 1:378, drop = FALSE], na.rm = TRUE),
  tetramer = colSums(bonds_ACD[, 1:378, drop = FALSE], na.rm = TRUE)) %>%
  
  # 3. Calculate fractions and categorize groups using mean_fitness
  mutate(
    total_bonds       = monomer + dimer + tetramer,
    monomer_fraction  = ifelse(total_bonds > 0, monomer / total_bonds, 0),
    dimer_fraction    = ifelse(total_bonds > 0, dimer / total_bonds, 0),
    tetramer_fraction = ifelse(total_bonds > 0, tetramer / total_bonds, 0),
    
    # Classify based on your fractional rules
    bond_group = case_when(
      monomer_fraction > dimer_fraction+tetramer_fraction ~ "intra",
      dimer_fraction+tetramer_fraction >= monomer_fraction ~ "inter"
    )
  ) %>%
  left_join(lof_df, by = "residue") %>%
  mutate(group = case_when(
    is.na(group) ~ "Non-LOF",
    group=="IFNa & NFkB & UPR" ~ "Total LOF",
    TRUE ~ "Selective LOF"
  ) )
