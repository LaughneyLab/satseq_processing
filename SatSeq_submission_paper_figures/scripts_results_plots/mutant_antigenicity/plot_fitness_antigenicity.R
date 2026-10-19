rm(list=ls())#Clear environment

#################
## libraries ####
#################
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(readxl)
library(viridis)
library(stringr)
library(cowplot)
library(forcats)
library(stringdist)
library(stringr)
library(purrr)
library(ggrepel)
library(ggpmisc)
library(patchwork)
library(svglite)

### Gather and process data ####
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# Plot antigenicity and primary fitness data for 4T1 and EO771-----------------
antigenicity_df <- data.frame()

for (cell_line in c('EO771','4T1')){

# Mouse strain for correct MHC for each cell line
if (cell_line=='4T1'){
  df <- read.csv('mhcflurry_outputs_BALBc.csv',row.names = 1)
}
  
if (cell_line=='EO771'){
  df <- read.csv('mhcflurry_outputs_B6.csv',row.names = 1)
}

# Add a col for peptide length
df$length <- nchar(df$peptide)

# Add WT values to all rows
# Columns we want to pull from the WT row
cols_to_copy <- c(
  "sequence_name", "pos", "peptide", "n_flank", "c_flank",
  "affinity", "best_allele", "affinity_percentile",
  "processing_score", "presentation_score", "presentation_percentile","length"
)

# Make WT and MUT subsets
df_WT  <- df[df$sequence_name == "WT", ]
df_MUT <- df[df$sequence_name != "WT", ]

# Rename WT columns so they get WT_ prefix later
names(df_WT)[names(df_WT) %in% cols_to_copy] <- paste0("WT_", cols_to_copy)

# Perform merge on flanks + peptide
epitope_data <- merge(
  df_MUT,
  df_WT,
  by.x = c("pos","n_flank", "c_flank","length"),
  by.y = c("WT_pos","WT_n_flank", "WT_c_flank","WT_length"),
  all.x = TRUE
)

# Add a column on difference between presentation score
epitope_data$diff_presentation_score <- 
  -log(epitope_data$affinity/epitope_data$WT_affinity)

# Add Primary fitness for mutants
fitness_data <- read_csv(paste0("../../datasets/bulk_fitness/Mammary_Fat_Pad_Initial_vs_Primary_boot_",
                                cell_line,
                                ".csv"))
# Average over mice
fitness_data <- fitness_data %>%
  group_by(variant) %>%
  summarise(mean_fitness = mean(observed_fitness, na.rm = TRUE),
            position = unique(position, na.rm = TRUE)) 

rownames(fitness_data) <- fitness_data$variant


epitope_data$fitness <- unlist(array(fitness_data[epitope_data$sequence_name,'mean_fitness']))

# Take max over all epitopes from one mutation
epitope_data <- epitope_data %>%
  group_by(sequence_name) %>%
  summarize(max_diff_presentation = max(diff_presentation_score),
            fitness = unique(fitness))

# Add to df
epitope_data$cell_line <- cell_line
antigenicity_df <- rbind(antigenicity_df,epitope_data)
}


## Plotting
p1 <- ggplot(antigenicity_df,
            aes(x = fitness,
                y = max_diff_presentation)) +
  # 1. Filled 2D Density Contours with a Power Transformation (The "Blob")
  stat_density_2d(
    geom = "polygon",
    contour = TRUE,
    aes(alpha = after_stat((level)^2)), 
    bins = 100, 
    fill = "gray40",
    show.legend = FALSE
  ) +
  
  # 2. Extreme outer contour lines to highlight the overall data boundary
  geom_density_2d(
    bins = 10, 
    linewidth = 0.5, 
    color = "black", 
    alpha = 1, 
    show.legend = FALSE
  ) +
  # 5. Multi-line Correlation showing BOTH signed R and P value
  stat_correlation(
    method = "pearson",
    label.y = "top",
    label.x = "right",
    output.type = "text", # Disables complex plotmath to prevent syntax crashes
    mapping = aes(label = paste(after_stat(r.label), 
                                after_stat(p.value.label), 
                                sep = "\n")),
    size = 4
  ) +
  theme_cowplot() +
  guides(fill = "none") + # Turn off density fill legend
  theme(legend.title = element_blank()) +
  facet_wrap(vars(cell_line),scales = 'free',
             nrow=1,strip.position = "top")+
  theme(strip.text = element_text(margin = margin(t = 1, b = 1))) +
  labs(x = "Primary tumor fitness",
       y = "Predicted antigenicity") +
  theme(axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15))

# Plot correlation between 4T1 and EO771 antigenicity

# Change to wide form
antigenicity_df <- 
  antigenicity_df[c('max_diff_presentation','cell_line','sequence_name')] %>%
  pivot_wider(names_from = 'cell_line',
              values_from = 'max_diff_presentation')

# Plot antigenicity correlation
p2 <- ggplot(antigenicity_df,
            aes(x = `4T1`,
                y = EO771)) +
  geom_smooth(method = "lm",linewidth = 0.5,color='black') +
  # 1. Filled 2D Density Contours with a Power Transformation (The "Blob")
  stat_density_2d(
    geom = "polygon",
    contour = TRUE,
    aes(alpha = after_stat((level)^2)), 
    bins = 100, 
    fill = "gray40",
    show.legend = FALSE
  ) +
  
  # 2. Extreme outer contour lines to highlight the overall data boundary
  geom_density_2d(
    bins = 10, 
    linewidth = 0.5, 
    color = "black", 
    alpha = 1, 
    show.legend = FALSE
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
  theme(legend.title = element_blank()) +
  labs(x = "4T1 antigenicity",
       y = "EO771\nantigenicity") +
  theme(axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15))

# Plot antigenicity and mets fitness data for 4T1-----------------
# Add Mets fitness for mutants
fitness_data <- read_csv("../../datasets/bulk_fitness/Mammary_Fat_Pad_Primary_vs_Met_boot_4T1.csv")
# Average over mice
fitness_data <- fitness_data %>%
  group_by(variant) %>%
  summarise(mean_fitness = mean(observed_fitness, na.rm = TRUE),
            position = unique(position, na.rm = TRUE)) 

rownames(fitness_data) <- fitness_data$variant


epitope_data$fitness <- unlist(array(fitness_data[epitope_data$sequence_name,'mean_fitness']))

## Add a threshold for labeling in plots
epitope_data$is_labeled <- 10*epitope_data$fitness * epitope_data$max_diff_presentation
epitope_data$label <- NA
epitope_data$label <- ifelse(epitope_data$is_labeled > 7.5, 
                             epitope_data$sequence_name, 
                             NA)  

epitope_data$is_immunogenic <- "No"
epitope_data$is_immunogenic[epitope_data$is_labeled>15.5] <- "Yes"

# Add Ref and STOP codons
epitope_data <- bind_rows(epitope_data,
                          data.frame(sequence_name=c('Ref','Stop'),
                                     max_diff_presentation=c(0,0),
                                     fitness=c(0,fitness_data$mean_fitness[fitness_data$variant=='STOP']),
                                     label=c('Ref','Stop')))
# Plot
p3 <- ggplot(epitope_data,
            aes(x = fitness,
                y = max_diff_presentation)) +
  # 1. Filled 2D Density Contours with a Power Transformation (The "Blob")
  stat_density_2d(
    geom = "polygon",
    contour = TRUE,
    aes(alpha = after_stat((level)^2)), 
    bins = 100, 
    fill = "gray40",
    show.legend = FALSE
  ) +
  
  # 2. Extreme outer contour lines to highlight the overall data boundary
  geom_density_2d(
    bins = 10, 
    linewidth = 0.5, 
    color = "black", 
    alpha = 1, 
    show.legend = FALSE
  ) +
  # 5. Multi-line Correlation showing BOTH signed R and P value
  stat_correlation(
    method = "pearson",
    label.y = "top",
    label.x = 'right',
    mapping = aes(label = paste("atop(",
                                after_stat(r.label),
                                ",", after_stat(p.value.label), ")")),
    size = 4,
    parse = TRUE # Evaluates atop() dynamically as stacked text
  )  +
  theme_cowplot() +
  guides(fill = "none") + # Turn off density fill legend
  theme(legend.title = element_blank()) +
  labs(x = "Metastatic fitness",
       y = "Predicted antigenicity") +
  theme(axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15))

# Save all plots together
# Combine plots
p <- ((p1 | p2 | p3) + plot_layout(widths = c(2, 1, 1))) + plot_annotation(
  tag_levels = list(c('d', 'e', 'f')),
  theme = theme(plot_tag = element_text(size = 15, face = "bold")))

# Save plot
ggsave(plot=p, 'fitness_antigenicity_EO771_4T1.pdf',
       width=12, height=3)

saveRDS(p, file = "fitness_antigenicity_EO771_4T1.rds")


# Save density plot
epitope_data_filter <- epitope_data %>%
  filter(cell_line=='4T1') %>%
  mutate(position = as.numeric(str_extract(sequence_name, "[0-9.]+"))) %>%
  mutate(label = case_when(
   position %in% c(13,337)| sequence_name %in% c('Q78Y','E148L','V146Y','C88P',
                                                 'D204P','C205V','E149I','K151F',
                                                 'L356F','N236Y','T139F','S194F',
                                                 'P316M')~ sequence_name,
    TRUE            ~ NA_character_))

p4 <- ggplot(epitope_data_filter, aes(x = fitness, y = max_diff_presentation)) + 
  
  # 1. Filled 2D Density Contours with a Power Transformation (The "Blob")
  stat_density_2d(
    geom = "polygon",
    contour = TRUE,
    aes(alpha = after_stat((level)^2)), 
    bins = 100, 
    fill = "gray40",
    show.legend = FALSE
  ) +
  
  # 2. Extreme outer contour lines to highlight the overall data boundary
  geom_density_2d(
    bins = 10, 
    linewidth = 0.5, 
    color = "black", 
    alpha = 1, 
    show.legend = FALSE
  ) +
   
  # 3. Layer individual highlighted points explicitly over the density contours
  geom_point(
    data = epitope_data_filter[!is.na(epitope_data_filter$label), ], 
    aes(x = fitness, y = max_diff_presentation), 
    size = 1, 
    shape = 16 
  ) +
  
  # 4. Text labels
  geom_text_repel(
    aes(label = label), 
    size = 3, 
    max.overlaps = 30, 
    alpha = 0.8, 
    fill = "transparent", 
    label.size = 0, 
    min.segment.length = unit(0.1, 'lines'), 
    show.legend = FALSE
  ) +
  
  theme_cowplot() + 
  labs(x = "Metastatic fitness", y = "Predicted antigenicity", color = NULL) + 
  theme(
    axis.text.x = element_text(size = 15), 
    axis.title.y = element_text(size = 15), 
    axis.title.x = element_text(size = 15), 
    axis.text.y = element_text(size = 15), 
    legend.text = element_text(size = 15)
  ) +
  ylim(-0.2, 4.8) +
  xlim(-0.2, 0.65)

# Save plot
svglite("4T1_Fitness_vs_antigenicity.svg", width=4, height=4)
print(p4)
dev.off()

