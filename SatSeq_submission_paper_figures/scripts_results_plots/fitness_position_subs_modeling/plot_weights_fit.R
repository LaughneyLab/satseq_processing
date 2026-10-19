rm(list=ls())#Clear environment

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
library(ggrastr)
library(patchwork)
library(ggrepel)
library(svglite)


# Process and plot Primary vs Initial weight data -----------------------------
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

## 4T1 weights
primary_weights_4T1 <- read.csv('inferred_weights_Initial_vs_Primary_4T1.csv')

colnames(primary_weights_4T1) <- c('position', 'additive','multiplicative')

primary_weights_4T1 <- primary_weights_4T1 %>%
  pivot_longer(!position, names_to = 'effect', values_to = 'weight')

primary_weights_4T1$position <- primary_weights_4T1$position +1

primary_weights_4T1$cell_line <- '4T1'

## EO771 weights
primary_weights_EO771 <- read.csv('inferred_weights_Initial_vs_Primary_EO771.csv')

colnames(primary_weights_EO771) <- c('position', 'additive','multiplicative')

primary_weights_EO771 <- primary_weights_EO771 %>%
  pivot_longer(!position, names_to = 'effect', values_to = 'weight')

primary_weights_EO771$position <- primary_weights_EO771$position +1

primary_weights_EO771$cell_line <- 'EO771'

# bind both cell types
primary_weights <- rbind(primary_weights_EO771,primary_weights_4T1)

# Discard interaction weights
primary_weights <- primary_weights[primary_weights$effect=='additive',]

### Plot weight barplots
p <- ggplot(primary_weights,
            aes(x = position,
                y = weight,
                fill = weight)) +
  geom_col() +
  scale_fill_distiller(
    name = NULL,
    type = "div",
    palette = "RdBu",
    direction = -1,
    aesthetics = "fill",
    limits = c(-0.1, 0.1),oob = scales::squish) +
  facet_grid(cell_line~effect,scales = 'free_y')+
  theme_cowplot() +
  theme(strip.text = element_text(margin = margin(t = 1, b = 1))) +
  theme(strip.text.x = element_blank()) +
  labs(x = "Position",
       y = "Weight") +
  theme(axis.text.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=15))

# Save plot
ggsave(plot=p, 'weight_barplot_primary_EO771_4T1.pdf',
        width=8, height=4)

# Plot correlation between cell lines for Primary --------------------------------- 
# Change to wide form
primary_weights <- primary_weights[c('position','weight','cell_line')] %>%
  pivot_wider(
    names_from = cell_line, 
    values_from = weight
  )

# take moving average of 5 residues up and downstream
primary_weights$moving_average_EO771 <- stats::filter(primary_weights$EO771,
                                                      filter = rep(1/11, 11),
                                                      sides = 2)

primary_weights$moving_average_4T1 <- stats::filter(primary_weights$`4T1`,
                                                      filter = rep(1/11, 11),
                                                      sides = 2)

# Plot correlation between cell lines
p <- ggplot(primary_weights,
            aes(x = moving_average_4T1,
                y = moving_average_EO771)) +
  geom_smooth(method = "lm",linewidth = 0.5,color='black') +
  geom_point(alpha=0.65,size=0.8) +
  stat_fit_glance(method = "lm",
                  label.y = "bottom",
                  label.x = 'right',
                  mapping = aes(label = sprintf('R~"="~%.3f~","~italic(P)~"="~%.2g',
                                                after_stat(r.squared^0.5), after_stat(p.value)),
                                size=20),
                  parse = TRUE) +
  theme_cowplot() +
  labs(x = "4T1 weight",
       y = "EO771 weight") +
  theme(axis.text.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=15))

# Save plot
ggsave(plot=p, 'weight_correlation_primary_EO771_4T1.pdf',
       width=3, height=4)



# Plot fit and additive/interaction effects in Primary -----------------------------

fitness_true_fit <- data.frame()

for (cell_line in c('4T1','EO771')){

# Load true fitness data
fitness_data <- read_csv(paste0("../../datasets/bulk_fitness/Mammary_Fat_Pad_Initial_vs_Primary_boot_",
                                cell_line,
                              ".csv"))

# Exclude controls
fitness_data <- fitness_data %>%
  filter(is_control == FALSE)

# Parse variant column (e.g., A123V)
fitness_data <- fitness_data %>%
  mutate(
    wt_residue = str_extract(variant, "^[A-Za-z]+"),
    position = as.numeric(str_extract(variant, "[0-9.]+")),
    mutant_residue = str_extract(variant, "[A-Za-z]+$")
  )

# Drop invalid parses (e.g., "_379X")
fitness_data <- fitness_data %>%
  filter(!is.na(wt_residue),
         !is.na(position),
         !is.na(mutant_residue))

# Amino acid list
aa_list <- c('R','H','K','D','E','S','T','N','Q','C',
             'G','P','A','V','I','L','M','F','Y','W')

# Keep only valid amino acids
fitness_data <- fitness_data %>%
  filter(wt_residue %in% aa_list,
         mutant_residue %in% aa_list)

# Drop unmutated residues
fitness_data <- fitness_data %>%
  filter(wt_residue != mutant_residue)

# Load inferred weight and AA subs files
weights <- read.csv(paste0('inferred_weights_Initial_vs_Primary_',cell_line,'.csv'))

colnames(weights) <- c('position', 'additive','multiplicative')

additive_sub_matrix <- read.csv(paste0("inferred_aamat_additive_Initial_vs_Primary_",cell_line,".csv"),
                                       row.names = 1)
multiplicative_sub_matrix <- read.csv(paste0("inferred_aamat_interaction_Initial_vs_Primary_",cell_line,".csv"),
                                      row.names = 1)
# Load experiment-specific intercept
experiment_intercept <- read.csv(paste0('inferred_intercept_experiment_Initial_vs_Primary_',cell_line,'.csv'),
                                 row.names = 1)
colnames(experiment_intercept) <- c('effect')
# Compute fitted fitness

fitness_data <- fitness_data %>%
  rowwise() %>%
  mutate(
    add_subs = additive_sub_matrix[wt_residue, mutant_residue],
    mult_subs = multiplicative_sub_matrix[wt_residue, mutant_residue],
    
    add_pos = weights$additive[position],
    mult_pos = weights$multiplicative[position],
    
    mult_total = mult_subs * mult_pos,
    
    experiment_effect = experiment_intercept[experiment,],
    
    fitted_fitness = experiment_effect + add_subs + add_pos + (mult_subs * mult_pos)
  ) %>%
  ungroup()

fitness_data$cell_line <- cell_line

fitness_true_fit <- rbind(fitness_true_fit,
                          fitness_data[c('observed_fitness',
                                         'fitted_fitness',
                                         'cell_line','add_pos',
                                         'mult_total')])

}

# Plot fits
### Plot true and fitted weights
p <- ggplot(fitness_true_fit,
            aes(x = observed_fitness,
                y = fitted_fitness)) +
  rasterise(geom_point(alpha=0.1,size=0.02), dpi = 300) +
  stat_fit_glance(method = "lm",
                  label.y = "bottom",
                  label.x = 'right',
                  mapping = aes(label = sprintf('R~"="~%.3f',
                                                after_stat(r.squared^0.5)),
                                size=20),
                  parse = TRUE) +
  facet_wrap(vars(cell_line),scales = 'free',
             ncol=1,strip.position = "right")+
  theme_cowplot() +
  theme(strip.text = element_text(margin = margin(t = 1, b = 1))) +
  labs(x = "Observed fitness",
       y = "Fitted fitness") +
  theme(axis.text.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=18))

# Save plot
ggsave(plot=p, 'fits_primary_EO771_4T1.pdf',
       width=3, height=3)

# Plot histograms of additive and multiplicative effects
# Pivot to long format
fitness_long <- fitness_true_fit %>%
  select(cell_line, add_pos, mult_total) %>%
  rename(additive = add_pos,
         interaction = mult_total) %>%
  pivot_longer(
    cols = c(additive, interaction),
    names_to = "term",
    values_to = "value"
  )

# Plot histogram
p <- ggplot(fitness_long, aes(x = value, fill = term)) +
  geom_histogram(aes(y = after_stat(density)),
    alpha = 0.5,
    position = "identity",
    bins = 100
  ) +
  facet_wrap(vars(cell_line),scales = 'free',
             ncol=1,strip.position = "right")+
  theme_cowplot() +
  theme(strip.text = element_text(margin = margin(t = 1, b = 1))) +
  scale_fill_manual(values = c("additive" = "#7b3294",
                               "interaction" = "#008837")) +
  labs(
    x = "",
    y = "Density",
    fill = "Term"
  ) +
  theme(axis.text.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=15))

# Save plot
ggsave(plot=p, 'term_histogram_primary_EO771_4T1.pdf',
       width=4, height=4)





# Plot weights and aa matrix for Primary vs Mets data for 4T1------------------------
primary_mets_weights <- read.csv('inferred_weights_Primary_vs_Met_4T1.csv')
colnames(primary_mets_weights) <- c('position', 'additive','interaction')

primary_mets_weights$position <- primary_mets_weights$position +1

# Add annotation for peak mets
primary_mets_weights$label <- NA
primary_mets_weights$label[c(13,337,146,78,209,181,107,198)] <- 
  c('P13','K337','V146','Q78','D209','F181','F107','Y198')

# Plot
p <- ggplot(primary_mets_weights,
             aes(x = additive,
                 y = abs(interaction))) +
  geom_smooth(method = "lm",linewidth = 0.5,color='black') +
  geom_point(alpha=0.65,size=0.8) +
  stat_fit_glance(method = "lm",
                  label.y = "top",
                  label.x = 'left',
                  mapping = aes(label = sprintf('R~"="~%.3f~","~italic(P)~"="~%.2g',
                                                after_stat(r.squared^0.5), after_stat(p.value)),
                                size=20),
                  parse = TRUE) +
  geom_text_repel(aes(label = label),
                  size = 5,max.overlaps = 20,alpha=0.8,fill = "transparent",
                  label.size = 0,
                  min.segment.length = unit(0.1, 'lines'),show.legend=FALSE) +
  theme_cowplot() +
  labs(x = 'Additive',
       y = "|Interaction|") +
  theme(axis.text.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=15))

# Save plot
ggsave(plot=p, 'weight_additive_vs_interaction_Primary_vs_Met_4T1.pdf',
       width=4, height=4)


# Pivot longer
primary_mets_weights <- read.csv('inferred_weights_Primary_vs_Met_4T1.csv')
colnames(primary_mets_weights) <- c('position', 'additive','interaction')

primary_mets_weights$position <- primary_mets_weights$position +1

primary_mets_weights <- primary_mets_weights %>%
  pivot_longer(!position, names_to = 'effect', values_to = 'weight')
primary_mets_weights$label <- NA
primary_mets_weights$label[2*c(13,337,146,78,209,181,107,198)-1] <- 
  c('P13','K337','V146','Q78','D209','F181','F107','Y198')

p1 <- ggplot(primary_mets_weights[primary_mets_weights$effect=='additive',],
            aes(x = position,
                y = weight,
                fill = weight)) +
  geom_col() +
  geom_text(aes(label = label), vjust = -0.04) +
  scale_fill_distiller(
    name = NULL,
    type = "div",
    palette = "RdBu",
    direction = -1,
    aesthetics = "fill",
    limits = c(-0.1, 0.1),oob = scales::squish
    ) +
  theme_cowplot() +
  labs(x = NULL,
       y = "Additive") +
  theme(axis.text.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=15))

p2 <- ggplot(primary_mets_weights[primary_mets_weights$effect=='interaction',],
             aes(x = position,
                 y = weight,
                 fill = weight)) +
  geom_col() +
  scale_fill_distiller(
    name = NULL,
    type = "div",
    palette = "RdBu",
    direction = -1,
    aesthetics = "fill",
    limits = c(-0.5, 0.5),oob = scales::squish
  ) +
  theme_cowplot() +
  labs(x = "Position",
       y = "Interaction") +
  theme(axis.text.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=15))

# Stack plots
p <- p1 / p2

# Save plot
ggsave(plot=p, 'weight_barplot_Primary_vs_Met_4T1.pdf',
       width=8, height=4)

## Plot AA matrix now
# Load AA matrix
aa_matrix <- read.csv('inferred_aamat_additive_Primary_vs_Met_4T1.csv',
                   row.names = 1)

heat_df <- as.data.frame(aa_matrix) %>%
  mutate(Row = rownames(.)) %>%
  pivot_longer(-Row, names_to = "Col", values_to = "Value")

# Keep row order fixed
heat_df$Row <- factor(heat_df$Row, levels = rev(aa_list))
heat_df$Col <- factor(heat_df$Col, levels = aa_list)

heatmap_plot <- ggplot(heat_df, aes(x = Col, y = Row, fill = Value)) +
geom_tile() +
  scale_fill_distiller(
    type = "div",
    palette = "RdBu",
    direction = -1,
    limits = c(-0.06, 0.06),oob = scales::squish,
    name = NULL
  ) +
  coord_equal() +  # makes it square
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
  theme(legend.position = "left") +
  labs(x = "Mutated AA",
       y = "WT AA") +
  theme(axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12))

# Add row means barplot
row_means <- rowMeans(aa_matrix, na.rm = TRUE)

bar_df <- data.frame(
  Row = factor(names(row_means),
               levels = levels(heat_df$Row)),
  Mean = row_means)

bar_plot <- ggplot(bar_df, aes(x = Mean, y = Row, fill=Mean)) +
  geom_col() +
  scale_fill_distiller(
    name = NULL,
    type = "div",
    palette = "RdBu",
    direction = -1,
    aesthetics = "fill",
    limits = c(-0.01, 0.01),oob = scales::squish,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(-0.01, 0.01),      # Keeps ticks exactly at these two points
    limits = c(-0.012, 0.012)     # Sets the outer borders of the x-axis
  ) +
  theme_cowplot() +
  theme(axis.text.x = element_text(size=12),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.y = element_blank()
  )

# Combine with patchwork
p <- (heatmap_plot | bar_plot) + plot_layout(widths = c(6, 1))

# Save plot
ggsave(plot=p, 'aa_matrix_Primary_vs_Met_4T1.pdf',
       width=6, height=4)

svglite("aa_matrix_Primary_vs_Met_4T1.svg", width=6, height=4)
print(p)
dev.off()



# Plot fit, and additive/interaction effects for 4T1 Primary vs Mets-------
# Load true fitness data
fitness_data <- read.csv("../../datasets/bulk_fitness/Mammary_Fat_Pad_Primary_vs_Met_boot_4T1.csv")

# Exclude controls
fitness_data <- fitness_data %>%
  filter(is_control == 'False')

# Parse variant column (e.g., A123V)
fitness_data <- fitness_data %>%
  mutate(
    wt_residue = str_extract(variant, "^[A-Za-z]+"),
    position = as.numeric(str_extract(variant, "[0-9.]+")),
    mutant_residue = str_extract(variant, "[A-Za-z]+$")
  )

# Drop invalid parses (e.g., "_379X")
fitness_data <- fitness_data %>%
  filter(!is.na(wt_residue),
         !is.na(position),
         !is.na(mutant_residue))

# Amino acid list
aa_list <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
             'S','T','V','W','Y')

# Keep only valid amino acids
fitness_data <- fitness_data %>%
  filter(wt_residue %in% aa_list,
         mutant_residue %in% aa_list)

# Drop unmutated residues
fitness_data <- fitness_data %>%
  filter(wt_residue != mutant_residue)

# Load inferred weight and AA subs files
weights <- read.csv('inferred_weights_Primary_vs_Met_4T1.csv')

colnames(weights) <- c('position', 'additive','multiplicative')

additive_sub_matrix <- read.csv("inferred_aamat_additive_Primary_vs_Met_4T1.csv",
                                row.names = 1)
multiplicative_sub_matrix <- read.csv("inferred_aamat_interaction_Primary_vs_Met_4T1.csv",
                                      row.names = 1)
# Load experiment-specific intercept
experiment_intercept <- read.csv('inferred_intercept_experiment_Primary_vs_Met_4T1.csv',
                                 row.names = 1)
colnames(experiment_intercept) <- c('effect')
# Compute fitted fitness

fitness_data <- fitness_data %>%
  rowwise() %>%
  mutate(
    add_subs = additive_sub_matrix[wt_residue, mutant_residue],
    mult_subs = multiplicative_sub_matrix[wt_residue, mutant_residue],
    
    add_pos = weights$additive[position],
    mult_pos = weights$multiplicative[position],
    
    mult_total = mult_subs * mult_pos,
    
    experiment_effect = experiment_intercept[experiment,],
    
    fitted_fitness_add = experiment_effect + add_subs + add_pos,
    
    fitted_fitness = experiment_effect + add_subs + add_pos + (mult_subs * mult_pos)
  ) %>%
  ungroup()


# Calculate variance explained by additive and full models
var_df <- fitness_data %>%
  group_by(position, residue) %>%
  summarise(
    across(c(observed_fitness, fitted_fitness_add, fitted_fitness),
           \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

print(c(var(var_df$fitted_fitness_add)/var(var_df$observed_fitness),
        var(var_df$fitted_fitness)/var(var_df$observed_fitness)))

print(cor.test(var_df$fitted_fitness_add,var_df$observed_fitness))
print(cor.test(var_df$fitted_fitness,var_df$observed_fitness))


# Plot fits
### Plot true and fitted weights
p <- ggplot(fitness_data,
            aes(x = observed_fitness,
                y = fitted_fitness)) +
  rasterise(geom_point(alpha=0.1,size=0.02), dpi = 600) +
  stat_fit_glance(method = "lm",
                  label.y = "bottom",
                  label.x = 'right',
                  mapping = aes(label = sprintf('R~"="~%.3f',
                                                after_stat(r.squared^0.5)),
                                size=20),
                  parse = TRUE) +
  theme_cowplot() +
  labs(x = "Observed fitness",
       y = "Fitted fitness") +
  theme(axis.text.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=15))

# Save plot
ggsave(plot=p, 'fits_Primary_vs_Met_4T1.pdf',
       width=4, height=2)

# Plot histograms of additive and multiplicative effects
# Pivot to long format
fitness_long <- fitness_data %>%
  select(add_pos, mult_total) %>%
  rename(additive = add_pos,
         interaction = mult_total) %>%
  pivot_longer(
    cols = c(additive, interaction),
    names_to = "term",
    values_to = "value"
  )

# Plot histogram
p <- ggplot(fitness_long, aes(x = value, fill = term)) +
  geom_histogram(aes(y = after_stat(density)),
                 alpha = 0.5,
                 position = "identity",
                 bins = 100
  ) +
  theme_cowplot() +
  scale_fill_manual(values = c("additive" = "#7b3294",
                               "interaction" = "#008837")) +
  labs(
    x = "",
    y = "Density",
    fill = "Term"
  ) +
  theme(axis.text.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=15))

# Save plot
ggsave(plot=p, 'term_histogram_Primary_vs_Met_4T1.pdf',
       width=4, height=2)







