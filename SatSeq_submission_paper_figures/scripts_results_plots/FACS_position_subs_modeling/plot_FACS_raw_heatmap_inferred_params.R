rm(list=ls())#Clear environment

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(stringr)
library(patchwork)
library(ggrepel)
library(cowplot)
library(ggnewscale)
library(RColorBrewer)
library(scales)
library(svglite)
library(ggpmisc)
library(tidyverse)


setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

aa_list <- c('R','H','K','D','E','S','T','N','Q','C',
             'G','P','A','V','I','L','M','F','Y','W')

# 1. SETUP DATA -----------------------------------------------------------
wt_seq <- paste0("MPYSNLHPAIPRPRGHRSKYVALIFLVASLMILWVAKDPPNHTLKYLALHLASHELGLLL",
                 "KNLCCLAEELCHVQSRYQGSYWKAVRACLGCPIHCMAMILLSSYFYFLQNTADIYLSWMF",
                 "GLLVLYKSLSMLLGLQSLTPAEVSAVCEEKKLNVAHGLAWSYYIGYLRLILPGLQARIRM",
                 "FNQLHNNMLSGAGSRRLYILFPLDCGVPDNLSVVDPNIRFRDMLPQQNIDRAGIKNRVYS",
                 "NSVYEILENGQPAGVCILEYATPLQTLFAMSQDAKAGFSREDRLEQAKLFCRTLEEILED",
                 "VPESRNNCRLIVYQEPTDGNSFSLSQEVLRHIRQEEKEEVTMNAPMTSVAPPPSVLSQEP",
                 "RLLISGMDQPLPLRTDLI")
seq_vector <- strsplit(wt_seq, "")[[1]]


for (readout in c('STING','ATF3','Dual_IFN','NFkB')){

# Load weights and format labels
weights_raw <- read.csv(paste0('inferred_weights_FACS_4T1_',readout,'.csv'))
colnames(weights_raw) <- c('position', 'additive', 'interaction')

weights_long <- weights_raw %>%
  mutate(position = position + 1) %>%
  pivot_longer(cols = c(additive, interaction), names_to = "effect", values_to = "weight") %>%
  group_by(effect) %>%
  mutate(
    rank_desc = min_rank(desc(weight)),
    rank_asc = min_rank(weight),
    label = ifelse(rank_desc <= 10 | rank_asc <= 10, 
                   paste0(seq_vector[position], position), "")
  ) %>%
  ungroup()

# 2. PROCESS HEATMAP DATA -------------------------------------------------
# (Assumes df_up and df_down are loaded from previous steps)
# 1. LOAD RAW DATA --------------------------------------------------------
sheet_up   <- paste0("4T1_", readout, "_FACS_All_vs_Up")
sheet_down <- paste0("4T1_", readout, "_FACS_All_vs_Down")

df_up   <- read_excel("../../datasets/FACS/FACS_ATF3_IFN_NFkB_STING_Proliferation_4T1_EO771_all.xlsx", 
                      sheet = sheet_up)
df_down <- read_excel("../../datasets/FACS/FACS_ATF3_IFN_NFkB_STING_Proliferation_4T1_EO771_all.xlsx", 
                      sheet = sheet_down)

# 2. Assign regex patterns dynamically
patterns <- list(
  "ATF3"     = c("ATF3_R(\\d+)_Up", "ATF3_R(\\d+)_Dn"),
  "Dual_IFN" = c("R(\\d+)_Up",      "R(\\d+)_Down"),
  "STING"    = c("Up_(\\d+)",       "Down_(\\d+)"),
  "NFkB"     = c("NFkB_UP(\\d+)",   "NFkB_DN(\\d+)")
)

pattern_up   <- patterns[[readout]][1]
pattern_down <- patterns[[readout]][2]

# 2. MATCH PYTHON EXTRACTION, FILLNA, AND REPLICA MERGING -----------------
merged_fitness <- df_up %>%
  mutate(
    # Extract replicate dynamically using the loop variables
    rep = as.integer(str_extract(experiment, pattern_up, group = 1)),
    rep = coalesce(rep, 1L) # Use coalesce to match Python's fillna safely
  ) %>%
  inner_join(
    df_down %>% 
      mutate(
        rep = as.integer(str_extract(experiment, pattern_down, group = 1)),
        rep = coalesce(rep, 1L)
      ),
    by = c("variant", "rep", "is_control"), 
    suffix = c("_Up", "_Down")
  ) %>%
  # Filter out controls early to reduce downstream parsing workload
  filter(is_control == FALSE) %>%
  mutate(
    fitness_Up       = observed_fitness_Up,
    fitness_Down     = observed_fitness_Down,
    experiment       = paste0("Up_", rep, "_Down_", rep),
    observed_fitness = fitness_Up - fitness_Down,
    
    # Regex extraction for variant parsing
    variant_match  = str_match(variant, "^([A-Za-z])([0-9.]+)([A-Za-z])$")
  ) %>%
  # Critical Fix: Filter out failed regex matches before indexing matrix rows
  filter(!is.na(variant_match[, 1])) %>%
  mutate(
    wt             = variant_match[, 2],
    position       = as.numeric(variant_match[, 3]),
    mutant_residue = variant_match[, 4]
  ) %>%
  # Filter out invalid biological entries
  filter(
    wt != "X", 
    mutant_residue != "X",
    wt != mutant_residue,
    position >= 1,
    position <= 378
  ) %>%
  select(-variant_match) # Drop the temporary matrix helper column

# 3. CONVERT TO HEATMAP MATRIX FORMAT ------------------------------------
# Average over replicate experiments
sting_fitness <- merged_fitness  %>%
  group_by(position, mutant_residue) %>%
  summarize(
    mean_fit = mean(observed_fitness, na.rm = TRUE), 
    .groups = 'drop'
  ) 


# 3. GENERATE PLOTS -------------------------------------------------------

# --- TOP PLOT: Additive Weights (Bright Orange/Purple) ---
rd_bu_colors <- rev(brewer.pal(7, "RdBu"))
p_add <- ggplot(filter(weights_long, effect == "additive"), 
                aes(x = position, y = weight, fill = weight)) +
  geom_col(show.legend = TRUE) +
  geom_text_repel(aes(label = label), size = 3.5, max.overlaps = 15) +
  scale_fill_gradientn(colors = rd_bu_colors,
                       limits = c(-1, 1),
                         oob = scales::squish,
                       name = NULL) +
  labs(y = "Additive", x = NULL) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# --- MIDDLE PLOT: Interaction Weights (Bright Pink/Green) ---
p_int <- ggplot(filter(weights_long, effect == "interaction"), 
                aes(x = position, y = weight, fill = weight)) +
  geom_col(show.legend = TRUE) +
  geom_text_repel(aes(label = label), size = 3.5, max.overlaps = 15) +
  scale_fill_gradientn(colors = rd_bu_colors,
                       limits = c(-2.5, 2.5),
                       oob = scales::squish,
                       name = NULL) +
  labs(y = "Interaction", x = NULL) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# --- BOTTOM PLOT: Symmetric Heatmap (Blue-White-Red) ---
# Explicitly apply factor levels to mutant_residue
sting_fitness_ordered <- sting_fitness %>%
  mutate(mutant_residue = factor(mutant_residue, levels = aa_list))

p_heat <- ggplot(sting_fitness_ordered,
                 aes(x = position, y = mutant_residue, fill = mean_fit)) +
  geom_tile() +
  scale_fill_gradientn(colors = rd_bu_colors, 
                       limits = c(-5, 5), # Force symmetry
                       oob = scales::squish,
                       name = readout) +
  scale_x_continuous(expand = c(0,0)) +
  theme_cowplot() +
  labs(x = "Position", y = "Mutation")

# 4. ASSEMBLE -------------------------------------------------------------
final_layout <- (p_heat / p_add / p_int) + 
  plot_layout(heights = c(2.5, 1, 1)) & 
  scale_x_continuous(limits = c(0, length(seq_vector)), expand = c(0.01, 0.01))

# Save plot
ggsave(plot=final_layout, paste0('FACS_weights_heatmap_',readout,'.pdf'),
       width=12, height=7)

if (readout=='STING'){
svglite("FACS_weights_heatmap_STING.svg", width = 12, height = 7)
print(final_layout)
dev.off()
}


# Plot aa matrices
# 1. SETUP COLORS AND AMINO ACID LIST -------------------------------------
rd_bu_colors <- rev(brewer.pal(7, "RdBu"))

# 2. FUNCTION TO PLOT MATRIX ----------------------------------------------
plot_aa_matrix <- function(file_path) {
  
  # Read the CSV (assuming first column is row names)
  df <- read.csv(file_path, row.names = 1)
  
  # Ensure it is a full 20x20 by creating a long-form template
  matrix_long <- expand.grid(wt = aa_list, mut = aa_list) %>%
    left_join(
      df %>% 
        as.data.frame() %>% 
        mutate(wt = rownames(.)) %>% 
        pivot_longer(-wt, names_to = "mut", values_to = "value"),
      by = c("wt", "mut")
    ) %>%
    mutate(value = replace_na(value, 0)) # Set diagonal/missing to 0
  
  # Calculate symmetric limits based on 95th percentile to ensure "brightness"
  limit_val <- quantile(abs(matrix_long$value), 0.95, na.rm = TRUE)
  if(limit_val == 0) limit_val <- max(abs(matrix_long$value)) # Fallback
  
  # Plot
  ggplot(matrix_long, aes(x = factor(mut, levels = aa_list), 
                          y = factor(wt, levels = rev(aa_list)), 
                          fill = value)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_gradientn(
      colors = rd_bu_colors,
      limits = c(-limit_val, limit_val),
      oob = scales::squish,
      name = NULL
    ) +
    coord_fixed() + # Ensure cells are perfect squares
    labs(title = NULL, x = "Mutant Residue", y = "WT Residue") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 8),
      panel.grid = element_blank(),
      legend.title = element_text(angle = -90, hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_colorbar(title.position = "right",
                                 title.hjust = 0.5,
                                 barheight = unit(3, "cm"),
                                 barwidth = unit(0.3, "cm")))
}

# 3. GENERATE THE TWO PLOTS -----------------------------------------------

# Plot Additive Matrix
p_mat_add <- plot_aa_matrix(
  paste0("inferred_aamat_additive_FACS_4T1_",readout,".csv"))

# Plot Interaction Matrix
p_mat_int <- plot_aa_matrix(
  paste0("inferred_aamat_interaction_FACS_4T1_",readout,".csv"))

# 4. DISPLAY --------------------------------------------------------------
# Side-by-side display
p <- (p_mat_add / p_mat_int)

# Save plot
ggsave(plot=p, paste0('aa_mats_',readout,'.pdf'),
       width=3, height=6)

if (readout=='STING'){
svglite("aa_mats_STING.svg", width = 3, height = 6)
print(p)
dev.off()
}


# Align weights column names to match your mathematical formula
additive_sub_matrix <- read.csv(
  paste0("inferred_aamat_additive_FACS_4T1_",readout,".csv"),
                                row.names = 1)

multiplicative_sub_matrix <- read.csv(
  paste0("inferred_aamat_interaction_FACS_4T1_",readout,".csv"),
                                      row.names = 1)

# Get experiment intercept
experiment_intercept <- read.csv(
  paste0("inferred_intercept_experiment_FACS_4T1_",readout,".csv"),
                                 row.names = 1)


weights_clean <- weights_raw %>%
  mutate(position = position + 1) %>% # Ensure 1-based indexing
  rename(multiplicative = interaction)

# Add wild-type residue mapping and compute fits/residuals row-by-row
sting_models <- sting_fitness %>%
  mutate(
    position = as.integer(position),
    wt_residue = seq_vector[position]
  ) %>%
  rowwise() %>%
  mutate(
    # Matrix lookups
    add_subs  = additive_sub_matrix[wt_residue, mutant_residue],
    mult_subs = multiplicative_sub_matrix[wt_residue, mutant_residue],
    
    # Position weight lookups
    add_pos   = weights_clean$additive[weights_clean$position == position][1],
    mult_pos  = weights_clean$multiplicative[weights_clean$position == position][1],
    
    # Model 1: Additive Only
    fit_add   = mean(experiment_intercept[,1]) + add_subs + add_pos,
    resid_add = mean_fit - fit_add,
    
    # Model 2: Additive + Interaction (Multiplicative)
    fit_full   = mean(experiment_intercept[,1]) + add_subs + add_pos + (mult_subs * mult_pos),
    resid_full = mean_fit - fit_full
  ) %>%
  ungroup() %>%
  mutate(mutant_residue = factor(mutant_residue, levels = c(
  'R','H','K','D','E','S','T','N','Q','C',
  'G','P','A','V','I','L','M','F','Y','W')))

# 2. RESHAPE FOR GGPLOT FACETING -------------------------------------------

heatmap_data <- sting_models %>% 
  select(position, mutant_residue, fit_add, resid_add, fit_full, resid_full, mean_fit) %>% 
  pivot_longer(
    cols = c(fit_add, resid_add, fit_full, resid_full, mean_fit),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  mutate(
    # Create clean grouping variables
    model = case_when(
      metric == "mean_fit" ~ "True data",
      str_detect(metric, "_add") ~ "Additive",
      TRUE ~ "Additive + Interaction"
    ),
    type = case_when(
      metric == "mean_fit" ~ "True data",
      str_detect(metric, "fit_") ~ "Fitted Fitness",
      TRUE ~ "Residual"
    )
  ) %>%
  # FIX: Explicitly factorize to control the exact layout order in facet_wrap
  mutate(
    model = factor(model, levels = c("Additive", "Additive + Interaction", "True data")),
    type = factor(type, levels = c("True data", "Fitted Fitness", "Residual"))
  )

# 3. GENERATE THE 2x2 GRID OF HEATMAPS -------------------------------------

rd_bu_colors <- rev(brewer.pal(7, "RdBu"))

# Helper function to generate standardized heatmaps
plot_grid_heatmap <- function(data, type_filter, max_scale) {
  ggplot(filter(data, type == type_filter),
         aes(x = position, y = mutant_residue, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = rd_bu_colors, 
      limits = c(-5, 5),
      oob = scales::squish, 
      name = type_filter
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    facet_wrap(~ model, ncol = 1) +
    theme_cowplot() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      panel.spacing = unit(1, "lines")
    ) +
    labs(x = "Position", y = "Mutation")
}

# Generate both rows
p_fitted    <- plot_grid_heatmap(heatmap_data, "Fitted Fitness", max_fit)
p_residual  <- plot_grid_heatmap(heatmap_data, "Residual", max_res)

p_true <- ggplot(filter(heatmap_data, type == "True data"),
                 aes(x = position, y = mutant_residue, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rd_bu_colors, 
    limits = c(-5, 5),
    oob = scales::squish, 
    name = "Full data"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1, "lines")
  ) +
  labs(x = "Position", y = "Mutation")

# Stack rows together vertically using patchwork
final_grid <- (p_true / p_fitted + plot_layout(heights = c(1, 2.2))) | (p_true / p_residual + plot_layout(heights = c(1, 2.2)))

# Save plot
ggsave(plot=final_grid, paste0('FACS_fit_residual_heatmaps_',readout,'.pdf'),
       width=18, height=11)

if (readout=='STING'){
svglite("FACS_fit_residual_heatmaps_STING.svg", width=18, height=11)
print(final_grid)
dev.off()
}

# Plot dotplots of fit
p1 <- ggplot(sting_models,
             aes(x = mean_fit,
                 y = fit_add)) +
  geom_smooth(method = "lm", linewidth = 0.5, color = "#7393B3") +
  geom_point(alpha=0.5, size=0.8) +
  stat_fit_glance(method = "lm",
                  label.y = "bottom",
                  label.x = 'right',
                  mapping = aes(label = sprintf('R~"="~%.3f~","~italic(P)~"="~%.2g',
                                                after_stat(r.squared^0.5), after_stat(p.value)),
                                size=6),
                  parse = TRUE) +
  theme_cowplot() +
  labs(
    x = "True",
    y = "additive fit")

p2 <- ggplot(sting_models,
             aes(x = mean_fit,
                 y = fit_full)) +
  geom_smooth(method = "lm", linewidth = 0.5, color = "#7393B3") +
  geom_point(alpha=0.5, size=0.8) +
  stat_fit_glance(method = "lm",
                  label.y = "bottom",
                  label.x = 'right',
                  mapping = aes(label = sprintf('R~"="~%.3f~","~italic(P)~"="~%.2g',
                                                after_stat(r.squared^0.5), after_stat(p.value)),
                                size=6),
                  parse = TRUE) +
  theme_cowplot() +
  labs(
    x = "True",
    y = "full fit")

p <- p1|p2

print(c(var(sting_models$fit_add)/var(sting_models$mean_fit),
      var(sting_models$fit_full)/var(sting_models$mean_fit)))

print(cor.test(sting_models$fit_add,sting_models$mean_fit))
print(cor.test(sting_models$fit_full,sting_models$mean_fit))


}
