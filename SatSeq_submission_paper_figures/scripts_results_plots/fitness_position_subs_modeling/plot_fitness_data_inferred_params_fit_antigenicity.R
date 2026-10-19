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
library(magick)
library(grid)
library(svglite)



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


plot_fitness_data_fit <- function(readout, panel_labels) {

# Load weights and format labels
weights_raw <- read.csv(paste0('inferred_weights_',readout,'.csv'))
colnames(weights_raw) <- c('position', 'additive', 'interaction')

weights_long <- weights_raw %>%
  mutate(position = position + 1) %>%
  pivot_longer(cols = c(additive, interaction), names_to = "effect", values_to = "weight") %>%
  group_by(effect) %>%
  mutate(
    rank_desc = min_rank(dplyr::desc(weight)),
    rank_asc = min_rank(weight),
    label = ifelse(rank_desc <= 10 | rank_asc <= 10, 
                   paste0(seq_vector[position], position), "")
  ) %>%
  ungroup()

# 2. PROCESS HEATMAP DATA -------------------------------------------------
# 1. MAP READOUTS TO CORRESPONDING BULK FITNESS FILES ----------------------
bulk_file_mapping <- list(
  "Initial_vs_Primary_4T1"   = "Mammary_Fat_Pad_Initial_vs_Primary_boot_4T1.csv",
  "Initial_vs_Primary_EO771" = "Mammary_Fat_Pad_Initial_vs_Primary_boot_EO771.csv",
  "Primary_vs_Met_4T1"       = "Mammary_Fat_Pad_Primary_vs_Met_boot_4T1.csv",
  "Tail_Vein_vs_Met_4T1" = "Tail_Vein_vs_Met_boot_4T1.csv",
  "Initial_vs_Tail_Vein_EO771" = "Initial_vs_Tail_Vein_boot_EO771.csv"
)

bulk_file <- bulk_file_mapping[[readout]]

# Handle invalid or missing readout names gracefully
if (is.null(bulk_file)) {
  stop(paste("Invalid readout provided. Choose from:", 
             paste(names(bulk_file_mapping), collapse = ", ")))
}

# 2. LOAD RAW BULK FITNESS DATA ---------------------------------------------
file_path <- paste0("../../datasets/bulk_fitness/", bulk_file)
fitness_data <- read_csv(file_path)

# 3. CLEAN AND PARSE DATASET -----------------------------------------------
# Valid 20 amino acids list
aa_list <- c('R','H','K','D','E','S','T','N','Q','C', 
             'G','P','A','V','I','L','M','F','Y','W')

fitness_data_clean <- fitness_data %>% 
  # Exclude control sequences
  filter(is_control == FALSE) %>% 
  # Parse variant syntax (e.g., "A123V")
  mutate(
    wt_residue = str_extract(variant, "^[A-Za-z]+"),
    position = as.numeric(str_extract(variant, "[0-9.]+")),
    mutant_residue = str_extract(variant, "[A-Za-z]+$")
  ) %>% 
  # Filter out failed regex parsing
  filter(!is.na(wt_residue), !is.na(position), !is.na(mutant_residue)) %>% 
  # Enforce standard amino acid constraints
  filter(wt_residue %in% aa_list, mutant_residue %in% aa_list) %>% 
  # Drop unmutated wild-type references
  filter(wt_residue != mutant_residue)

# 3. CONVERT TO HEATMAP MATRIX FORMAT ------------------------------------
# Average over replicate experiments
fitness <- fitness_data_clean  %>%
  group_by(position, mutant_residue) %>%
  summarize(
    mean_fit = mean(observed_fitness, na.rm = TRUE), 
    .groups = 'drop'
  ) 


# 3. GENERATE PLOTS -------------------------------------------------------

# --- TOP PLOT: Additive Weights ---
rd_bu_colors <- rev(brewer.pal(7, "RdBu"))
# Calculate symmetric limits based on 95th percentile to ensure "brightness"
heat_lim <- quantile(abs(fitness$mean_fit),
                     0.95, na.rm = TRUE)

p_add <- ggplot(filter(weights_long, effect == "additive"), 
                aes(x = position, y = weight, fill = weight)) +
  geom_col(show.legend = TRUE) +
  geom_text_repel(aes(label = label), size = 4, max.overlaps = Inf) +
  scale_fill_gradientn(colors = rd_bu_colors,
                       limits = c(-heat_lim, heat_lim),
                       oob = scales::squish,
                       name = NULL,
                       labels = function(x) sprintf("%.1f", x)) +
  labs(y = "Additive", x = NULL) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# --- MIDDLE PLOT: Interaction Weights ---
# Calculate symmetric limits based on 95th percentile to ensure "brightness"
w_int_lim <- quantile(abs(weights_long$weight[weights_long$effect == "interaction"]),
                      0.95, na.rm = TRUE)

p_int <- ggplot(filter(weights_long, effect == "interaction"), 
                aes(x = position, y = weight, fill = weight)) +
  geom_col(show.legend = TRUE) +
  geom_text_repel(aes(label = label), size = 4, max.overlaps = Inf) +
  scale_fill_gradientn(colors = rd_bu_colors,
                       limits = c(-heat_lim, heat_lim),
                       oob = scales::squish,
                       labels = function(x) sprintf("%.1f", x),
                       name = NULL) +
  labs(y = "Interaction", x = NULL) +
  theme_cowplot() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# --- BOTTOM PLOT: Symmetric Heatmap (Blue-White-Red) ---
# 1. Explicitly apply factor levels AND expand the grid to a full rectangle
wt_ref <- tibble(
  position = 1:length(seq_vector),
  mutant_residue = factor(seq_vector, levels = rev(aa_list))
)

# 2. Complete the grid for ALL explicit positions (1:378) and overwrite WT to 0
fitness_ordered <- fitness %>%
  mutate(mutant_residue = factor(mutant_residue, levels = rev(aa_list))) %>%
  # Explicitly state 1:378 to generate rows for entirely missing position columns
  complete(position = 1:378, mutant_residue, fill = list(mean_fit = NA)) %>%
  # Add a temporary flag indicating if the row matches the WT residue at that position
  left_join(wt_ref %>% mutate(is_wt = TRUE), by = c("position", "mutant_residue")) %>%
  # Overwrite mean_fit to 0 if it is a WT residue row, otherwise keep the existing mean_fit
  mutate(mean_fit = if_else(!is.na(is_wt) & is_wt, 0, mean_fit)) %>%
  select(-is_wt) # Clean up the temporary flag column



p_heat <- ggplot(fitness_ordered,
                 aes(x = position, y = mutant_residue, fill = mean_fit)) +
  geom_tile() +
  scale_fill_gradientn(colors = rd_bu_colors, 
                       limits = c(-heat_lim, heat_lim), # Force symmetry
                       oob = scales::squish,
                       labels = function(x) sprintf("%.1f", x),
                       name = ifelse(readout == "Initial_vs_Primary_4T1", "4T1 primary tumor\nfitness",
                                            ifelse(readout == "Initial_vs_Primary_EO771", "EO771 primary tumor\nfitness",
                                                   ifelse(readout == "Primary_vs_Met_4T1", "4T1 primary tumor\nmetastatic fitness",
                                                          ifelse(readout == "Tail_Vein_vs_Met_4T1", "4T1 Tail vein\nmetastatic fitness",
                                                                 ifelse(readout == "Initial_vs_Tail_Vein_EO771", "EO771 Tail vein", 
                                                                        readout))))),
                       # Color for missing data points
                       na.value = "gray80") +
  scale_x_continuous(expand = c(0,0)) +
  theme_cowplot() +
  labs(x = "Position", y = "Amino acid substitution") +
  # Customize the colorbar guide properties
  guides(fill = guide_colorbar(
    title.position = "top",                    # Position title to the side of the bar [1]
    title.theme = element_text(angle = 90,        # Rotate text 90 degrees [1]
                               hjust = 0.5,       # Center the text horizontally [1]
                               vjust = 0.5)       # Center the text vertically
  ))

# 4. ASSEMBLE -------------------------------------------------------------
p_positional <- (p_heat / p_add / p_int) + 
  plot_layout(heights = c(2, 1, 1)) & 
  scale_x_continuous(limits = c(0, 379),
                     breaks = seq(5, 375, 10),
                     expand = c(0.01, 0.01))

# Import 3d plots of weights
image_add <- image_read(paste0("inferred_weights_3dplot_", readout, "_additive.png")) |> 
  image_crop(geometry_area(width = 930, height = 1000))

image_add <- wrap_elements(rasterGrob(image_add)) + 
  theme(plot.margin = margin(0, 0, 0, 0, "pt"))

image_int <- image_read(paste0("inferred_weights_3dplot_", readout, "_interaction.png")) |> 
  image_crop(geometry_area(width = 930, height = 1000))

image_int <- wrap_elements(rasterGrob(image_int)) + 
  theme(plot.margin = margin(0, 0, 0, 0, "pt"))

# Plot aa matrices
# 1. SETUP COLORS AND AMINO ACID LIST -------------------------------------
rd_bu_colors <- rev(brewer.pal(7, "RdBu"))

# 2. FUNCTION TO PLOT MATRIX ----------------------------------------------
plot_aa_matrix <- function(file_path,heat_lim) {
  
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
  heat_lim <- quantile(abs(matrix_long$value),
                       0.95, na.rm = TRUE)
  # Plot
  ggplot(matrix_long, aes(x = factor(mut, levels = aa_list), 
                          y = factor(wt, levels = rev(aa_list)), 
                          fill = value)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_gradientn(
      colors = rd_bu_colors,
      limits = c(-heat_lim, heat_lim),
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

# 3. GENERATE THE TWO AA PLOTS -----------------------------------------------
# Plot Additive Matrix
p_mat_add <- plot_aa_matrix(
  paste0("inferred_aamat_additive_",readout,".csv"), heat_lim)

# Plot Interaction Matrix
p_mat_int <- plot_aa_matrix(
  paste0("inferred_aamat_interaction_",readout,".csv"), heat_lim)

# plot additive interaction histogram and fit
# Align weights column names to match your mathematical formula
additive_sub_matrix <- read.csv(
  paste0("inferred_aamat_additive_",readout,".csv"),
  row.names = 1)

multiplicative_sub_matrix <- read.csv(
  paste0("inferred_aamat_interaction_",readout,".csv"),
  row.names = 1)

# Get experiment intercept
experiment_intercept <- read.csv(
  paste0("inferred_intercept_experiment_",readout,".csv"),
  row.names = 1)


weights_clean <- weights_raw %>%
  mutate(position = position + 1) %>% # Ensure 1-based indexing
  dplyr::rename(multiplicative = interaction)

# Add wild-type residue mapping and compute fits/residuals row-by-row
fit <- fitness %>%
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
    additive_total = add_subs + add_pos,
    fit_add   = mean(experiment_intercept[,1]) + add_subs + add_pos,
    resid_add = mean_fit - fit_add,
    
    # Model 2: Additive + Interaction (Multiplicative)
    interaction_total = mult_subs * mult_pos,
    fit_full   = mean(experiment_intercept[,1]) + add_subs + add_pos + (mult_subs * mult_pos),
    resid_full = mean_fit - fit_full
  ) %>%
  ungroup() %>%
  mutate(mutant_residue = factor(mutant_residue, levels = c(
    'R','H','K','D','E','S','T','N','Q','C',
    'G','P','A','V','I','L','M','F','Y','W')))

# Plot histograms of additive and multiplicative effects
# Pivot to long format
fitness_long <- fit %>%
  select(additive_total, interaction_total) %>%
  dplyr::rename(additive = additive_total,
         interaction = interaction_total) %>%
  pivot_longer(
    cols = c(additive, interaction),
    names_to = "term",
    values_to = "value"
  )

# Plot histogram
p_hist <- ggplot(fitness_long, aes(x = value, fill = term)) +
  geom_histogram(aes(y = after_stat(density)),
                 alpha = 0.5,
                 position = "identity",
                 bins = 100
  ) +
  theme_cowplot() +
  theme(
    strip.text = element_text(margin = margin(t = 1, b = 1)),
    # Move legend inside the panel area
    legend.position = "inside",
    # Coordinates from 0 to 1: c(x, y). This puts it in the top right.
    legend.position.inside = c(0.95, 0.95), 
    # Anchor point of the legend box (1, 1 means top-right corner of the legend box)
    legend.justification = c(1, 1),
    # Optional: Removes the gray background box from the legend if it blocks data
    legend.background = element_blank() 
  ) +
  scale_fill_manual(values = c("additive" = "#880808",
                               "interaction" = "#0BDA51")) +
  labs(
    x = "",
    y = "Density",
    fill = ""
  )


# Plot fit
p_fit <- ggplot(fit, aes(x = mean_fit, y = fit_full)) +
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
  labs(
    x = "True fitness",
    y = "Fit")


## Assemble full plot per readout
# Step 1: Stack your square elements tightly
# Using heights = c(1.2, 1, 1) balances the histogram/fit plots with the squares
p_add_int <- (p_hist / image_add / image_int) + 
  plot_layout(heights = c(2, 1, 1))

p_aa_fit <- (p_fit / p_mat_add / p_mat_int) + 
  plot_layout(heights = c(2, 1, 1))

# Step 2: Combine the two square columns into a separate sub-grid
p_squares <- (p_add_int | p_aa_fit) + 
  plot_layout(widths = c(1, 1))

# Step 3: Bind the long plot to the square grid
# Adjust the widths ratio (e.g., 2:1 or 3:1) so the squares have breathing room
(p_positional | p_squares) + 
  plot_layout(widths = c(2.5, 1)) +
  # Adds sequential lowercase letters (a, b, c, d) to all plots
  plot_annotation(tag_levels = list(panel_labels)) &
  theme(
    plot.margin = margin(2, 2, 2, 2, "pt"),
    # Forces all subplot tags to be bold and slightly larger
    plot.tag = element_text(face = "bold", size = 14) 
  )
}

# Make figs per readout
p_primary_4T1 <- plot_fitness_data_fit("Initial_vs_Primary_4T1", c("a", "b", "", "c", "", "", "d", "e"))
p_primary_EO771 <- plot_fitness_data_fit("Initial_vs_Primary_EO771", c("f", "g", "", "h", "", "", "i", "j"))
p_mets_4T1 <- plot_fitness_data_fit("Primary_vs_Met_4T1", c("k", "l", "", "m", "", "", "n", "p"))
p_tailvein_4T1 <- plot_fitness_data_fit("Tail_Vein_vs_Met_4T1", c("q", "r", "", "s", "", "", "t", "u"))
p_tailvein_EO771 <- plot_fitness_data_fit("Initial_vs_Tail_Vein_EO771", c("v", "w", "", "x", "", "", "y", "z"))

## Make track of LOF marks
# 1. Load your LOF data
lof_data <- read.csv("../cGAMP_sc/LOF_variants.csv")

# Create a master dataframe for all 378 residues
# 1. Prepare data (ensure all 378 residues are present)
df_tracks <- data.frame(residue = 1:378) %>%
  left_join(lof_data, by = c("residue" = "residue")) %>%
  mutate(
    `IFN LOF`  = ifelse(grepl("IFNa", group), "TRUE", "FALSE"),
    `NFkB LOF` = ifelse(grepl("NFkB", group), "TRUE", "FALSE"),
    `UPR LOF`  = ifelse(grepl("UPR", group), "TRUE", "FALSE")
  ) %>%
  mutate(across(c(`IFN LOF`, `NFkB LOF`, `UPR LOF`), ~replace_na(., "FALSE"))) %>%
  pivot_longer(cols = c(`IFN LOF`, `NFkB LOF`, `UPR LOF`), names_to = "Pathway", values_to = "Status")

# Create the top panel containing just the 3 pixel rows
p_tracks <- ggplot(df_tracks, aes(x = residue, y = Pathway,
                                  fill = interaction(Pathway, Status))) +
  geom_tile(color = "white", linewidth = 0) + # Light grid border for crisp pixels
  scale_fill_manual(
    values = c(
      "UPR LOF.TRUE"  = "#35b6e3", "UPR LOF.FALSE"  = "white",
      "IFN LOF.TRUE"  = "#e3312a", "IFN LOF.FALSE"  = "white",
      "NFkB LOF.TRUE" = "#7ac943", "NFkB LOF.FALSE" = "white"
    ),
    guide = "none"
  ) +
  scale_x_continuous(limits = c(0, 379),
                     breaks = seq(5, 375, 10),
                     expand = c(0.01, 0.01)) +
  theme_cowplot() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(), # Hide X text since the bottom plot covers it
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )


p_tracks <- (p_tracks | plot_spacer()) + plot_layout(widths = c(1.95, 1))

## Load antigenicty data
p_fitness <- readRDS('../mutant_antigenicity/fitness_antigenicity_EO771_4T1.rds')

p_full <- (p_tracks / p_primary_4T1 / p_primary_EO771 / p_mets_4T1 / p_tailvein_4T1 / p_tailvein_EO771 / p_fitness) +
  plot_layout(heights = c(0.1, 1, 1, 1, 1, 1, 0.5)) +
  
  plot_annotation(tag_levels = list(c("a","", "b", "", "c", "", "", "d", "e","",
                                      "f", "g", "", "h", "", "", "i", "j", "",
                                      "k", "l", "", "m", "", "", "n", "p", "",
                                      "q", "r", "", "s", "", "", "t", "u", "",
                                      "u", "v", "", "w", "", "", "x", "y", "",
                                      "z", "aa", "bb"))) &
  theme(
    plot.margin = margin(2, 2, 2, 2, "pt"),
    # Forces all subplot tags to be bold and slightly larger
    plot.tag = element_text(face = "bold", size = 20) 
  )


# Save plot
ggsave(plot=p_full, 'Fitness_heatmap_weights_aamat_fit_antigenicity.pdf',
       width=20, height=36)


svglite("Fitness_heatmap_weights_aamat_fit_antigenicity.svg", width=20, height=36)
print(p_full)
dev.off()



