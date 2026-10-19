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
library(ggseqlogo)
library(Biostrings)
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


plot_FACS_data_fit <- function(readout, panel_labels) {

# Load weights and format labels
weights_raw <- read.csv(paste0('inferred_weights_FACS_4T1_',readout,'.csv'))
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
facs_fitness <- merged_fitness  %>%
  group_by(position, mutant_residue) %>%
  summarize(
    mean_fit = mean(observed_fitness, na.rm = TRUE), 
    .groups = 'drop'
  ) 


# 3. GENERATE PLOTS -------------------------------------------------------

# --- TOP PLOT: Additive Weights ---
rd_bu_colors <- rev(brewer.pal(7, "RdBu"))
# Calculate symmetric limits based on 95th percentile to ensure "brightness"
heat_lim <- quantile(abs(facs_fitness$mean_fit),
                     0.95, na.rm = TRUE)

p_add <- ggplot(filter(weights_long, effect == "additive"), 
                aes(x = position, y = weight, fill = weight)) +
  geom_col(show.legend = TRUE) +
  geom_text_repel(aes(label = label), size = 4, max.overlaps = 15) +
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
  geom_text_repel(aes(label = label), size = 4, max.overlaps = 15) +
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
facs_fitness_ordered <- facs_fitness %>%
  mutate(mutant_residue = factor(mutant_residue, levels = rev(aa_list))) %>%
  # Explicitly state 1:378 to generate rows for entirely missing position columns
  complete(position = 1:378, mutant_residue, fill = list(mean_fit = NA)) %>%
  # Add a temporary flag indicating if the row matches the WT residue at that position
  left_join(wt_ref %>% mutate(is_wt = TRUE), by = c("position", "mutant_residue")) %>%
  # Overwrite mean_fit to 0 if it is a WT residue row, otherwise keep the existing mean_fit
  mutate(mean_fit = if_else(!is.na(is_wt) & is_wt, 0, mean_fit)) %>%
  select(-is_wt) # Clean up the temporary flag column



p_heat <- ggplot(facs_fitness_ordered,
                 aes(x = position, y = mutant_residue, fill = mean_fit)) +
  geom_tile() +
  scale_fill_gradientn(colors = rd_bu_colors, 
                       limits = c(-heat_lim, heat_lim), # Force symmetry
                       oob = scales::squish,
                       labels = function(x) sprintf("%.1f", x),
                       # If readout is "Dual_IFN", use "IFN"; otherwise use readout value
                       name = ifelse(readout == "Dual_IFN", "IFN", readout),
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
image_add <- image_read(paste0("inferred_weights_FACS_4T1_3dplot_", readout, "_additive.png")) |> 
  image_crop(geometry_area(width = 930, height = 1000))

image_add <- wrap_elements(rasterGrob(image_add)) + 
  theme(plot.margin = margin(0, 0, 0, 0, "pt"))

image_int <- image_read(paste0("inferred_weights_FACS_4T1_3dplot_", readout, "_interaction.png")) |> 
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
  paste0("inferred_aamat_additive_FACS_4T1_",readout,".csv"), heat_lim)

# Plot Interaction Matrix
p_mat_int <- plot_aa_matrix(
  paste0("inferred_aamat_interaction_FACS_4T1_",readout,".csv"), heat_lim)

# plot additive interaction histogram and fit
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
  dplyr::rename(multiplicative = interaction)

# Add wild-type residue mapping and compute fits/residuals row-by-row
facs_fit <- facs_fitness %>%
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
fitness_long <- facs_fit %>%
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
p_fit <- ggplot(facs_fit, aes(x = mean_fit, y = fit_full)) +
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
    x = "True FACS enrichment",
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
p_STING <- plot_FACS_data_fit("STING", c("a", "b", "", "c", "", "", "d", "e"))
p_IFN <- plot_FACS_data_fit("Dual_IFN", c("f", "g", "", "h", "", "", "i", "j"))
p_NFkB <- plot_FACS_data_fit("NFkB", c("k", "l", "", "m", "", "", "n", "p"))
p_ATF3 <- plot_FACS_data_fit("ATF3", c("q", "r", "", "s", "", "", "t", "u"))

# Create seequence logo
## Plot Seq Logo --------------
# 1. Load the alignment file
# Ensure the .fa file is in your current working directory
alignment <- readAAStringSet("../common_use_functions/STING_Alignment.fa")

# 2. Identify the mouse sequence
# Use grep to find the index of the mouse sequence based on unique identifier
mouse_idx <- grep("10090_0:0052ff", names(alignment))
mouse_seq <- alignment[mouse_idx]

# 3. Convert to character vector to work with string positions
mouse_char <- as.character(mouse_seq)

# 4. Find positions where the mouse sequence is NOT a gap ('-')
# This creates a logical vector: TRUE where no gap, FALSE where gap
keep_pos <- strsplit(mouse_char, "")[[1]] != "-"

# 5. Subset all sequences to those non-gap positions
# Convert full alignment to character to split it
seq_data <- as.character(alignment)
# Split strings into characters
split_seqs <- strsplit(seq_data, "")

# Subset every sequence using the keep_pos vector
# Use lapply to apply the subsetting to each sequence in the list
subset_split <- lapply(split_seqs, function(x) x[keep_pos])

# Re-collapse the sequences
final_seqs <- sapply(subset_split, paste, collapse = "")

# 6. Generate the sequence logo
# Convert to character vector for ggseqlogo
final_seqs_vec <- as.character(final_seqs)


# 3. Generate the Amino Acid Sequence Logo
color_scheme = make_col_scheme(chars=c('R','H','K','D','E',
                                       'S','T','N','Q','C',
                                       'W','F','Y','L','I','M','V',
                                       'A','G','P'),
                               
                               groups=c('Charged','Charged','Charged','Charged','Charged',
                                        'Polar','Polar','Polar','Polar','Polar',
                                        'Hydrophobic','Hydrophobic','Hydrophobic','Hydrophobic','Hydrophobic','Hydrophobic','Hydrophobic',
                                        'Special','Special','Special'), 
                               
                               cols=c('#00A99D','#00A99D','#00A99D','#00A99D','#00A99D',
                                               '#FBB03B','#FBB03B','#FBB03B','#FBB03B','#FBB03B',
                                               '#FF7BAC','#FF7BAC','#FF7BAC','#FF7BAC','#FF7BAC','#FF7BAC','#FF7BAC',
                                               '#BF8FFF','#BF8FFF','#BF8FFF'))
                                               
p_seq <- ggseqlogo(final_seqs_vec,seq_type='aa', method = 'bits', 
                   col_scheme=color_scheme) + 
  xlab("") +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    strip.text = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.title.x = element_blank(),
    legend.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"),
    legend.spacing.x = unit(5, "pt"),
    legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"),
    plot.margin = margin(t = 2, l = 2, b = 2, r = 2, unit = "pt")) +
    scale_x_continuous(limits = c(0, 379),
                     expand = c(0.01, 0.01))

p_seq_row <- (p_seq | plot_spacer()) + plot_layout(widths = c(2, 1))

p_full <- (p_seq_row / p_STING / p_IFN / p_NFkB / p_ATF3) +
  plot_layout(heights = c(0.05, 1, 1, 1, 1)) +
  
  plot_annotation(tag_levels = list(c("","a", "b", "", "c", "", "", "d", "e","",
                                      "f", "g", "", "h", "", "", "i", "j", "",
                                      "k", "l", "", "m", "", "", "n", "p", "",
                                      "q", "r", "", "s", "", "", "t", "u", ""))) &
  theme(
    plot.margin = margin(2, 2, 2, 2, "pt"),
    # Forces all subplot tags to be bold and slightly larger
    plot.tag = element_text(face = "bold", size = 20) 
  )


# Save plot
ggsave(plot=p_full, 'FACS_heatmap_weights_aamat_fit_conservation.pdf',
       width=20, height=29)

svglite("FACS_heatmap_weights_aamat_fit_conservation.svg", width=20, height=29)
print(p_full)
dev.off()


