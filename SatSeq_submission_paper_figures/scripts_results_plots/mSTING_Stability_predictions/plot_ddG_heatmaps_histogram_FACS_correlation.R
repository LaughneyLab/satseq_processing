rm(list=ls())#Clear environment
setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(RColorBrewer)
library(tidyr)
library(readxl)
library(stringr)
library(patchwork)
library(ggseqlogo)
library(Biostrings)
library(ggpmisc)
library(svglite)



plot_mutation_heatmap <- function(filename, 
                                  model_name, 
                                  pos_col, 
                                  mutant_col, 
                                  ddg_col,
                                  scale=1) {
  
  # 1. Read the data
  df <- read.csv(filename, check.names = FALSE)
  
  # 2. Clean and format columns using standard evaluation
  df_clean <- df %>%
    mutate(
      position = as.numeric(.data[[pos_col]]),
      mutant_residue = as.character(.data[[mutant_col]]),
      mean_fit = as.numeric(.data[[ddg_col]])
    ) %>%
    filter(position >= 1 & position <= 378)
  
  # 3. Define and flip the custom amino acid order
  # Reversing the order so 'W' is at the top and 'R' is at the bottom of the Y-axis
  y_order_flipped <- rev(c('R','H','K','D','E','S','T','N','Q','C', 
                           'G','P','A','V','I','L','M','F','Y','W'))
  
  df_clean$mutant_residue <- factor(df_clean$mutant_residue, levels = y_order_flipped)
  
  # convert score to ddG etc
  df_clean$mean_fit <- scale*df_clean$mean_fit
  
  # 4. Define the 95th percentile absolute value limit for symmetry
  heat_lim <- quantile(abs(df_clean$mean_fit), 0.95, na.rm = TRUE)
  
  # 5. Define Red-Blue color palette (Diverging Red-White-Blue)
  rd_bu_colors <- rev(brewer.pal(7, "RdBu"))
  
  # 6. Generate the heatmap using your exact structure
  p_heat <- ggplot(df_clean, aes(x = position, y = mutant_residue, fill = mean_fit)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = rd_bu_colors, 
      limits = c(-heat_lim, heat_lim), 
      oob = scales::squish,
      labels = function(x) sprintf("%.1f", x),
      name = model_name,
      na.value = "gray80"
    ) +
    scale_x_continuous(expand = c(0,0), breaks = seq(1, 378, by = 10)) +
    theme_cowplot() +
    labs(x = "Position", y = "AA substitution") +
    guides(fill = guide_colorbar(
      title.position = "left",                    # Places title to the side of the bar
      title.theme = element_text(angle = 90,       # Rotates the text vertically 90 degrees
                                 hjust = 0.5,      
                                 vjust = 0.5)      
    ))
  
  return(list(plot = p_heat, data = df_clean))
}

# Do it for longform data
# 1. Read your wide-form data (replace with your actual file path)
# Using check.names = FALSE to preserve exact uppercase amino acid column names
df_wide <- read.csv("ESMC/job-files_mSTING_AF_dimer_ESMC_esmc-6b-2024-12_single_chain_seq_input/llr_matrix.csv",
                    check.names = FALSE)

# 2. Pivot the wide data to long-form format
df_long <- df_wide %>%
  pivot_longer(
    cols = -c(Position, Wildtype),  # Keep Position and Wildtype, pivot all AA columns
    names_to = "Mutant AA",         # Matches the column expected by your function
    values_to = "ddG (kcal/mol)"    # Matches the column expected by your function
  )

# 3. Create a virtual text connection (tmp CSV) in memory
# This mimics writing a file and reading it back
write.csv(df_long, "tmp.csv", row.names = FALSE)

# 4. Generate the heatmap using your function and the virtual file
p_ESMC <- plot_mutation_heatmap(
  filename = "tmp.csv", 
  model_name = "ESMC", 
  pos_col = "Position", 
  mutant_col = "Mutant AA", 
  ddg_col = "ddG (kcal/mol)"
)

p_ThermoMPNN <- plot_mutation_heatmap(filename='ThermoMPNN/output_mSTING_AF_dimer_ThermoMPNN_Single_threshold_100_Distance_5_Penalize.csv', 
                                      model_name='ThermoMPNN', 
                                      pos_col='Position', 
                                      mutant_col='Mutant AA', 
                                      ddg_col='ddG (kcal/mol)', 
                                      scale = -1)

p_ProteinMPNN <- plot_mutation_heatmap(filename='ProteinMPNN-ddG/output_mSTING_AF_dimer_ProteinMPNN-ddG_model_v_48_020_rep_5_seed_123.csv', 
                                       model_name='ProteinMPNN', 
                                       pos_col='Position', 
                                       mutant_col='Post', 
                                       ddg_col='ProteinMPNN-ddG Score')

p_DeepDDG <- plot_mutation_heatmap(filename='DeepDDG/output_cols_split.csv', 
                                   model_name='DeepDDG', 
                                   pos_col='ResID', 
                                   mutant_col='Mut', 
                                   ddg_col='ddG')


#Read true STING data
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

readout <- "STING"


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
    mean_fit_true = mean(observed_fitness, na.rm = TRUE), 
    .groups = 'drop'
  ) 

# Reversing the order so 'W' is at the top and 'R' is at the bottom of the Y-axis
y_order_flipped <- rev(c('R','H','K','D','E','S','T','N','Q','C', 
                         'G','P','A','V','I','L','M','F','Y','W'))

facs_fitness$mutant_residue <- factor(facs_fitness$mutant_residue,
                                      levels = y_order_flipped)

# 4. Define the 95th percentile absolute value limit for symmetry
heat_lim <- quantile(abs(facs_fitness$mean_fit_true), 0.95, na.rm = TRUE)

# 5. Define Red-Blue color palette (Diverging Red-White-Blue)
rd_bu_colors <- rev(brewer.pal(7, "RdBu"))

# 6. Generate the heatmap using your exact structure
p_heat_true <- ggplot(facs_fitness,
                      aes(x = position,
                          y = mutant_residue,
                          fill = mean_fit_true)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rd_bu_colors, 
    limits = c(-heat_lim, heat_lim), 
    oob = scales::squish,
    labels = function(x) sprintf("%.1f", x),
    name = "STING FACS",
    na.value = "gray80"
  ) +
  scale_x_continuous(expand = c(0,0), breaks = seq(1, 378, by = 10)) +
  theme_cowplot() +
  labs(x = "Position", y = "AA substitution") +
  guides(fill = guide_colorbar(
    title.position = "left",                    # Places title to the side of the bar
    title.theme = element_text(angle = 90,       # Rotates the text vertically 90 degrees
                               hjust = 0.5,      
                               vjust = 0.5)      
  ))

# All 4 heatmaps
p_heatmaps <- p_heat_true/p_ThermoMPNN$plot/p_ProteinMPNN$plot/p_DeepDDG$plot



# Define a quick reusable function to save space
make_hist <- function(data_vec) {
  ggplot(data.frame(x = data_vec), aes(x = x)) +
    geom_histogram(bins = 150, fill = "#C0C0C0", color = "#C0C0C0") +
    labs(x = "Protein stability", y = "Frequency", title = NULL) +
    theme_cowplot() +
    xlim(-9,7) +
    geom_vline(xintercept = 0,linewidth = 0.5)
}


# Now the division operator works perfectly!
p_histograms <- make_hist(facs_fitness$mean_fit_true) /
  make_hist(p_ThermoMPNN$data$mean_fit) /
  make_hist(p_ProteinMPNN$data$mean_fit) /
  make_hist(p_DeepDDG$data$mean_fit)

# Make plots of true vs predicted weights and aa matrix
make_weights_aa_fit <- function(ddG_model) {
  
# Load model weights
weights_model <- read.csv(paste0('inferred_weights_',ddG_model,'.csv'))
colnames(weights_model) <- c('position', 'additive_model','interaction_model') 

# Load FACS true positional weights
weights_FACS <- read.csv('../FACS_position_subs_modeling/inferred_weights_FACS_4T1_STING.csv')
colnames(weights_FACS) <- c('position', 'additive_true','interaction_true') 

weights_FACS <- weights_FACS %>%
  left_join(weights_model,by='position')

# Plotting weights
p1 <- ggplot(weights_FACS,
             aes(x = additive_model,
                 y = additive_true)) +
  geom_point(alpha=0.65, size=1.5, color="gray50") +
  geom_smooth(method = "lm",linewidth = 0.75,color='black') +
  
  # 5. Multi-line Correlation showing BOTH signed R and P value
  stat_correlation(
    method = "pearson",
    label.y = "top",
    label.x = "left",
    output.type = "text", # Disables complex plotmath to prevent syntax crashes
    mapping = aes(label = paste(after_stat(r.label), 
                                after_stat(p.value.label), 
                                sep = "\n")),
    size = 4
  ) +
  geom_vline(xintercept = 0,linewidth = 0.5) +
  geom_hline(yintercept = 0,linewidth = 0.5) +
  theme_cowplot() +
  theme(legend.title = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")) +
  labs(x = "Predicted positional weight",
       y = "FACS positional weight")

# Load aa matrices
aa_FACS <- stack(read.csv('../FACS_position_subs_modeling/inferred_aamat_additive_FACS_4T1_STING.csv',
                    row.names = 1)) 

aa_model <- stack(read.csv(paste0('inferred_aamat_additive_',ddG_model,'.csv'),
                           row.names = 1))

aa_FACS$values_model <- aa_model$values

# Plotting weights
p2 <- ggplot(aa_FACS,
             aes(x = values_model,
                 y = values)) +
  geom_point(alpha=0.65, size=1.5, color="gray50") +
  geom_smooth(method = "lm",linewidth = 0.75,color='black') +
  
  # 5. Multi-line Correlation showing BOTH signed R and P value
  stat_correlation(
    method = "pearson",
    label.y = "top",
    label.x = "left",
    output.type = "text", # Disables complex plotmath to prevent syntax crashes
    mapping = aes(label = paste(after_stat(r.label), 
                                after_stat(p.value.label), 
                                sep = "\n")),
    size = 4
  ) +
  geom_vline(xintercept = 0,linewidth = 0.5) +
  geom_hline(yintercept = 0,linewidth = 0.5) +
  theme_cowplot() +
  theme(legend.title = element_blank(),
        plot.margin = margin(t = 50, r = 50, b = 50, l = 50, unit = "pt")) +
  labs(x = "Predicted substitution effect",
       y = "FACS substitution effect")

return(list(p1 = p1, p2 = p2))

}

p_ThermoMPNN <- make_weights_aa_fit('ThermoMPNN')
p_ProteinMPNN <- make_weights_aa_fit('ProteinMPNN')
p_DeepDDG <- make_weights_aa_fit('DeepDDG')


p_fit_w <-  plot_spacer() /
  p_ThermoMPNN$p1 /
  p_ProteinMPNN$p1 /
  p_DeepDDG$p1 +
  # Control heights: Spacer gets 0.5 units, each plot gets 1 unit
  plot_layout(heights = c(1.2, 1, 1, 1)) 

p_fit_aa <-  plot_spacer() /
  p_ThermoMPNN$p2 /
  p_ProteinMPNN$p2 /
  p_DeepDDG$p2 +
  # Control heights: Spacer gets 0.5 units, each plot gets 1 unit
  plot_layout(heights = c(1.2, 1, 1, 1)) 


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
  scale_x_continuous(expand = c(0,0)) +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    strip.text = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    axis.title.x = element_blank(),
    legend.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"),
    legend.spacing.x = unit(5, "pt"),
    legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0, unit = "pt"),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

p_seq_row <- (p_seq | plot_spacer()) + plot_layout(widths = c(4.65, 2.5))


# Assemble full fig
p_stability <- p_seq_row / ((p_heatmaps | p_histograms | p_fit_w | p_fit_aa) +  
  plot_layout(widths = c(7, 0.5, 1, 1))) + 
  plot_layout(heights = c(1, 20)) +
  plot_annotation(tag_levels = list(c("","a", "b", "c", "d", 
                                      "e", "f", "g", "h",
                                      "i", "j", "k"))) & 
  theme(
    plot.margin = margin(2, 2, 2, 2, "pt"),
    plot.tag = element_text(face = "bold", size = 20) 
  )

# Save plot
svglite("ddG_heatmap_histogram_fit.svg", width=25, height=14)
print(p_stability)
dev.off()




