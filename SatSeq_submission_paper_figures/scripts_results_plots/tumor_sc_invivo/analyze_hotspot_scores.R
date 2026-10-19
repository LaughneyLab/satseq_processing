rm(list=ls())#Clear environment

# 1. Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggsignif)
library(ggpubr)
library(scales)
library(patchwork)
library(tidyr)
library(readxl)
library(stringr)
library(svglite)
library(stringr)
library(RColorBrewer)

# Set working directory to where this script is
setwd(getSrcDirectory(function(){})[1]) 

# 2. Load contact matrices and calculate monomer/dimer sums
bonds_AA <- read.csv('bonds_monomer.csv', row.names = 1)
bonds_AB <- read.csv('bonds_dimer.csv', row.names = 1)
bonds_ACD <- read.csv('bonds_tetramer.csv', row.names = 1)

# 3. Load hotspot data & assign structural groups
hotspot_df <- read.csv('hotspot_module_scores_tumor_4T1.csv',row.names = 1) %>%
  mutate(
    pos_raw = case_when(
      variant_substitution == "mSTING_WT" ~ "REF",
      variant_substitution == "STOP" ~ "STOP",
      TRUE ~ str_extract(variant_substitution, "\\d+")
    ),
    pos_num = as.numeric(ifelse(pos_raw %in% c("REF", "STOP"), NA, pos_raw))
  )


# # Load per variant mets fitness data
# library(dplyr)
# library(stringr)
# library(ggplot2)
# library(ggrepel)
# library(cowplot)
# 
# # 0. Load and process data based on your specific pipeline
# hotspot_df <- read.csv('hotspot_module_scores_tumor_4T1.csv', row.names = 1) %>% 
#   mutate(
#     pos_raw = case_when(
#       variant_substitution == "mSTING_WT" ~ "REF",
#       variant_substitution == "STOP" ~ "STOP",
#       TRUE ~ str_extract(variant_substitution, "\\d+")
#     ),
#     pos_num = as.numeric(ifelse(pos_raw %in% c("REF", "STOP"), NA, pos_raw))
#   )
# 
# # 1. Calculate the baseline average for X1, X2, and X3 where pos_raw == "REF"
# ref_means <- hotspot_df %>%
#   filter(pos_raw == "REF") %>%
#   summarise(
#     ref_X1 = mean(X1, na.rm = TRUE),
#     ref_X2 = mean(X2, na.rm = TRUE),
#     ref_X3 = mean(X3, na.rm = TRUE)
#   )
# 
# ref_x1_val <- ifelse(nrow(ref_means) > 0, ref_means$ref_X1, 0)
# ref_x2_val <- ifelse(nrow(ref_means) > 0, ref_means$ref_X2, 0)
# ref_x3_val <- ifelse(nrow(ref_means) > 0, ref_means$ref_X3, 0)
# 
# # 2. Group and average per pos_num, constructing the combined WT + position label
# agg_df <- hotspot_df %>%
#   group_by(pos_num) %>%
#   summarise(
#     X1 = mean(X1, na.rm = TRUE) - ref_x1_val,
#     X2 = mean(X2, na.rm = TRUE) - ref_x2_val,
#     X3 = mean(X3, na.rm = TRUE) - ref_x3_val,
#     
#     # Extracts the leading letters (WT AA) from the first variant string in the group,
#     # then pastes it directly to the numeric position (e.g., "P" + "13" = "P13")
#     pos_label = {
#       wt_aa <- str_extract(first(variant_substitution), "^[A-Za-z]+")
#       raw_val <- first(pos_raw)
#       if (!is.na(wt_aa) && raw_val != "REF" && raw_val != "STOP") {
#         paste0(wt_aa, raw_val)
#       } else {
#         raw_val
#       }
#     },
#     .groups = "drop"
#   )
# 
# # 3. Merge with fitness_data matching position to pos_num
# merged_df <- agg_df %>%
#   left_join(fitness_data, by = c("pos_num" = "position"))
# 
# # 4. Calculate threshold and identify points with fitness > 1 SD of fitness
# fitness_mean <- mean(merged_df$additive, na.rm = TRUE)
# fitness_sd   <- sd(merged_df$additive, na.rm = TRUE)
# sd_threshold <- fitness_mean + fitness_sd
# 
# high_fitness_points <- merged_df %>%
#   filter(additive > sd_threshold)
# 
# # 5. Make the dotplot
# ggplot(merged_df, aes(x = X2, y = X3, color = additive, size = additive)) +
#   geom_point(alpha = 0.9) +
#   scale_color_gradientn(
#     colors = rev(RColorBrewer::brewer.pal(11, "RdBu")),
#     limits = c(-0.1, 0.1),
#     oob = scales::squish,
#     name = "Fitness"
#   ) +
#   scale_size_continuous(
#     range = c(0.25, 3.5), 
#     name = "Fitness"
#   ) +
#   geom_text_repel(
#     data = high_fitness_points,
#     aes(label = pos_label),
#     color = "black",
#     segment.color = "black",
#     size = 3,
#     box.padding = 0.5,
#     max.overlaps = Inf
#   ) +
#   labs(
#     x = "IFN Module score",
#     y = "NFkB Module score"
#   ) +
#   theme_cowplot()

# Load per position mets fitness
fitness_data <- read.csv('../fitness_position_subs_modeling/inferred_weights_Primary_vs_Met_4T1.csv') %>%
  mutate(mean_fitness = additive,
         position = X+1)
mets_sd <- sd(fitness_data$mean_fitness)

# 2. Build the structural bond classes dataframe
mets_bonds_variant <- data.frame(
  position = 1:378,
  monomer  = colSums(bonds_AA[, 1:378, drop = FALSE], na.rm = TRUE),
  dimer    = colSums(bonds_AB[, 1:378, drop = FALSE], na.rm = TRUE),
  tetramer = colSums(bonds_ACD[, 1:378, drop = FALSE], na.rm = TRUE)) %>%
  # Join with your processed fitness and position data
  left_join(
    select(fitness_data, position, mean_fitness),
    by = "position"
  ) %>%
  # 3. Calculate fractions and categorize groups using mean_fitness
  mutate(
    total_bonds       = monomer + dimer + tetramer,
    monomer_fraction  = ifelse(total_bonds > 0, monomer / total_bonds, 0),
    dimer_fraction    = ifelse(total_bonds > 0, dimer / total_bonds, 0),
    tetramer_fraction = ifelse(total_bonds > 0, tetramer / total_bonds, 0),
    
    # Classify based on your fractional rules and your SD cutoff
    bond_group = case_when(
      monomer_fraction > dimer_fraction & mean_fitness > mets_sd ~ "Metastatic (intra)",
      dimer_fraction >= monomer_fraction & mean_fitness > mets_sd ~ "Metastatic (inter)",
      mean_fitness <= mets_sd ~ "Non-metastatic"
    )
  )

# Join hotspot data with bond/structural data to get 'bond_group'
plot_df <- hotspot_df %>%
  #mutate(pos_raw = as.integer(pos_raw)) %>%
  left_join(select(mets_bonds_variant, position, mean_fitness,bond_group),
            by = c("pos_num" = "position")) %>%
  mutate(
    group = case_when(
      pos_raw == "REF" ~ "REF",
      pos_raw == "STOP" ~ "STOP",
      !is.na(bond_group) ~ bond_group,
      TRUE ~ "Other"
    )) %>%
  filter(group %in% c("REF", "STOP", "Metastatic (intra)",
                      "Metastatic (inter)", "Non-metastatic")) %>%
  mutate(group = factor(group, levels = c("REF", "STOP", "Metastatic (intra)",
                                          "Metastatic (inter)", "Non-metastatic"))) %>%
  pivot_longer(cols = starts_with("X"), names_to = "module", values_to = "score") %>%
  mutate(
    # Extract the number from 'X1', 'X2', etc.
    module_num = as.numeric(str_extract(module, "\\d+")),
    # Create the new clean name
    module_clean = paste0("Module ", module_num),
    # Reorder the factor levels based on the numeric order
    module_clean = factor(module_clean, levels = paste0("Module ", 1:4))
  ) %>%
  mutate(module_clean = fct_recode(
    module_clean,
    "Proliferation" = "Module 1",
    "IFN"  = "Module 2",
    "NFkB" = "Module 3"
  )) %>% 
  
  # 3. Filter down to your 3 specific analysis groups and rename them
  filter(group %in% c("Metastatic (inter)",
                      "Metastatic (intra)", "Non-metastatic")) %>% 
  mutate(
    # Explicitly set the x-axis display order for your plot
    group = fct_relevel(group, "Metastatic (inter)",
                        "Metastatic (intra)",
                        "Non-metastatic"))


# 4. Define specific Wilcoxon comparisons
my_comparisons <- list(
  c("Metastatic (inter)", "Metastatic (intra)"),
  c("Metastatic (inter)", "Non-metastatic"),
  c("Metastatic (intra)", "Non-metastatic")
)

# 5. Plot the Violin/Boxplot with Significance
p1 <- ggplot(plot_df, aes(x = group, y = score, color = group)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) + 
  #stat_summary(fun = mean, geom = "crossbar", width = 0.3,
  #color = "black", linewidth = 0.3, show.legend = FALSE) +
  geom_hline(yintercept = 0, color="black", linetype="dashed", linewidth=0.4) +
  facet_wrap(~module_clean, scales = "free_y", ncol = 4) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format", 
    bracket.size = 0.5,
    size = 6) +
  scale_color_manual(
    values = c("Metastatic (intra)" = "#5e3c99",
               "Metastatic (inter)" = "#d01c8b",
               "Non-metastatic" = "black")
  ) +
  theme_cowplot() +
  labs(x = NULL, y = "Module Score", fill = "Variant Group") +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "none",
    strip.text = element_text(margin = margin(t = 1, b = 1))
  )

# Save plot
svglite("Hotspot_score_mets_bonds.svg", width = 10, height = 6)
print(p1)
dev.off()

# Plot sc 2d density contours
# Step A: Pair IFN and NFkB scores for each individual cell using a row index
plot_wide <- plot_df %>%
  filter(module_clean %in% c("NFkB", "IFN")) %>%
  group_by(variant_substitution, pos_raw, pos_num, mean_fitness, bond_group, group, module_clean) %>%
  mutate(cell_index = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c(variant_substitution, pos_raw, pos_num, mean_fitness, bond_group, group, cell_index),
    names_from = module_clean,
    values_from = score
  ) %>%
  drop_na(NFkB, IFN) %>%
  mutate(label = paste0(pos_raw, pos_num, variant_substitution))


# Define your base color mapping
group_colors <- c(
  "Metastatic (intra)" = "#5e3c99",
  "Metastatic (inter)" = "#d01c8b",
  "Non-metastatic"     = "black"
)

group_colors_gray <- c(
  "Metastatic (intra)" = "black",
  "Metastatic (inter)" = "black",
  "Non-metastatic"     = "black"
)

p4 <- ggplot(plot_wide,aes(x = IFN, y = NFkB)) +
  # 1. Filled 2D density contour bands (Gradients from white to group color)
  geom_density_2d_filled(
    aes(fill = group,alpha = after_stat(level)), # Map fill to group, alpha to density depth
    bins = 8,
    show.legend = F
  ) +
  # Control the alpha range so the lowest density level is white/transparent
  scale_alpha_manual(values = seq(0, 0.85, length.out = 8)) + 
  
  # 2. Extreme outer contour lines matching the group color
  geom_density_2d(
    aes(color = group),
    bins = 8, 
    linewidth = 0.5, 
    show.legend = FALSE
  ) +
  
  # Intersecting baseline grids at 0
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) + 
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) + 
  
  # Separate the plot into 3 distinct panels based on group
  facet_wrap(~group, ncol = 3) +
  
  # Apply your specific group colors to fill and lines
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  
  # Labels and themes
  labs(
    x = "IFN Module score", 
    y = "NFkB Module score"
  )  +
  theme_cowplot(font_size = 20) +
  xlim(-15, 15) +
  ylim(-8, 8) +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(size = 20),
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75)
  )


# Calculate Earth mover's distance
library(emdist)

# 1. Clean data and split into groups
clean_data <- plot_wide[!is.na(plot_wide$IFN) & !is.na(plot_wide$NFkB), ]
groups <- split(clean_data, clean_data$group)

# 2. Ultra-fast 1D histogram builder 
get_1d_bins <- function(values, limits, n_bins = 100) {
  v_clamped <- pmax(limits[1], pmin(limits[2], values))
  bin_width <- (limits[2] - limits[1]) / n_bins
  bin_idx <- floor((v_clamped - limits[1]) / bin_width) + 1
  bin_idx <- pmin(bin_idx, n_bins)
  
  counts <- tabulate(bin_idx, nbins = n_bins)
  weights <- counts / sum(counts)
  bin_centers <- seq(limits[1] + bin_width/2, limits[2] - bin_width/2, length.out = n_bins)
  return(matrix(c(weights, bin_centers), ncol = 2))
}

# Generate 1D signatures using your exact coordinate boundaries
ifn_intra <- get_1d_bins(groups[["Metastatic (intra)"]]$IFN, limits = c(-15, 15))
ifn_inter <- get_1d_bins(groups[["Metastatic (inter)"]]$IFN, limits = c(-15, 15))
ifn_non   <- get_1d_bins(groups[["Non-metastatic"]]$IFN,     limits = c(-15, 15))

nfkb_intra <- get_1d_bins(groups[["Metastatic (intra)"]]$NFkB, limits = c(-8, 8))
nfkb_inter <- get_1d_bins(groups[["Metastatic (inter)"]]$NFkB, limits = c(-8, 8))
nfkb_non   <- get_1d_bins(groups[["Non-metastatic"]]$NFkB,     limits = c(-8, 8))


# 4. Helper function to compute EMD and inject direction sign
# If Group A has a HIGHER mean than Group B, the result is POSITIVE (+)
get_signed_emd <- function(sig_former, sig_latter, vals_former, vals_latter) {
  magnitude <- emdist::emd(sig_former, sig_latter)
  direction <- sign(mean(vals_latter) - mean(vals_former)) # Latter - Former
  return(magnitude * direction)
}

# Compute signed values
ifn_comp1 <- get_signed_emd(ifn_non,   ifn_inter,   groups[["Non-metastatic"]]$IFN,     groups[["Metastatic (inter)"]]$IFN)
ifn_comp2 <- get_signed_emd(ifn_non,   ifn_intra,   groups[["Non-metastatic"]]$IFN,     groups[["Metastatic (intra)"]]$IFN)
ifn_comp3 <- get_signed_emd(ifn_intra, ifn_inter,   groups[["Metastatic (intra)"]]$IFN, groups[["Metastatic (inter)"]]$IFN)

nfkb_comp1 <- get_signed_emd(nfkb_non,   nfkb_inter,   groups[["Non-metastatic"]]$NFkB,     groups[["Metastatic (inter)"]]$NFkB)
nfkb_comp2 <- get_signed_emd(nfkb_non,   nfkb_intra,   groups[["Non-metastatic"]]$NFkB,     groups[["Metastatic (intra)"]]$NFkB)
nfkb_comp3 <- get_signed_emd(nfkb_intra, nfkb_inter,   groups[["Metastatic (intra)"]]$NFkB, groups[["Metastatic (inter)"]]$NFkB)

# ==============================================================================
# 4. TIDY DATA STRUCTURE WITH FORCED FACTOR ORDER
# ==============================================================================
comparison_order <- c(
  "Non-metastatic to Metastatic (inter)",
  "Non-metastatic to Metastatic (intra)",
  "Metastatic (intra) to Metastatic (inter)"
)

# Extract the destination group for each comparison
get_destination <- function(x) {
  sapply(strsplit(x, " to "), `[`, 2)
}

# Build data frame and scale alpha based on bar height (using absolute values for fairness)
plot_df_emd <- data.frame(
  Comparison = factor(rep(comparison_order, times = 2), levels = rev(comparison_order)),
  Variable = rep(c("IFN", "NFkB"), each = 3),
  EMD = c(ifn_comp1, ifn_comp2, ifn_comp3, nfkb_comp1, nfkb_comp2, nfkb_comp3)
)

# Extract destination group from the comparison string
plot_df_emd$Destination <- sapply(strsplit(as.character(plot_df_emd$Comparison), " to "), `[`, 2)

# Calculate absolute heights for alpha scaling (0 to 1 range)
max_emd <- max(abs(plot_df_emd$EMD), na.rm = TRUE)
plot_df_emd$BarAlpha <- abs(plot_df_emd$EMD) / max_emd

# ==============================================================================
# 5. HORIZONTAL BAR PLOT GENERATION WITH LINE TYPE & ALPHA
# ==============================================================================
p_emd <- ggplot(plot_df_emd, aes(
  x = Comparison, 
  y = EMD, 
  fill = Destination, 
  alpha = BarAlpha,     # Re-added the proportional alpha mapping
  linetype = Variable
)) +
  # Black borders ensure the solid and dashed patterns remain visible over translucent fills
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  coord_flip() +
  scale_fill_manual(values = group_colors) +
  scale_alpha_identity() + # Maps exact calculated BarAlpha values directly
  scale_linetype_manual(values = c("IFN" = "solid", "NFkB" = "dashed")) +
  labs(
    x = NULL, 
    y = "Signed earth mover's distance", 
    linetype = "Module"
  ) +
  theme_cowplot(font_size = 20) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  ) +
  # Keeping only the Module legend visible
  guides(
    fill = "none",
    alpha = "none",
    linetype = guide_legend(override.aes = list(fill = "white", color = "black"))
  )

p_full <- (p4 | p_emd) + plot_layout(widths = c(2.8, 1))

svglite("Hotspot_score_2d_density_mets_bonds.svg", width = 18, height = 5)
print(p_full)
dev.off()

saveRDS(p_full,"Hotspot_score_2d_density_mets_bonds.rds")

# Extract variants in NFkB mid IFN lo population

# 2. Count cells at each position and fill missing positions (1:378) with 0
plot_counts <- plot_wide %>%
  filter(group %in% c("Metastatic (inter)",
         "Metastatic (intra)",
         "Non-metastatic")) %>%
  group_by(group, pos_num) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  complete(group, pos_num = 1:378, fill = list(cell_count = 0))

# 3. Plot as a single column of stacked bar charts
p_bar <- ggplot(plot_counts[plot_counts$group %in% c("Metastatic (inter)",
                                                     "Metastatic (intra)",
                                                     "Non-metastatic"),], 
                aes(x = pos_num, y = cell_count)) +
  geom_col(fill = "black", width = 1) +
  scale_x_continuous(limits = c(0, 379),
                     breaks = seq(5, 378, by = 10),
                     expand = c(0, 0)) +
  facet_wrap(~group, ncol = 1, scales = "free_y",  strip.position = "right") +
  labs(
    x = "Position",
    y = "Cell count in tumor") +
  theme_cowplot() +
  theme(
    strip.text = element_text(size = 10.5),
    panel.spacing = unit(0.8, "lines")
  )

# Add LOF marks
# 1. Load LOF data
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
                     breaks = seq(5, 378, by = 10),
                     expand = c(0, 0)) +
  theme_cowplot() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(), # Hide X text since the bottom plot covers it
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )

# Combine plots
p_full <- p_tracks / p_bar + plot_layout(heights = c(1, 8))

svglite("tumor_cell_count_by_mets_bond.svg", width = 15, height = 5)
print(p_full)
dev.off()

# Combine dimer and monomer
plot_df <- plot_df %>%
  mutate(group = fct_collapse(group,
                              "Metastatic" = c("Metastatic (inter)",
                                               "Metastatic (intra)"),
                              "Non-metastatic" = "Non-metastatic"))

my_comparisons <- list(c("Metastatic", "Non-metastatic"))

# 5. Plot the Violin/Boxplot with Significance
p2 <- ggplot(plot_df, aes(x = group, y = score, color = group)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) + 
  #stat_summary(fun = mean, geom = "crossbar", width = 0.3,
  #color = "black", linewidth = 0.3, show.legend = FALSE) +
  geom_hline(yintercept = 0, color="black", linetype="dashed", linewidth=0.4) +
  facet_wrap(~module_clean, scales = "free_y", ncol = 4) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format", 
    bracket.size = 0.5,
    size = 6) +
  scale_color_manual(
    values = c("Metastatic" = "#d73027",
               "Non-metastatic" = "black")
  ) +
  theme_cowplot() +
  labs(x = NULL, y = "Module Score", fill = "Variant Group") +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "none",
    strip.text = element_text(margin = margin(t = 1, b = 1))
  )

# Save plot
svglite("Hotspot_score_mets.svg", width = 10, height = 6)
print(p2)
dev.off()

# Plot sc 2d density contours
# Step A: Pair IFN and NFkB scores for each individual cell using a row index
plot_wide <- plot_df %>%
  filter(module_clean %in% c("NFkB", "IFN")) %>%
  group_by(variant_substitution, pos_raw, pos_num, mean_fitness, bond_group, group, module_clean) %>%
  mutate(cell_index = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c(variant_substitution, pos_raw, pos_num, mean_fitness, bond_group, group, cell_index),
    names_from = module_clean,
    values_from = score
  ) %>%
  drop_na(NFkB, IFN) %>%
  mutate(label = paste0(pos_raw, pos_num, variant_substitution))

# Step C: Plot individual cells with filled contours, faceted by group
p4 <- ggplot(plot_wide, aes(x = IFN, y = NFkB)) +
  # 1. Filled 2D density contour bands (Calculates density topography automatically)
  geom_density_2d_filled(
    bins = 12,
    alpha = 0.65,
    show.legend = FALSE # Hides the default multi-shade gradient legend
  ) +
  # Intersecting baseline grids at 0
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) + 
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) + 
  # Separate the plot into 3 distinct panels based on group
  facet_wrap(~group, ncol = 3) +
  # Use a standard clean color palette for the density steps
  scale_fill_manual(
    values = colorRampPalette(c("white", "black"))(12)
  ) +
  # Labels and themes
  labs(
    x = "IFN Module score", 
    y = "NFkB Module score")  +
  theme_cowplot() +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
  )

svglite("Hotspot_score_2d_density_mets.svg", width = 8, height = 4)
print(p4)
dev.off()

# Do module gsea dotplot

# --- 1. Define Variables ---
file_path <- "hotspot_modules_gsea_tumor.xlsx"

# --- 2. Load Data from Excel Sheets ---
df_pval <- read_excel(file_path, sheet = "Adjusted P-value")
df_over <- read_excel(file_path, sheet = "Overlap")

# --- 3. Reshape and Filter P-values ---
df_pval_long <- df_pval %>%
  pivot_longer(cols = -Term, names_to = "Module", values_to = "P_val") %>%
  mutate(Module = factor(Module, levels = unique(Module)))

# Identify terms where at least one module meets the significance threshold
keep_terms <- df_pval_long %>%
  filter(P_val < 1e-5) %>%
  pull(Term) %>%
  unique()

# --- 4. Reshape Overlap Matrix and Calculate Fractions ---
df_over_long <- df_over %>%
  pivot_longer(cols = -Term, names_to = "Module", values_to = "Overlap_Str") %>%
  separate(Overlap_Str, into = c("Overlap_Genes", "Total_Genes"), sep = "/", convert = TRUE) %>%
  mutate(
    Fraction = ifelse(Total_Genes > 0, Overlap_Genes / Total_Genes, 0)
  )

# --- 5. Merge Datasets and Apply Filter ---
df_combined <- df_pval_long %>%
  inner_join(df_over_long, by = c("Term", "Module")) %>%
  filter(Term %in% keep_terms)

# --- 6. Perform Hierarchical Clustering for Term Order ---
mat_clustering <- df_combined %>%
  mutate(NegLogP = -log10(P_val)) %>% 
  select(Term, Module, NegLogP) %>%
  pivot_wider(names_from = Module, values_from = NegLogP, values_fill = 0) %>%
  column_to_rownames("Term")

hr <- hclust(dist(mat_clustering), method = "complete")
pathway_order <- rownames(mat_clustering)[hr$order]

# --- 7. Clean Text Names and Apply Factor Order ---
df_plot <- df_combined %>%
  mutate(
    Term_Clean = Term %>% 
      str_remove_all("^HALLMARK_|^GOBP_") %>% 
      str_replace_all("_", " ") %>% 
      str_to_title() %>% 
      str_replace("Tnfa Signaling Via NfkB", "TNFa Signaling via NFkB") %>% 
      str_replace("Il2 Stat5", "IL2 STAT5") %>% 
      str_replace("Dna", "DNA") %>% 
      str_replace("Tgf Beta", "TGFb") %>% 
      str_replace("E2f", "E2F") %>% 
      str_replace("G2m", "G2M") %>% 
      str_replace("Uv", "UV") %>% 
      str_replace("NfKb", "NFkB") %>% 
      str_replace("Il6 Jak Stat3", "IL6 JAK STAT3") %>% 
      str_replace("Reactive Oxygen Species", "ROS"),
    Term = factor(Term, levels = pathway_order)
  )


df_plot$Module <- paste0("Module ", df_plot$Module)

p_final <- ggplot(df_plot, aes(x = -log10(P_val), y = Term)) +
  # Create solid black bars with uniform heights
  geom_col(fill = "black", color = "white", linewidth = 0.2, width = 0.8) + 
  # Apply the cleaned up labels directly to the shared y-axis
  scale_y_discrete(labels = setNames(df_plot$Term_Clean, df_plot$Term)) +
  # Force layout into a single row with a globally shared y-axis
  facet_wrap(~Module, nrow = 1, scales = "fixed") + 
  theme_cowplot() +
  labs(
    x = "-log10(P value)", 
    y = "Hallmark Pathway"
  ) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    # Style the module facet banner labels
    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "gray85", colour = "black", linewidth = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10))
  )

svglite("Hotspot_module_enrichR_barplot.svg", width = 6, height = 4)
print(p_final)
dev.off()
