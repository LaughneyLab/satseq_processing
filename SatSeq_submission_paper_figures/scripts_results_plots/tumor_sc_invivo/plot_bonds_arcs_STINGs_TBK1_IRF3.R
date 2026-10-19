rm(list=ls()) # Clear environment
setwd(getSrcDirectory(function(){})) # Set working directory to script location

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(cowplot)
library(svglite)
library(patchwork)

num_res_STING <- 378
num_res_TBK1 <- 729
num_res_IRF3 <- 419

# ==============================================================================
# 1. LOAD METADATA AND LOF PATHWAY TRACKS
# ==============================================================================
lof_data     <- read.csv("../cGAMP_sc/LOF_variants.csv")
fitness_data <- read.csv("../cGAMP_sc/FACS_weights_all.csv")
mets_data    <- read.csv('../fitness_position_subs_modeling/inferred_weights_Primary_vs_Met_4T1.csv')

# Calculate the metric threshold limit
mets_cutoff  <- sd(mets_data$additive)

# Structural dataframes for bars
df_bars  <- data.frame(Residue = 1:378, Fitness = mets_data$additive)
df_sting <- data.frame(Residue = 1:378, Fitness = fitness_data$STING_expression[1:378])

# Build the discrete pixel matrix tracks for the very top of the system
df_tracks <- data.frame(residue = 1:378) %>% 
  left_join(lof_data, by = c("residue" = "residue")) %>% 
  mutate(
    `IFN LOF`  = ifelse(grepl("IFNa", group), "TRUE", "FALSE"),
    `NFkB LOF` = ifelse(grepl("NFkB", group), "TRUE", "FALSE"),
    `UPR LOF`  = ifelse(grepl("UPR", group), "TRUE", "FALSE")
  ) %>% 
  mutate(across(c(`IFN LOF`, `NFkB LOF`, `UPR LOF`), ~replace_na(., "FALSE"))) %>% 
  pivot_longer(cols = c(`IFN LOF`, `NFkB LOF`, `UPR LOF`), names_to = "Pathway", values_to = "Status")

# Top Panel Plot: Pixel matrices (aligned layout)
p_tracks <- ggplot(df_tracks, aes(x = residue, y = Pathway, fill = interaction(Pathway, Status))) + 
  geom_tile(color = "white", linewidth = 0) + 
  scale_fill_manual(
    values = c(
      "UPR LOF.TRUE"   = "#35b6e3", "UPR LOF.FALSE"  = "white", 
      "IFN LOF.TRUE"   = "#e3312a", "IFN LOF.FALSE"  = "white", 
      "NFkB LOF.TRUE"  = "#7ac943", "NFkB LOF.FALSE" = "white"
    ), 
    guide = "none"
  ) + 
  scale_x_continuous(limits = c(0, 379), breaks = seq(5, 378, by = 10), expand = c(0, 0)) + 
  theme_cowplot() + 
  theme(
    axis.title   = element_blank(), 
    axis.text.x  = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(), 
    panel.grid   = element_blank(), 
    axis.line.y  = element_blank(), 
    axis.line.x  = element_blank(),
    plot.margin  = margin(t = 5, r = 10, b = 0, l = 10, unit = "pt")
  )

# ==============================================================================
# 2. INITIALIZE INTERACTION MATRICES
# ==============================================================================
mat_monomer  <- matrix(0, nrow = num_res_STING, ncol = num_res_STING)
mat_dimer    <- matrix(0, nrow = num_res_STING, ncol = num_res_STING)
mat_tetramer <- matrix(0, nrow = num_res_STING, ncol = num_res_STING)
mat_tbk1     <- matrix(0, nrow = num_res_STING, ncol = num_res_TBK1)
mat_irf3     <- matrix(0, nrow = num_res_STING, ncol = num_res_IRF3)

# ==============================================================================
# 3. FUNCTION: PROCESS INTERACTION FILE WITH METS CONDITION FILTERS
# ==============================================================================
process_interaction_file <- function(file_path) {
  lines <- readLines(file_path)
  
  # Define chain groups
  sting_chains <- c("C", "D", "E", "F")
  irf3_chains <- c("G", "H")
  tbk1_chains <- c("A", "B")
  
  for (line in lines) {
    parts <- unlist(strsplit(trimws(line), "\\s+"))
    slash_indices <- grep("^/", parts)
    if (length(slash_indices) < 2) next
    
    # Extract chains and residue values
    chain1 <- gsub("/", "", parts[slash_indices[1]])
    res1 <- as.numeric(parts[slash_indices[1] + 2])
    chain2 <- gsub("/", "", parts[slash_indices[2]])
    res2 <- as.numeric(parts[slash_indices[2] + 2])
    
    if (is.na(res1) || is.na(res2) || res1 < 1 || res2 < 1) next
    
    # ---------------------------------------------------------
    # 1. PROCESS IRF3 TO STING INTERACTIONS (ALL STING CHAINS)
    # ---------------------------------------------------------
    is_irf3_interaction <- (chain1 %in% sting_chains && chain2 %in% irf3_chains) || 
      (chain2 %in% sting_chains && chain1 %in% irf3_chains)
    
    if (is_irf3_interaction) {
      res_sting <- if (chain1 %in% sting_chains) res1 else res2
      res_irf3  <- if (chain1 %in% irf3_chains) res1 else res2
      
      mat_irf3[res_sting, res_irf3] <<- mat_irf3[res_sting, res_irf3] + 1
      next 
    }
    
    # ---------------------------------------------------------
    # 2. PROCESS TBK1 TO STING INTERACTIONS (ALL STING CHAINS)
    # ---------------------------------------------------------
    is_tbk1_interaction <- (chain1 %in% sting_chains && chain2 %in% tbk1_chains) || 
      (chain2 %in% sting_chains && chain1 %in% tbk1_chains)
    
    if (is_tbk1_interaction) {
      res_sting <- if (chain1 %in% sting_chains) res1 else res2
      res_tbk1  <- if (chain1 %in% tbk1_chains) res1 else res2
      
      mat_tbk1[res_sting, res_tbk1] <<- mat_tbk1[res_sting, res_tbk1] + 1
      next 
    }
    
    # ---------------------------------------------------------
    # 3. PROCESS INTRA-STING (STRICTLY VIA CHAIN C)
    # ---------------------------------------------------------
    if (chain1 != "C" && chain2 != "C") next
    
    res_C <- if (chain1 == "C") res1 else res2
    partner_ch <- if (chain1 == "C") chain2 else chain1
    res_partner <- if (chain1 == "C") res2 else res1
    
    if (partner_ch == "C") {
      mat_monomer[res_C, res_partner] <<- mat_monomer[res_C, res_partner] + 1
      mat_monomer[res_partner, res_C] <<- mat_monomer[res_partner, res_C] + 1
    } else if (partner_ch == "D") {
      mat_dimer[res_C, res_partner] <<- mat_dimer[res_C, res_partner] + 1
      mat_dimer[res_partner, res_C] <<- mat_dimer[res_partner, res_C] + 1
    } else if (partner_ch %in% c("E", "F")) {
      mat_tetramer[res_C, res_partner] <<- mat_tetramer[res_C, res_partner] + 1
      mat_tetramer[res_partner, res_C] <<- mat_tetramer[res_partner, res_C] + 1
    }
  }
}

process_interaction_file("tbk1_sting_irf3_full_data_0_bonds.txt")
process_interaction_file("tbk1_sting_irf3_full_data_0_hbonds.txt")

# ==============================================================================
# 4. GENERATE DESIGNATED GEOMETRY DATAFRAMES
# ==============================================================================
generate_arc_dataframe <- function(contact_matrix, df_bars, mets_cutoff, check_triangle = FALSE) {
  if (check_triangle) {
    pairs <- which(contact_matrix > 0 & lower.tri(contact_matrix), arr.ind = TRUE)
  } else {
    pairs <- which(contact_matrix > 0, arr.ind = TRUE)
  }
  
  df_list <- list()
  if (nrow(pairs) == 0) return(data.frame())
  
  for (k in 1:nrow(pairs)) {
    r1 <- pairs[k, 1]
    r2 <- pairs[k, 2]
    bonds <- contact_matrix[r1, r2]
    
    # CHANGE: Only generate arc segment if at least one residue passes the threshold
    if (df_bars$Fitness[r1] >= mets_cutoff || df_bars$Fitness[r2] >= mets_cutoff) {
      span <- abs(r1 - r2)
      midpoint <- (r1 + r2) / 2
      dip_depth <- -0.01 - (span * 0.003)
      arc_id <- paste0("arc_", r1, "_", r2, "_", k)
      
      df_list[[k]] <- data.frame(
        x = c(r1, midpoint, r2),
        y = c(0, dip_depth, 0),
        Bonds = bonds,
        Group = arc_id
      )
    }
  }
  return(do.call(rbind, df_list))
}

# Negative offsets drop tracking loops directly under the shared X-axis baseline
# CHANGE: Pass fitness thresholds to filter arcs while keeping matrices complete
df_arcs_monomer  <- generate_arc_dataframe(mat_monomer, df_bars, mets_cutoff, check_triangle = TRUE)
df_arcs_dimer    <- generate_arc_dataframe(mat_dimer, df_bars, mets_cutoff, check_triangle = TRUE)
df_arcs_tetramer <- generate_arc_dataframe(mat_tetramer, df_bars, mets_cutoff, check_triangle = TRUE)

generate_converging_curves <- function(contact_matrix, target_x_center, target_y) {
  c_counts <- rowSums(contact_matrix)
  active_res <- which(c_counts > 0)
  
  if (length(active_res) == 0) return(data.frame())
  
  df_list <- list()
  for (k in seq_along(active_res)) {
    r_C <- active_res[k]
    bonds <- c_counts[r_C]
    curve_id <- paste0("converge_", r_C, "_", target_y)
    
    df_list[[k]] <- data.frame(
      x = c(r_C, r_C, target_x_center, target_x_center),
      y = c(0, target_y * 0.3, target_y * 0.7, target_y),
      Bonds = bonds,
      Group = curve_id
    )
  }
  return(do.call(rbind, df_list))
}

df_curves_tbk1 <- generate_converging_curves(mat_tbk1,
                                             target_x_center = 370,
                                             target_y = -0.25)

df_curves_irf3 <- generate_converging_curves(mat_irf3,
                                             target_x_center = 200,
                                             target_y = -0.25)

# ==============================================================================
# 5. BOTTOM PANEL: INTERACTION TRACKS ALIGNED ON THE BASAL X LAYOUT
# ==============================================================================
p_bottom <- ggplot() +
  geom_hline(yintercept = 0, color = "black") +
  annotate("text", x = 370, y = -0.28, label = "TBK1",
           color = "#377eb8", fontface = "bold", size = 5) +
  annotate("text", x = 200, y = -0.28, label = "IRF3",
           color = "#FF7F00", fontface = "bold", size = 5)

if (!is.null(df_arcs_monomer) && nrow(df_arcs_monomer) > 0) {
  p_bottom <- p_bottom + geom_bezier(data = df_arcs_monomer, aes(x = x, y = y, group = Group, linewidth = Bonds, color = "Monomer"), alpha = 0.5)
}
if (!is.null(df_arcs_dimer) && nrow(df_arcs_dimer) > 0) {
  p_bottom <- p_bottom + geom_bezier(data = df_arcs_dimer, aes(x = x, y = y, group = Group, linewidth = Bonds, color = "Dimer"), alpha = 0.5)
}
if (!is.null(df_arcs_tetramer) && nrow(df_arcs_tetramer) > 0) {
  p_bottom <- p_bottom + geom_bezier(data = df_arcs_tetramer, aes(x = x, y = y, group = Group, linewidth = Bonds, color = "Tetramer"), alpha = 0.8)
}
if (!is.null(df_curves_tbk1) && nrow(df_curves_tbk1) > 0) {
  p_bottom <- p_bottom + geom_bezier(data = df_curves_tbk1, aes(x = x, y = y, group = Group, linewidth = Bonds, color = "TBK1"), alpha = 0.5)
}
if (!is.null(df_curves_irf3) && nrow(df_curves_irf3) > 0) {
  p_bottom <- p_bottom + geom_bezier(data = df_curves_irf3, aes(x = x, y = y, group = Group, linewidth = Bonds, color = "IRF3"), alpha = 0.5)
}

p_bottom <- p_bottom + 
  scale_color_manual(values = c("Monomer" = "#5e3c99",
                                "Dimer" = "#d01c8b", 
                                "Tetramer" = "#018571",
                                "TBK1" = "#377eb8", 
                                "IRF3" = "#ff7f00"),
                     breaks = c("Monomer", "Dimer", "Tetramer", "TBK1", "IRF3"),
                     name = "Bond type") + 
  scale_linewidth(range = c(0.3, 1), name = "# Bonds", guide="none") + 
  scale_x_continuous(limits = c(0, 379),
                     breaks = seq(5, 378, by = 10),
                     expand = c(0, 0)) +
  theme_cowplot(font_size = 15) +
  theme(legend.position   = "right",
        plot.margin = margin(t = 0, r = 10, b = 10, l = 10, unit = "pt"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y      = element_blank()) +
  labs(x = "Residue")

p_middle <- ggplot() +geom_col(data = df_bars,
                               aes(x = Residue,
                                   y = Fitness, fill = Fitness), width = 0.8) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       limits = c(-0.1, 0.1), 
                       oob = scales::squish, name = "Fitness") +
  scale_x_continuous(limits = c(0, 379), breaks = seq(5, 378, by = 10),
                     expand = c(0, 0)) +
  theme_cowplot(font_size = 15) +
  theme(legend.position = "right",
        axis.text.x     = element_blank(),
        axis.title.x    = element_blank(),
        plot.margin     = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt")) +
  labs(y = "Metastatic fitness")

# ==============================================================================
# NEW: GENERATE PROPORTIONAL STACKED BAR PLOT FOR STING INTERACTION SHARE
# ==============================================================================
# Aggregate total bonds per residue position straight from the processed matrices
df_share <- data.frame(
  position  = rep(1:378, each = 3),
  type      = rep(c("Monomer", "Dimer", "Tetramer"), times = 378),
  nbonds    = c(rbind(rowSums(mat_monomer), rowSums(mat_dimer), rowSums(mat_tetramer)))
) %>%
  filter(nbonds > 0) %>%
  group_by(position, type) %>%
  summarise(total_bonds = sum(nbonds), .groups = "drop") %>%
  mutate(type = factor(type, levels = c("Monomer", "Dimer", "Tetramer")))

# Build the fractional percentage stacked bar plot
p_share_bar <- ggplot(df_share, aes(x = position, y = total_bonds, fill = type)) +
  geom_col(position = "fill", linewidth = 0.1, color = "white") +
  scale_fill_manual(
    values = c("Monomer" = "#5e3c99", "Dimer" = "#d01c8b", "Tetramer" = "#018571")
  ) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 379), breaks = seq(5, 378, by = 10), expand = c(0, 0)) +
  labs(
    y = "% Bonds"
  ) +
  theme_cowplot() +
  theme(
    axis.text.x      = element_blank(),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(size = 15),
    legend.position  = "none",
    plot.margin      = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt"),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    axis.line        = element_blank(),
    axis.ticks.x     = element_blank()
  )

# ==============================================================================
# UPDATE: RE-ARRANGE THE COMPOSITE GRID LAYOUT
# ==============================================================================
# Add the new percentage plot right between the fitness barplot and the bottom arcs
p_final <- p_tracks / p_middle / p_bottom / p_share_bar + 
  plot_layout(heights = c(0.10, 0.25, 0.40, 0.2))

ggsave(filename = "mets_LOF_bonds_complex.svg",
       plot = p_final, device = "svg", width = 18, height = 5,
       limitsize = FALSE)
