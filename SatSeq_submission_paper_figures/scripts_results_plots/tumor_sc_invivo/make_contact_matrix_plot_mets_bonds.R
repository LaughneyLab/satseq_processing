rm(list=ls())#Clear environment

library(stringr)
library(igraph)
library(ggplot2)
library(ggforce)
library(scales)
library(dplyr)
library(tidyverse)
library(cowplot)
library(patchwork)
library(readr)
library(purrr)
library(ggpubr)
library(rstatix)
library(svglite)

setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is


# Process dimer File to get contact matrix------------

# --- Parameters ---
file_contacts <- "mSTING_AF_dimer_bonds.txt"
file_hbonds   <- "mSTING_AF_dimer_hbonds.txt"
num_res       <- 378

# Initialize empty 378x378 matrices
mat_A_A <- matrix(0, nrow = num_res, ncol = num_res)
mat_A_B <- matrix(0, nrow = num_res, ncol = num_res)

# ==============================================================================
# FUNCTION: Process dimer File
# ==============================================================================
process_dimer_file <- function(file_path) {
  lines <- readLines(file_path)
  
  # Process line by line
  for (line in lines) {
    # Split by whitespace
    parts <- unlist(strsplit(trimws(line), "\\s+"))
    
    # Find positions of all chain identifiers starting with '/'
    slash_indices <- grep("^/", parts)
    
    # We need at least two chain references to define an interaction pair
    if (length(slash_indices) < 2) next
    
    # Extract structural positions (using first two distinct chain references)
    chain1 <- gsub("/", "", parts[slash_indices[1]])
    res1   <- as.numeric(parts[slash_indices[1] + 2])
    
    chain2 <- gsub("/", "", parts[slash_indices[2]])
    res2   <- as.numeric(parts[slash_indices[2] + 2])
    
    # Safety validation check
    if (is.na(res1) || is.na(res2) || res1 > num_res || res2 > num_res || res1 < 1 || res2 < 1) next
    
    # Populate Global Matrices (Superassignment inside function)
    if (chain1 == "A" && chain2 == "A") {
      mat_A_A[res1, res2] <<- mat_A_A[res1, res2] + 1
      mat_A_A[res2, res1] <<- mat_A_A[res2, res1] + 1
      
    } else if ((chain1 == "A" && chain2 == "B") || (chain1 == "B" && chain2 == "A")) {
      res_A <- if (chain1 == "A") res1 else res2
      res_B <- if (chain2 == "B") res2 else res1
      
      mat_A_B[res_A, res_B] <<- mat_A_B[res_A, res_B] + 1
      mat_A_B[res_B, res_A] <<- mat_A_B[res_B, res_A] + 1
    }
  }
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

# 1. Parse contacts
process_dimer_file(file_contacts)

# 2. Parse H-bonds and add them to the existing matrices
process_dimer_file(file_hbonds)

# Process tetramer File to get contact matrix------------

# --- Parameters ---
file_contacts <- "mSTING_tetramer_bonds_8flk.txt"
file_hbonds   <- "mSTING_tetramer_hbonds_8flk.txt"
num_res       <- 378

# --- Define all possible unique chain combinations ---
chains <- c("A", "B", "C", "D")
matrix_names <- c()

# Generate the 10 distinct matrix names (4 intra-chain + 6 inter-chain)
for (i in 1:length(chains)) {
  for (j in i:length(chains)) {
    matrix_names <- c(matrix_names, paste0(chains[i], "_", chains[j]))
  }
}

# Initialize a named list of 378x378 zero matrices
matrix_list <- list()
for (name in matrix_names) {
  matrix_list[[name]] <- matrix(0, nrow = num_res, ncol = num_res)
}

# ==============================================================================
# FUNCTION: Process Tetramer Files
# ==============================================================================
process_tetramer_file <- function(file_path) {
  lines <- readLines(file_path)
  
  for (line in lines) {
    parts <- unlist(strsplit(trimws(line), "\\s+"))
    
    # Locate positions of all chain identifiers (items starting with /)
    slash_indices <- grep("^/", parts)
    if (length(slash_indices) < 2) next
    
    # Extract the interacting pair structural data
    chain1 <- gsub("/", "", parts[slash_indices[1]])
    res1 <- as.numeric(parts[slash_indices[1] + 2])
    
    chain2 <- gsub("/", "", parts[slash_indices[2]])
    res2 <- as.numeric(parts[slash_indices[2] + 2])
    
    # Apply mouse conversion rules: Discard 108, shift >= 109 down by 1
    if (is.na(res1) || is.na(res2)) next
    if (res1 == 108 || res2 == 108) next
    
    if (res1 >= 109) res1 <- res1 - 1
    if (res2 >= 109) res2 <- res2 - 1
    
    # Safety Check: Skip if residue extraction fails or falls out of bounds
    # Note: Handles non-standard entries safely (e.g., ligand lines like 1SY 401)
    if (res1 > num_res || res2 > num_res || res1 < 1 || res2 < 1) next
    
    # Verify that both extracted chains belong to our target set A, B, C, D
    if (!(chain1 %in% chains) || !(chain2 %in% chains)) next
    
    # Sort the chain characters alphabetically to match our list keys (e.g., B_A becomes A_B)
    sorted_chains <- sort(c(chain1, chain2))
    mat_key <- paste0(sorted_chains[1], "_", sorted_chains[2])
    
    # Map the interaction into the matrix ensuring complete symmetry
    if (chain1 == chain2) {
      # Intra-chain interaction
      matrix_list[[mat_key]][res1, res2] <<- matrix_list[[mat_key]][res1, res2] + 1
      matrix_list[[mat_key]][res2, res1] <<- matrix_list[[mat_key]][res2, res1] + 1
    } else {
      # Inter-chain interaction
      # Match the correct residue index to the sorted matrix label
      res_1st <- if (chain1 == sorted_chains[1]) res1 else res2
      res_2nd <- if (chain2 == sorted_chains[2]) res2 else res1
      
      matrix_list[[mat_key]][res_1st, res_2nd] <<- matrix_list[[mat_key]][res_1st, res_2nd] + 1
      matrix_list[[mat_key]][res_2nd, res_1st] <<- matrix_list[[mat_key]][res_2nd, res_1st] + 1
    }
  }
}

process_tetramer_file(file_contacts)
process_tetramer_file(file_hbonds)

# Save connectivity matrices
write.csv(data.frame(mat_A_A),'bonds_monomer.csv')
write.csv(data.frame(mat_A_B),'bonds_dimer.csv')
write.csv(data.frame(matrix_list$A_C+matrix_list$A_D),'bonds_tetramer.csv')

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

# --- 1. Load and Prepare Mets and STING Barplot Data ---
fitness_data <- read.csv("../cGAMP_sc/FACS_weights_all.csv")
mets_data <- read.csv('../fitness_position_subs_modeling/inferred_weights_Primary_vs_Met_4T1.csv')

df_bars <- data.frame(
  Residue = 1:378,
  Fitness = mets_data$additive
)

df_sting <- data.frame(
  Residue = 1:378,
  Fitness = fitness_data$STING_expression[1:378]
)

# --- Helper Function to Generate Arc Data ---
mets_cutoff <- sd(mets_data$additive)

generate_arc_dataframe <- function(contact_matrix, check_triangle = FALSE) {
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
    
    if (mets_data$additive[r1]>mets_cutoff | mets_data$additive[r2]>mets_cutoff){ 
    span <- abs(r1 - r2)
      midpoint <- (r1 + r2) / 2
      
      dip_depth <- -0.01 - (span * 0.003) 
      arc_id <- paste0("arc_", r1, "_", r2)
      
      df_list[[k]] <- data.frame(
        x = c(r1, midpoint, r2),
        y = c(0, dip_depth, 0) -0.15,# a -ve y offset to take it down
        Bonds = bonds,
        Group = arc_id)
    }
    
  }
  return(do.call(rbind, df_list))
}

# --- 2. Generate Arc Tracks ---
df_arcs_AA <- generate_arc_dataframe(mat_A_A, check_triangle = TRUE)
df_arcs_AB <- generate_arc_dataframe(mat_A_B, check_triangle = FALSE)
df_arcs_ACD <- generate_arc_dataframe(matrix_list$A_C+matrix_list$A_D,
                                      check_triangle = FALSE)

# --- 3. Build the Combined Plot Base ---
p <- ggplot() +
  # Draw the fitness bars
  geom_col(data = df_bars, aes(x = Residue, y = Fitness, fill = Fitness), width = 0.8) +
  
  # Continuous ColorBrewer Red-Blue (RdBu) palette
  scale_fill_distiller(
    palette = "RdBu",
    direction = -1,
    limits = c(-0.1, 0.1),
    oob = scales::squish,
    name = "Fitness",
    guide = "none"
  ) +
  
  # Format X-axis
  scale_x_continuous(
    limits = c(0, 379),
    breaks = seq(5, 378, by = 10),
    expand = c(0, 0)
  ) +
  
  # Minimal Presentation
  theme_cowplot(font_size=20) +
  theme(legend.position = "right",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
  labs(
    x = "Residue",
    y = "Metastatic fitness"
  )

# --- 4. Layer the Discrete Arc Classes Safely ---

# Layer 1: Monomer Arcs (A-A Interactions)
if (!is.null(df_arcs_AA) && nrow(df_arcs_AA) > 0) {
  p <- p + geom_bezier(
    data = df_arcs_AA,
    # Added 'color = "Monomer"' inside aes() to register it to the scale legend
    aes(x = x, y = y, group = Group, linewidth = Bonds, color = "Monomer"),
    alpha = 0.8
  )
}

# Layer 2: Dimer Arcs (A-B Interactions)
if (!is.null(df_arcs_AB) && nrow(df_arcs_AB) > 0) {
  p <- p + geom_bezier(
    data = df_arcs_AB,
    # Added 'color = "Dimer"' inside aes() to register it to the scale legend
    aes(x = x, y = y, group = Group, linewidth = Bonds, color = "Dimer"),
    alpha = 0.8
  )
}

# Layer 3: Tetramer Arcs (A-C and A-D Interactions)
if (!is.null(df_arcs_ACD) && nrow(df_arcs_ACD) > 0) {
  p <- p + geom_bezier(
    data = df_arcs_ACD,
    # Added 'color = "tetramer"' inside aes() to register it to the scale legend
    aes(x = x, y = y, group = Group, linewidth = Bonds, color = "Tetramer"),
    alpha = 0.8
  )
}

# --- 5. Final Scale Adjustments (Legends) ---

# Define the custom color scale and labels for the structural interactions
p <- p + scale_color_manual(
  values = c("Monomer" = "#5e3c99",
             "Dimer" = "#d01c8b",
             "Tetramer" = "#018571"),
  breaks = c("Monomer", "Dimer", "Tetramer"),
    name = "Bond type"
) 

# Apply global line weight scaling constraints
p <- p + scale_linewidth(
  range = c(0.075, 0.75),
  name = "# Bonds"
) 

# Ensure the arc alpha/transparency renders correctly inside the legend box
p <- p + guides(
  color = guide_legend(override.aes = list(alpha = 1, size = 1.0))
) + theme(legend.spacing.y = unit(1.5, "cm"),
          legend.justification = "bottom")

## Barplot showing residue-wise share of bond types
# Get monomer, dimer, and tetramer #bonds for hi mets
## Load contact matrices
bonds_AA <- read.csv('bonds_monomer.csv',row.names = 1)
bonds_AB <- read.csv('bonds_dimer.csv',row.names = 1)
bonds_ACD <- read.csv('bonds_tetramer.csv',row.names = 1)

mets_bonds <- data.frame()

for (pos in 1:378){
  mets_bonds <- rbind(mets_bonds,
                      data.frame("nbonds"=bonds_AA[,pos],
                                 "type" = "Monomer",
                                 "position" = pos),
                      
                      data.frame("nbonds"=bonds_AB[,pos],
                                 "type" = "Dimer",
                                 "position" = pos),
                      
                      data.frame("nbonds"=bonds_ACD[,pos],
                                 "type" = "Tetramer",
                                 "position" = pos)
  )
}

# Remove disconnected bonds (nbonds=0) and add fitness df
mets_bonds <- mets_bonds[mets_bonds$nbonds!=0,]

# 1. Aggregate total bonds per position and type
bonds_summary <- mets_bonds %>%
  group_by(position, type) %>%
  summarise(total_bonds = sum(nbonds), .groups = "drop") %>%
  # 3. Explicitly set factor levels for consistent legend/stack ordering
  mutate(type = factor(type, levels = c("Monomer", "Dimer", "Tetramer")))

# 4. Generate the proportional stacked bar plot (fractions)
p_bar <- ggplot(bonds_summary, aes(x = position, y = total_bonds, fill = type)) +
  geom_col(position = "fill",linewidth = 0.1, color = "white") + 
  scale_fill_manual(
    values = c("Monomer" = "#5e3c99",
               "Dimer" = "#d01c8b",
               "Tetramer" = "#018571")
  ) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 379), breaks = seq(5, 378, by = 10),
                     expand = c(0, 0)) +
  labs(
    x = "Residue Position",
    y = "% Bonds"
  ) +
  theme_cowplot() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=11),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    
    # Removes the gray background grid and panel
    panel.background = element_blank(), 
    # Removes the main outer background border/color
    plot.background = element_blank(),
    # Removes all axes completely
    axis.line = element_blank(),
  )
# Make bond # dotplots
bonds_summary <- mets_bonds %>% 
  group_by(position, type) %>% 
  summarise(total_bonds = sum(nbonds), .groups = "drop") %>% 
  mutate(type = factor(type, levels = c("Monomer", "Dimer", "Tetramer"))) %>% 
  # Fill in missing combinations from position 1 to 379 with 0
  complete(position = 1:378, type, fill = list(total_bonds = 0))

# 1. Create the Dotplot Function
# This creates a thin "strip" where dots appear if a bond exists at that residue
create_dotplot <- function(data, bond_type, color) {
  # Filter for the specific type and where bonds actually exist
  type_data <- data %>% filter(type == bond_type)
  
  ggplot(type_data, aes(x = position, y = total_bonds)) +
    geom_point(color = color, size = 0.8) +
    # Match the exact x-axis limits and expansion of your bar plot
    scale_x_continuous(limits = c(0, 379), breaks = seq(5, 378, by = 10),
                       expand = c(0, 0)) +
    ylab(bond_type) +
    theme_cowplot() +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.text.x = element_blank(),
        axis.title.x = element_blank(),
      axis.title.y = element_text(size=11),
        legend.position = "none",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
    )
}

# 2. Generate the individual dot strips
p_mono <- create_dotplot(bonds_summary, "Monomer", "#5e3c99")
p_dimer <- create_dotplot(bonds_summary, "Dimer", "#d01c8b")
p_tetra <- create_dotplot(bonds_summary, "Tetramer", "#018571")

# Make a dotplot for sequence conservation
seq_conservation <- read.csv('conservation.csv')
seq_conservation$position <- 1:378

p_conservation <-  ggplot(seq_conservation, aes(x = position, y = conservation)) +
  geom_point(color = 'black', size = 0.8) +
  geom_line() +
  # Match the exact x-axis limits and expansion of your bar plot
  scale_x_continuous(limits = c(0, 379), breaks = seq(5, 378, by = 10),
                     expand = c(0, 0)) +
  ylab("Seq\nconserved") +
  theme_cowplot() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=11),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )


# 3. Ensure your upper plot 'p' also shares matching scale expansion rules 
# Add this line to make sure plot 'p' doesn't use default padding that shifts it from p_bar
p <- p + scale_x_continuous(limits = c(0, 379),
                            breaks = seq(5, 378, by = 10),
                            expand = c(0, 0))

# Add STING barplot on top
p_sting <- ggplot() +
  # Draw the fitness bars
  geom_col(data = df_sting,
           aes(x = Residue, y = Fitness, fill = Fitness), width = 0.8) +
  
  # Continuous ColorBrewer Red-Blue (RdBu) palette
  scale_fill_distiller(
    palette = "RdBu",
    direction = -1,
    limits = c(-2, 2),
    breaks = c(-2, 0, 2),
    oob = scales::squish,
    name = NULL,
    guide = guide_colorbar(barheight = unit(1, "cm"),
                           barwidth = unit(0.3, "cm"))
  ) +
  
  # Format X-axis
  scale_x_continuous(
    limits = c(0, 379),
    expand = c(0, 0),
    breaks = NULL
  ) +
  
  # Minimal Presentation
  theme_cowplot(font_size=20) +
  theme(legend.position = "right",
        plot.margin = margin(t = 1, r = 0, b = 0, l = 0, unit = "pt"),
        # Removes X axes completely
        axis.line.x = element_blank()) +
  labs(
    x = NULL,
    y = "STING\nprotein"
  )


# Combine plots
p_full <- p_sting / p_tracks / p / p_bar / p_mono / p_dimer / p_tetra / p_conservation + plot_layout(heights = c(1, 1, 8, 1,1,1,1,1))

svglite("mets_barplot_bonds_LOFs.svg", width = 15, height = 10)
print(p_full)
dev.off()

# Save selected plots
p_save <- p +
  theme_cowplot(font_size=30) +
  scale_x_continuous(limits = c(0, 378),
                     breaks = c(0, seq(5, 378, by = 5)),
                     labels = NULL,
                     expand = c(0.0008, 0.0008)) +
  labs(x=NULL,y="Metastatic\nfitness")
  plot_layout(heights = c(6,1))
saveRDS(p_save, file = "mets_bars_arcs.rds")

