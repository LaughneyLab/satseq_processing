rm(list=ls())#Clear environment

library(tidyverse)
library(ggrepel)
library(stringr)
library(ggpmisc)
library(cowplot)
library(scales)
library(svglite)
library(ggpubr)


setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# Load FACS data
facs_weights <- read.csv('FACS_weights_all.csv')

# List of STING-dependent genes which are down in STOP vs REF in 4T1 cGAMP
STING_df <- rbind(read.csv('4T1_interaction/fcHurdle_STOP.csv'),
                  read.csv('EO771_interaction/fcHurdle_STOP.csv'))

STING_genes <- unique(STING_df$primerid[STING_df$fdr < 0.5])

df_cgamp_stop <- read.csv('4T1_variant_vs_REF_cGAMP/fcHurdle_4T1_STOP.csv')
STING_genes <- intersect(STING_genes,df_cgamp_stop$primerid[df_cgamp_stop$coef<0])

# Read Hallmark genes and intersect with STING genes
hallmark_df <- read.csv('mouse_hallmarks_cGAMP_extended.csv')

# ER Stress
er_genes <- intersect(STING_genes,
                      unique(hallmark_df$Gene[hallmark_df$Pathway %in%
                                      c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                                        "HALLMARK_UNFOLDED_PROTEIN_RESPONSE_CORE"
                                        )]))
# IFNa response
ifna_genes <- intersect(STING_genes,
                        unique(hallmark_df$Gene[hallmark_df$Pathway %in%
                                        c("HALLMARK_INTERFERON_ALPHA_RESPONSE")]))
# NFkB
nfkb_genes <- intersect(STING_genes,
                        unique(hallmark_df$Gene[hallmark_df$Pathway %in%
                                        c("HALLMARK_TNFA_SIGNALING_VIA_NFKB")]))
# Combine to get all genes
all_hallmark_genes <- unique(c(er_genes, ifna_genes, nfkb_genes))

# Create a dataframe mapping genes to pathways
gene_pathway_map <- tibble(Gene = all_hallmark_genes) %>%
  mutate(
    `ER Stress` = Gene %in% er_genes,
    `IFN Response` = Gene %in% ifna_genes,
    `NFkB signaling` = Gene %in% nfkb_genes
  ) %>%
  # Create a combined label for genes in multiple pathways, or just one
  mutate(Pathway = case_when(
    `ER Stress` & `IFN Response` & `NFkB signaling` ~ "All 3",
    `ER Stress` & `IFN Response` ~ "ER stress + IFN",
    `IFN Response` & `NFkB signaling` ~ "IFN + NFkB",
    `ER Stress` & `NFkB signaling` ~ "ER stress + NFkB",
    `ER Stress` ~ "ER stress",
    `IFN Response` ~ "IFN response",
    `NFkB signaling` ~ "NFkB signaling",
    TRUE ~ "Other"
  ))

# -----------------------------
# Load gene by variant data
# -----------------------------
load_scores <- function(directory, cell_line) {
  files <- list.files(directory,
                      pattern = paste0("^fcHurdle_", cell_line, "_.*\\.csv$"),
                      full.names = TRUE)
  
  bind_rows(lapply(files, function(f) {
    variant_name <- gsub(paste0("fcHurdle_", cell_line, "_|.csv"), "",
                         basename(f))
    
    df <- read.csv(f)
    df$score <- -log10(df$fdr + 1e-300) * sign(df$coef)
    df$variant <- variant_name
    
    df[, c("primerid", "score", "variant")]
  }))
}

sc_4t1 <- load_scores("4T1_variant_vs_REF_cGAMP", "4T1") %>%
  filter(primerid %in% all_hallmark_genes) %>%
  pivot_wider(names_from = variant, values_from = score) %>%
  column_to_rownames("primerid") %>%
  as.matrix()

sc_EO771 <- load_scores("EO771_variant_vs_REF_cGAMP", "EO771") %>%
  filter(primerid %in% all_hallmark_genes) %>%
  pivot_wider(names_from = variant, values_from = score) %>%
  column_to_rownames("primerid") %>%
  as.matrix()

# 3. Function to average scores for a specific gene set
get_pathway_avg <- function(genes, label, heat_mat) {
  valid_genes <- intersect(genes, rownames(heat_mat))
  colMeans(heat_mat[valid_genes, , drop = FALSE], na.rm = TRUE) %>%
    enframe(name = "variant", value = label)
}

# 4. Combine into final dataframe
sc_4T1 <- get_pathway_avg(er_genes, "er_stress",sc_4t1) %>%
  left_join(get_pathway_avg(ifna_genes, "ifna_response",sc_4t1), by = "variant") %>%
  left_join(get_pathway_avg(nfkb_genes, "nfkb_response",sc_4t1), by = "variant")

sc_EO771 <- get_pathway_avg(er_genes, "er_stress",sc_EO771) %>%
  left_join(get_pathway_avg(ifna_genes, "ifna_response",sc_EO771), by = "variant") %>%
  left_join(get_pathway_avg(nfkb_genes, "nfkb_response",sc_EO771), by = "variant")

write.csv(sc_4T1, "sc_MAST_pathway_score_4T1.csv")

# 2. Reshape and Merge Single Cell Scores
# We pivot both to long format so we can join them by variant and pathway
sc_4T1_long <- sc_4T1 %>%
  pivot_longer(cols = c(er_stress, ifna_response, nfkb_response), 
               names_to = "pathway", values_to = "score_4T1")

sc_EO771_long <- sc_EO771 %>%
  pivot_longer(cols = c(er_stress, ifna_response, nfkb_response), 
               names_to = "pathway", values_to = "score_EO771")

merged_sc <- inner_join(sc_4T1_long, sc_EO771_long, by = c("variant", "pathway"))

# 3. Join FACS Weights and Map Pathways to specific weight columns
# We create a mapping to assign the correct weight to each pathway row
plot_data <- merged_sc %>%
  left_join(facs_weights, by = "variant") %>%
  mutate(
    weight = case_when(
      pathway == "nfkb_response" ~ NFkB_expression,
      pathway == "er_stress"     ~ ATF3_expression,
      pathway == "ifna_response" ~ Dual_IFN_expression
    ),
    abs_weight = abs(weight)
  )

# Rename for printing
plot_data$pathway = case_when(
  plot_data$pathway == "nfkb_response" ~ "NFkB signaling",
  plot_data$pathway == "er_stress"     ~ "Unfolded protein response",
  plot_data$pathway == "ifna_response" ~ "IFNa response"
)

# 1. Pre-calculate the stats so they don't change with the scale
stats_data <- plot_data %>%
  group_by(pathway) %>%
  do(model = lm(score_EO771 ~ score_4T1, data = .)) %>%
  summarize(
    pathway = pathway,
    r_val = sqrt(summary(model)$r.squared),
    p_val = summary(model)$coefficients[2,4],
    # Create the label string manually
    label = sprintf('R~"="~%.3f~","~italic(P)~"="~%.2g', r_val, p_val)
  )

# 4. Generate Plot

# Get LOF labels
plot_data <- plot_data %>%
  group_by(pathway) %>%
  mutate(
    # Calculate pathway-specific mean and SD
    mean_4T1 = mean(score_4T1, na.rm = TRUE),
    sd_4T1 = sd(score_4T1, na.rm = TRUE),
    mean_EO771 = mean(score_EO771, na.rm = TRUE),
    sd_EO771 = sd(score_EO771, na.rm = TRUE),
    
    # --- YOUR LOGIC ---
    # Criteria: Score is < 1SD below the mean in EITHER cell line
    # AND the direction of effect (sign) is the same in both
    is_hit_score = (score_4T1 < (mean_4T1 - 1 * sd_4T1) | score_EO771 < (mean_EO771 - 1 * sd_EO771)) & 
      (sign(score_4T1) == sign(score_EO771)),
    
    # Combined logic: (Hit score AND weight < -1) OR (weight < -2)
    # Note: Assuming 'weight' is the correct column based on your colnames, 
    # as 'facs_weight' was not in the provided list.
    final_hit = (is_hit_score & (weight < -1)) | (weight < -2)
  ) %>%
  ungroup() # Good practice to ungroup after mutate

p1 <- ggplot(plot_data, aes(x = score_4T1, y = score_EO771)) +
  geom_smooth(method = "lm", linewidth = 0.5, color='#7393B3') +
  geom_point(aes(color = weight, size = abs_weight), alpha = 0.7) + 
  geom_text(data = stats_data, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1.5, parse = TRUE, size = 4) +
  scale_x_continuous(trans = pseudo_log_trans(sigma = 0.05, base = 10),
                     breaks = scales::breaks_pretty(n = 4)) + 
  scale_y_continuous(trans = pseudo_log_trans(sigma = 0.05, base = 10),
                     breaks = scales::breaks_pretty(n = 4)) +
  scale_size_continuous(range = c(0.2, 3)) +
  facet_wrap(~pathway, scales = "free") +
  scale_color_distiller(
    palette = "RdBu",
    limits = c(-1.5, 1.5),
    oob = scales::squish) +
  theme_cowplot() +
  geom_text_repel(
    data = plot_data %>% filter(final_hit == TRUE),
    aes(label = variant),
    size = 3,
    max.overlaps = 30
  ) +
  geom_vline(xintercept = 0, color = "#7393B3", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "#7393B3", linetype = "dashed") +
  labs(
    x = "4T1 single cell score",
    y = "EO771 single cell score",
    color = "FACS score",
    size = "|FACS score|"
  ) +
  # Force the colorbar to order 1 (top) and size legend to order 2 (bottom)
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  )

# Save
# ggsave(plot=p1, 'sc_scores_cell_lines_FACS_colored.pdf',
#        width=15, height=4)

svglite("sc_scores_cell_lines_FACS_colored.svg", width =16, height = 4)
print(p1)
dev.off()

saveRDS(p1, file = "sc_scores_cell_lines_FACS_colored.rds")



# Define sign consistency across cell lines
# 2. Join and define color logic based on sign agreement
plot_data_facs <- sc_4T1_long %>%
  inner_join(sc_EO771_long, by = c("variant", "pathway")) %>%
  mutate(
    point_color = case_when(
      score_4T1 > 0 & score_EO771 > 0 ~ "Both up",
      score_4T1 < 0 & score_EO771 < 0 ~ "Both down",
      TRUE                           ~ "Different"
    )
  ) %>%
  # Join FACS weights
  left_join(facs_weights, by = "variant") %>%
  # Map specific FACS columns to the pathway in the row
  mutate(
    facs_weight = case_when(
      pathway == "nfkb_response" ~ NFkB_expression,
      pathway == "er_stress"     ~ ATF3_expression,
      pathway == "ifna_response" ~ Dual_IFN_expression
    )
  )

# 3. Generate the Plot
# Rename for printing
plot_data_facs$pathway = case_when(
  plot_data_facs$pathway == "nfkb_response" ~ "NFkB signaling",
  plot_data_facs$pathway == "er_stress"     ~ "Unfolded protein response",
  plot_data_facs$pathway == "ifna_response" ~ "IFNa response"
)

# Get LOF labels
plot_data_facs <- plot_data_facs %>%
  group_by(pathway) %>%
  mutate(
    # Calculate pathway-specific mean and SD
    mean_4T1 = mean(score_4T1, na.rm = TRUE),
    sd_4T1 = sd(score_4T1, na.rm = TRUE),
    mean_EO771 = mean(score_EO771, na.rm = TRUE),
    sd_EO771 = sd(score_EO771, na.rm = TRUE),
    
    # --- YOUR LOGIC ---
    # Criteria: Score is < 1SD below the mean in EITHER cell line
    # AND the direction of effect (sign) is the same in both
    is_hit_score = (score_4T1 < (mean_4T1 - 1 * sd_4T1) | score_EO771 < (mean_EO771 - 1 * sd_EO771)) & 
      (sign(score_4T1) == sign(score_EO771)),
    
    # Combined logic: (Hit score AND weight < -1) OR (weight < -2)
    # Note: Assuming 'weight' is the correct column based on your colnames, 
    # as 'facs_weight' was not in the provided list.
    final_hit = (is_hit_score & (facs_weight < -1)) | (facs_weight < -2)
  ) %>%
  ungroup() # Good practice to ungroup after mutate


# 1. Pre-calculate the stats so they don't change with the scale
stats_data <- plot_data_facs %>%
  group_by(pathway) %>%
  do(model = lm(facs_weight ~ score_4T1, data = .)) %>%
  summarize(
    pathway = pathway,
    r_val = sqrt(summary(model)$r.squared),
    p_val = summary(model)$coefficients[2,4],
    # Create the label string manually
    label = sprintf('R~"="~%.3f~","~italic(P)~"="~%.2g', r_val, p_val)
  )

# 2. Build the plot
p2 <- ggplot(plot_data_facs, aes(x = score_4T1, y = facs_weight)) + 
  geom_smooth(method = "lm", linewidth = 0.5, color='#7393B3') + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7393B3") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "#7393B3") + 
  geom_point(aes(color = point_color), size = 1, alpha = 0.6) + 
  # Use geom_text instead of stat_fit_glance to use our pre-calculated 'real' stats
  geom_text(data = stats_data, 
            aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1.5, parse = TRUE, size = 4) +
  # Apply pseudo_log here (this is what brings the point back for geom_text_repel)
  scale_x_continuous(trans = pseudo_log_trans(sigma = 0.05, base = 10),
                     breaks = scales::breaks_pretty(n = 4)) + 
  facet_wrap(~pathway, scales = "free") + 
  scale_color_manual(values = c(
    "Both up" = "#e41a1c", 
    "Both down" = "#377eb8", 
    "Different" = "#d9d9d9"
  )) + 
  theme_cowplot() + 
  # Add labels for high/low outliers (adjust quantiles as needed)
  geom_text_repel(
    # Filter using the newly created column
    data = plot_data_facs %>% filter(final_hit == TRUE),
    aes(label = variant),
    size = 3,
    max.overlaps = 30
  ) + 
  labs(
    x = "4T1 Single Cell Score", 
    y = "FACS score", 
    color = "4T1 & EO771"
  )

svglite("sc_vs_FACS.svg", width = 16, height = 4)
print(p2)
dev.off()

saveRDS(p2, file = "sc_vs_FACS.rds")


# We define the "Top" (bottom-ranked) variants based on your thresholds
sc_combined <- inner_join(
  sc_4T1 %>% pivot_longer(-variant, names_to = "pathway", values_to = "score_4T1"),
  sc_EO771 %>% pivot_longer(-variant, names_to = "pathway", values_to = "score_EO771"),
  by = c("variant", "pathway")
) %>%
  group_by(pathway) %>%
  mutate(
    # Calculate pathway-specific mean and SD for both cell lines
    mean_4T1 = mean(score_4T1, na.rm = TRUE),
    sd_4T1   = sd(score_4T1, na.rm = TRUE),
    mean_EO771 = mean(score_EO771, na.rm = TRUE),
    sd_EO771   = sd(score_EO771, na.rm = TRUE),
    
    # Criteria: Score is < 2SD below the mean in EITHER cell line 
    # AND the direction of effect (sign) is the same in both
    is_hit = (score_4T1 < (mean_4T1 - 1 * sd_4T1) | 
                score_EO771 < (mean_EO771 - 1 * sd_EO771)) & 
      (sign(score_4T1) == sign(score_EO771))
  ) %>%
  ungroup()


facs_long <- facs_weights %>%
  pivot_longer(
    cols = c(NFkB_expression, ATF3_expression, Dual_IFN_expression),
    names_to = "facs_type", values_to = "facs_weight"
  ) %>%
  mutate(pathway = case_when(
    facs_type == "NFkB_expression" ~ "nfkb_response",
    facs_type == "ATF3_expression" ~ "er_stress",
    facs_type == "Dual_IFN_expression" ~ "ifna_response"
  )) %>%
  group_by(pathway) %>%
  mutate(facs_rank = rank(facs_weight)) %>%
  ungroup()


# 3. Pivot FACS and join with SC hits
plot_data <- facs_long %>%
  left_join(sc_combined, by = c("variant", "pathway")) %>%
  mutate(
    # 1. Standardize pathway names first for easier logic mapping
    pathway = case_when(
      pathway == "nfkb_response" ~ "NFkB signaling",
      pathway == "er_stress"     ~ "Unfolded protein response",
      pathway == "ifna_response" ~ "IFNa response",
      TRUE ~ pathway
    ),
    
    # 2. Define the hit logic
    is_hit = ((is_hit == TRUE) & (facs_weight < -1)) | (facs_weight < -2),
    
    # 3. Assign specific LOF labels for hits, otherwise "Other (Gray)"
    highlight_color = case_when(
      !is_hit ~ "Other (Gray)",
      is_hit & pathway == "NFkB signaling" ~ "NFkB LOF",
      is_hit & pathway == "Unfolded protein response" ~ "UPR LOF",
      is_hit & pathway == "IFNa response" ~ "IFN LOF",
      TRUE ~ "Both -" # Fallback if a pathway doesn't match the three above
    )
  )

p3 <- ggplot(plot_data, aes(x = facs_rank, y = facs_weight)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7393B3") +
  geom_point(data = filter(plot_data, highlight_color == "Other (Gray)"), 
             color = "gray40", size = 1.5, alpha = 0.4,show.legend = FALSE) +
  geom_point(data = filter(plot_data, highlight_color != "Other (Gray)"), 
             aes(color = highlight_color), size = 1,show.legend = FALSE) +
  facet_wrap(~pathway, scales = "free") +
  scale_color_manual(values = c("IFN LOF" = "#e3312a",
                                "UPR LOF" = "#35b6e3",
                                "NFkB LOF" = "#7ac943")) +
  theme_cowplot() +
  geom_text_repel(
    data = filter(plot_data, highlight_color != "Other (Gray)"),
    aes(label = variant), color = 'black',
    size = 2.5, max.overlaps = 40, box.padding = 0.1,
    min.segment.length = 0.8,show.legend = FALSE
  ) +
  labs(
   x = "FACS Rank", y = "FACS Weight", color = "4T1 & EO771"
  )

svglite("FACS_rank_plot_LOF_colored.svg", width = 14, height = 4)
print(p3)
dev.off()

# 1. Identify the "hits" for each pathway from the existing plot_data
# This assumes plot_data exists from the previous code block
ifna_variants <- plot_data %>%
  filter(highlight_color == "IFN LOF") %>%
  pull(variant)

nfkb_variants <- plot_data %>%
  filter(highlight_color == "NFkB LOF") %>%
  pull(variant)

er_stress_variants <- plot_data %>%
  filter(highlight_color == "UPR LOF") %>%
  pull(variant)

# 1. Combine all unique variants from the three lists
all_relevant_variants <- unique(c(ifna_variants, nfkb_variants, er_stress_variants))

# 2. Create the data frame and map groups based on intersections
export_data <- data.frame(variant = all_relevant_variants) %>%
  mutate(
    group = case_when(
      # Triple Intersection
      variant %in% ifna_variants & variant %in% nfkb_variants & variant %in% er_stress_variants ~ "IFNa & NFkB & UPR",
      
      # Double Intersections
      variant %in% ifna_variants & variant %in% nfkb_variants ~ "IFNa & NFkB",
      variant %in% ifna_variants & variant %in% er_stress_variants ~ "IFNa & UPR",
      variant %in% nfkb_variants & variant %in% er_stress_variants ~ "NFkB & UPR",
      
      # Single Pathways
      variant %in% ifna_variants ~ "Only IFNa",
      variant %in% nfkb_variants ~ "Only NFkB",
      variant %in% er_stress_variants ~ "Only UPR"
    ),
    
    # 4. Extract numeric residue for sorting
    residue = as.numeric(gsub("[^0-9]", "", variant))
  ) %>%
  arrange(residue)

# 3. Save to CSV
export_data <- export_data[!export_data$variant %in% c('STOP','dCTT','RRDD'),]
write.csv(export_data, "LOF_variants.csv", row.names = FALSE)

## Sequence conservation vs pathway FACS LOF
# Load seq conservation df and merge with FACS positional weights
seq_conservation <- read.csv('conservation.csv')
facs_weights <- facs_weights[1:378,]
facs_weights$conservation <- seq_conservation$conservation

# 1. Transform to long form
facs_long <- facs_weights %>%
  pivot_longer(
    cols = c(Dual_IFN_expression, ATF3_expression, NFkB_expression),
    names_to = "pathway",
    values_to = "pathway_score"
  ) %>%
  mutate(pathway = recode(pathway,
                          "Dual_IFN_expression" = "IFN",
                          "ATF3_expression" = "ATF3",
                          "NFkB_expression" = "NFkB"
  ))

pathway_colors <- c(
  "IFN"  = "#e3312a", 
  "ATF3"  = "#35b6e3", 
  "NFkB" = "#7ac943"
)

p_cons <- ggplot(facs_long,
       aes(x = conservation, y = pathway_score, color = pathway, fill=pathway)) +
  geom_smooth(method = "lm", linewidth = 0.5) +
  geom_point(alpha=0.5, size=0.8) +
  stat_fit_glance(method = "lm",
                  label.y = "bottom",
                  label.x = 'right',
                  mapping = aes(label = sprintf('R~"="~%.3f~","~italic(P)~"="~%.2g',
                                                after_stat(r.squared^0.5), after_stat(p.value)),
                                size=6),
                  parse = TRUE) +
  theme_cowplot() +
  scale_color_manual(values = pathway_colors) +
  scale_fill_manual(values = pathway_colors) +
  geom_hline(yintercept=0, color = "#7393B3", linetype = "dashed") +
  guides(fill = "none") +
  labs(
    x = "Sequence conservation",
    y = "FACS score",
    color = NULL,
  )

svglite("FACS_vs_sequence_conservation.svg", width =6, height = 4)
print(p_cons)
dev.off()

saveRDS(p_cons, file = "FACS_vs_sequence_conservation.rds")


## Do boxplot of conservation LOF vs rest
facs_long <- facs_long[facs_long$pathway=='IFN',] # For unique values
facs_long$group <- 'Other\nvariants'
facs_long$group[facs_long$variant %in% export_data$variant] <- 'LOF\nvariants'

p <- ggplot(facs_long, aes(x = group, y = conservation)) + 
       geom_boxplot(outlier.shape = NA) + 
       geom_jitter(alpha = 0.45, width = 0.25, size=0.8) + 
       theme_cowplot() + 
       labs(x = NULL, 
            y = "STING sequence conservation") +
       stat_compare_means(comparisons = list( c("LOF\nvariants", "Other\nvariants") ),
                          method = "wilcox.test",
                          method.args = list(alternative = "greater"),
                          label = "p.format")

svglite("sequence_conservation_LOF.svg", width = 3, height = 4)
print(p)
dev.off()



