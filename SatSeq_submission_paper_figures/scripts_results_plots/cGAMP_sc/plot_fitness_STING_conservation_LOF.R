rm(list=ls()) # Clear environment
setwd(getSrcDirectory(function(){})) # set working directory to where this script is

library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpmisc)
library(ggvenn)
library(cowplot) # Added to ensure theme_cowplot() functions correctly

# ============================================================
# 1. Load and clean datasets
# ============================================================

# Metastatic fitness, 0-377 -> 1-378
mets_fitness <- read.csv("../fitness_position_subs_modeling/inferred_weights_Primary_vs_Met_4T1.csv")
names(mets_fitness)[1] <- "position" # Fixed dynamic renaming error
mets_fitness <- mets_fitness %>% 
  mutate(position = position + 1) %>% 
  select(position, additive, interaction)

# Sequence conservation, 0-377 -> 1-378
conservation <- read.csv("conservation.csv")
conservation <- conservation %>% 
  mutate(position = mouse_position + 1) %>% 
  select(position, mouse_AA, conservation)

# FACS weights: Aggregate by position to avoid row duplication during joins
facs_weights <- read.csv("FACS_weights_all.csv")
facs_weights_avg <- facs_weights %>% 
  mutate(position = suppressWarnings(as.numeric(gsub("^[A-Za-z]+", "", variant)))) %>% 
  filter(!is.na(position)) %>% 
  group_by(position) %>% 
  summarize(STING_expression = mean(STING_expression, na.rm = TRUE), .groups = "drop")

# LOF variants: Aggregate pathway flags by residue position
lof_df <- read.csv("LOF_variants.csv")
lof_flags_avg <- lof_df %>% 
  mutate(
    IFN  = grepl("IFNa", group),
    NFkB = grepl("NFkB", group),
    UPR  = grepl("UPR", group)
  ) %>% 
  group_by(residue) %>% 
  summarize(
    IFN  = any(IFN, na.rm = TRUE),
    NFkB = any(NFkB, na.rm = TRUE),
    UPR  = any(UPR, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 2. Master Merge (One row per residue position 1-378)
# ============================================================
per_residue <- data.frame(position = 1:378) %>% 
  left_join(conservation, by = "position") %>% 
  left_join(mets_fitness, by = "position") %>% 
  left_join(facs_weights_avg, by = "position") %>% 
  left_join(lof_flags_avg, by = c("position" = "residue")) %>% 
  mutate(
    IFN     = coalesce(IFN, FALSE),
    NFkB    = coalesce(NFkB, FALSE),
    UPR     = coalesce(UPR, FALSE),
    any_LOF = IFN | NFkB | UPR,
    # Construct clean variant string using WT amino acid + position (e.g., M1, S365)
    variant_label = paste0(mouse_AA, position) 
  )

# Separate dataset for ggrepel labeling based on user condition (additive > 0.1)
labeled_residues <- per_residue %>% 
  filter(additive > 0.1)

# ============================================================
# 3. Plot 1: metastatic additive fitness vs. STING protein level
# ============================================================
p1 <- ggplot(per_residue, aes(x = STING_expression, y = additive)) +
  geom_point(color = "grey30", alpha = 0.7, size = 1.8, na.rm = TRUE) +
  geom_text_repel(
    data = labeled_residues,
    aes(label = variant_label),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    na.rm = TRUE
  ) +
  stat_correlation(
    method = "pearson", 
    label.y = "top", 
    label.x = "left", 
    output.type = "text", 
    mapping = aes(label = paste(after_stat(r.label), after_stat(p.value.label), sep = "\n")), 
    size = 4
  ) +
  labs(x = "STING protein level", y = "Metastatic fitness") +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave("scatter_STING-vs-metastatic-fitness.svg", p1, width = 4, height = 4,
       device = svglite::svglite)

# ============================================================
# 4. Plot 2: metastatic additive fitness vs. sequence conservation
# ============================================================
p2 <- ggplot(per_residue, aes(x = conservation, y = additive)) +
  geom_point(color = "grey30", alpha = 0.7, size = 1.8, na.rm = TRUE) +
  geom_text_repel(
    data = labeled_residues,
    aes(label = variant_label),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    na.rm = TRUE
  ) +
  stat_correlation(
    method = "pearson", 
    label.y = "top", 
    label.x = "left", 
    output.type = "text", 
    mapping = aes(label = paste(after_stat(r.label), after_stat(p.value.label), sep = "\n")), 
    size = 4
  ) +
  labs(x = "Sequence conservation", y = "Metastatic fitness") +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave("scatter_conservation-vs-mets.svg", p2, width = 4, height = 4,
       device = svglite::svglite)

# ============================================================
# 5. Plot 3: metastatic additive fitness by LOF group
#
#    4 groups, overlapping: a residue counts toward IFN/NFkB/UPR LOF
#    independently (whichever pathways it's annotated as hitting), and
#    separately toward "Non-LOF" only if it's in none of them.
# ============================================================
box_df <- bind_rows(
  per_residue %>% filter(IFN)      %>% mutate(LOF_group = "IFN LOF"),
  per_residue %>% filter(NFkB)     %>% mutate(LOF_group = "NFkB LOF"),
  per_residue %>% filter(UPR)      %>% mutate(LOF_group = "UPR LOF"),
  per_residue %>% filter(!any_LOF) %>% mutate(LOF_group = "Non-LOF")
)
box_df$LOF_group <- factor(box_df$LOF_group,
                           levels = c("IFN LOF", "NFkB LOF", "UPR LOF", "Non-LOF"))

lof_colors <- c(
  "IFN LOF"  = "#e3312a",
  "NFkB LOF" = "#35b6e3",
  "UPR LOF"  = "#7ac943",
  "Non-LOF"  = "gray20"   # "deep gray" — this exact shade is a judgment call
)

p3 <- ggplot(box_df, aes(x = LOF_group, y = additive, color = LOF_group, fill = LOF_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25, na.rm = TRUE) +
  geom_jitter(width = 0.15, height = 0, size = 1.5, alpha = 0.8, na.rm = TRUE) +
  scale_color_manual(values = lof_colors) +
  scale_fill_manual(values = lof_colors) +
  labs(x = NULL, y = "Metastatic fitness") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
            legend.position = "none")  # width:height = 2:1

ggsave("boxplot_mets-by-LOF-group.svg",
       p3, width = 4, height = 4, device = svglite::svglite)

######### Venn diagram ############
# LOF pathway groups — a residue can belong to more than one (overlapping,
# same as the boxplot earlier)
ifn_residues  <- lof_df$residue[grepl("IFNa", lof_df$group)]
nfkb_residues <- lof_df$residue[grepl("NFkB", lof_df$group)]
upr_residues  <- lof_df$residue[grepl("UPR",  lof_df$group)]

# "Hi mets": additive metastatic fitness more than 1 SD above the mean
mets_sd   <- sd(mets_fitness$additive, na.rm = TRUE)
hi_mets_residues <- mets_fitness$position[mets_fitness$additive > (mets_sd)]

venn_sets <- list(
  "IFN LOF"  = ifn_residues,
  "NFkB LOF" = nfkb_residues,
  "UPR LOF"  = upr_residues,
  "Metastatic"  = hi_mets_residues
)

# ============================================================
# Venn diagram
#
# LOF circles keep the same fill colors used throughout (IFN/NFkB/UPR).
# "Hi mets" is unfilled (fill = NA) so only its black outline shows —
# stroke_color is one value applied to every circle and is already
# black by default, so no extra styling is needed to get that outline.
# ============================================================
p_venn <- ggvenn(
  venn_sets,
  fill_color = c("#e3312a", "#35b6e3", "#7ac943", NA),
  fill_alpha = 0.55,
  stroke_color = c("#e3312a", "#35b6e3", "#7ac943", "black"),
  stroke_size = 1,
  set_name_size = 4,
  text_size = 3.5,
  show_percentage = FALSE
  
)

ggsave("venn_mets-and-LOFs.svg",
       p_venn, width = 3, height = 3, device = svglite::svglite)