rm(list=ls())#Clear environment

library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(grid)
library(patchwork)
library(ggplot2)
library(svglite)
library(cowplot)
library(ggrepel)
library(ggpmisc)
library(scales)



setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# empty df
ranked_variants = data.frame()

for (cell_line in c("4T1","EO771")){

# --- 2. Load and Process Score Matrix ---
input_dir <- paste0(cell_line,"_variant_vs_REF_cGAMP")
files <- list.files(input_dir, pattern = paste0("^fcHurdle_",cell_line,
                                                "_.*\\.csv$"),
                    full.names = TRUE)


all_data_list <- lapply(files, function(filepath) {
  variant_name <- gsub(paste0("fcHurdle_",cell_line,"_|.csv"),
                       "", basename(filepath))
  df <- read.csv(filepath)
  df %>% 
    mutate(
      score = -log10(fdr + 1e-300) * sign(coef), 
      variant = variant_name,
    ) %>% 
    select(primerid, score, variant)
})

all_data_combined <- bind_rows(all_data_list) %>%
  select(primerid, score, variant) %>%
  pivot_wider(names_from = variant, values_from = score) %>%
  as.data.frame()

rownames(all_data_combined) <- all_data_combined$primerid
all_data_combined$primerid <- NULL

# --- 3. Cumsum Filtering (Variant Selection) ---
# Calculate absolute score sums per variant across ALL genes
variant_impact <- colSums(abs(all_data_combined), na.rm = TRUE)

ranked_variants <- rbind(ranked_variants, 
                         data.frame(
  variant = names(variant_impact),
  cell_line = cell_line,
  total_score = as.numeric(variant_impact)
) %>% 
  arrange(desc(total_score)) %>%
  mutate(cumsum_score = cumsum(total_score), variant_rank = row_number())
)

}

# Add FACS data
FACS_weights <- read.csv('FACS_weights_all.csv')

ranked_variants <- ranked_variants %>%
  left_join(FACS_weights, by = "variant") %>%
  mutate(
    label = ifelse(
      # Condition 1: Rank is within the top 130
      variant_rank < 130 & 
        # Condition 2: Expression is within the bottom 50
        min_rank(Dual_IFN_expression) <= 75,
      
      as.character(variant),
      NA_character_
    )
  )


# Plot elbows colored by FACS IFN showing specific variants
p_elbow <- ggplot(ranked_variants, aes(x = variant_rank, y = cumsum_score)) +
  geom_point(aes(color = Dual_IFN_expression),size=0.5,alpha=0.65,
             show.legend = TRUE) +
  geom_point(data=ranked_variants[abs(ranked_variants$Dual_IFN_expression)>1.5,], aes(color = Dual_IFN_expression),size=0.75,alpha=0.65,
             show.legend = TRUE) +
  scale_color_distiller(name = "IFN FACS score",
                        palette = "RdBu",
                        limits = c(-2, 2), oob = scales::squish) + 
  theme_cowplot() +
  facet_wrap(~cell_line, nrow = 1,scales='free') +
  geom_label_repel(aes(label = label),
                   size = 3, max.overlaps = 20, alpha = 0.8, 
                   fill = NA,
                   label.size = NA, # Use NA to remove the border box
                   min.segment.length = unit(0.1, 'lines'), show.legend = FALSE) +
  labs(x = "Variant rank", y = "Cumulative MAST scores") +
  # Move legend inside
  theme(
    legend.position.inside = c(0.85, 0.5), # Bottom right (x, y)
    legend.background = element_rect(fill = alpha("white", 0.5)), # Semi-transparent background
    legend.box.background = element_blank(),
    legend.title = element_text(size = 12), # Title size
    legend.text = element_text(size = 12),
  )

# Save plot
# ggsave(plot=p7,'DEG_cumsum_IFN_colored.pdf',
#        width=6, height=6)

svglite("MAST_cumsum_IFN_colored_4T1_EO771.svg", width = 8, height = 4)
print(p_elbow)
dev.off()

saveRDS(p_elbow, file = "MAST_cumsum_IFN_colored_4T1_EO771.rds")

## Do correlation plot of ranks

# Wide format
wide_variants <- ranked_variants %>%
  # 1. Select only needed columns for pivoting to keep it clean
  select(variant, cell_line, total_score, 
         Dual_IFN_expression, label) %>%
  
  # 2. Pivot variant_rank based on cell_line
  pivot_wider(
    names_from = cell_line, 
    values_from = total_score,
    names_prefix = "score_"
  ) %>%
  mutate(abs_IFN = abs(Dual_IFN_expression))

# 1. Pre-calculate the stats so they don't change with the scale
stats_data <- wide_variants %>% 
  do(model = lm(score_EO771 ~ score_4T1, data = .)) %>% 
  summarize( 
    r_val = sqrt(summary(model[[1]])$r.squared), 
    p_val = summary(model[[1]])$coefficients[2,4], 
    # Use atop() to split the expression cleanly onto 2 lines
    label = sprintf('atop(R == "%.3f", italic(P) == "%.2g")', r_val, p_val)
  )

p_corr <- ggplot(wide_variants, aes(x = score_4T1, y = score_EO771)) +
  geom_smooth(method = "lm", linewidth = 0.5, color='#7393B3') +
  geom_point(aes(color = Dual_IFN_expression,
                 size = squish(abs_IFN, range = c(0.5,2))), 
             alpha = 0.8) +
  geom_point(data=wide_variants[abs(wide_variants$Dual_IFN_expression)>1.5,],
             aes(color = Dual_IFN_expression,
                 size = squish(abs_IFN, range = c(0.5,2))), 
             alpha = 0.8,
             show.legend = TRUE) +
  scale_color_distiller(
    palette = "RdBu", 
    name = "IFN FACS score",
    limits = c(-2, 2),
    oob = squish
  ) +
  geom_text(data = stats_data, 
            aes(x = 0.1, y = 5000, label = label), 
            hjust = -0.1, vjust = 1.5, parse = TRUE, size = 4) +
  scale_x_continuous(trans = "log10") + 
  scale_y_continuous(trans = "log10") +
  scale_size_continuous(range = c(0.5, 2), name = "Score magnitude") +
  theme_cowplot() +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 15) +
  labs(x = "4T1 gene score", y = "EO771 gene score") +
  # Force colorbar legend (order = 1) above size legend (order = 2)
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  )
# Save svg for main fig
svglite("MAST_rank_correlation_IFN_colored.svg", width = 7, height = 5)
print(p_corr)
dev.off()

saveRDS(p_corr, file = "MAST_rank_correlation_IFN_colored.rds")
