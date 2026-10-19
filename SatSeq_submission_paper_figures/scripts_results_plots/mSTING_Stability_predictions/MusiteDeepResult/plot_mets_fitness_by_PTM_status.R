rm(list=ls()) # Clear environment
setwd(getSrcDirectory(function(){})) # set working directory to where this script is

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

ptm_file     <- "Prediction_results.txt"
fitness_file <- "inferred_weights_Primary_vs_Met_4T1.csv"

ptm_score_cutoff <- 0.0001

# ============================================================
# 1. Parse the musite.net PTM prediction output.
#
#    Format: a column-header line, then a FASTA description line
#    (starts with ">"), then one tab-delimited row per predicted site:
#      ID  Position  Residue  PTMscores  Cutoff=0.5
#    where PTMscores is a ";"-separated list of "PTMname:score" for
#    EVERY PTM musite scored at that site (not just the ones passing the
#    0.5 cutoff - the 5th column is the cutoff-filtered one, which we
#    don't use here). A site/PTM combo with no entry anywhere in the
#    file is assumed score = 0.
# ============================================================
raw_lines <- readLines(ptm_file)
raw_lines <- raw_lines[nzchar(raw_lines)]        # drop blank lines
raw_lines <- raw_lines[-1]                       # drop the column-header line
raw_lines <- raw_lines[!grepl("^>", raw_lines)]   # drop the FASTA description line

fields    <- strsplit(raw_lines, "\t")
positions <- as.integer(vapply(fields, `[`, character(1), 2))
ptm_field <- vapply(fields, `[`, character(1), 4)

ptm_long <- bind_rows(lapply(seq_along(positions), function(i) {
  entries <- strsplit(ptm_field[i], ";")[[1]]
  entries <- entries[nzchar(entries) & entries != "None"]
  if (length(entries) == 0) return(NULL)
  parts <- strsplit(entries, ":")
  data.frame(
    position = positions[i],
    ptm      = vapply(parts, `[`, character(1), 1),
    score    = as.numeric(vapply(parts, `[`, character(1), 2)),
    stringsAsFactors = FALSE
  )
}))

all_ptms <- sort(unique(ptm_long$ptm))
cat("PTM types found:", paste(all_ptms, collapse = ", "), "\n")

# Every site (1-378) x every PTM type, missing = 0
ptm_scores <- expand.grid(position = 1:378, ptm = all_ptms, stringsAsFactors = FALSE) %>%
  left_join(ptm_long, by = c("position", "ptm")) %>%
  mutate(score = coalesce(score, 0))

# ============================================================
# 2. Metastatic fitness (additive), 0-377 index -> 1-378 position.
# ============================================================
fitness <- read.csv(fitness_file)
names(fitness)[1] <- "position0"
fitness <- fitness %>%
  mutate(position = position0 + 1) %>%
  select(position, additive)

# ============================================================
# 3. This is the reverse of the metastatic-status plot: now, for each
#    PTM, sites are split into PTM (score >= 0.05) vs. Non-PTM
#    (score < 0.05) on the x-axis, and we compare their metastatic
#    additive fitness (y-axis).
# ============================================================
plot_df <- ptm_scores %>%
  mutate(status = ifelse(score >= ptm_score_cutoff, "PTM", "Non-PTM")) %>%
  left_join(fitness, by = "position") %>%
  filter(!is.na(additive)) %>%
  mutate(status = factor(status, levels = c("Non-PTM", "PTM")))

status_colors <- c("Non-PTM" = "grey40", "PTM" = "#d7301f")

# Two-sided Wilcoxon rank-sum test per PTM (PTM vs. Non-PTM sites).
# exact = FALSE avoids the "cannot compute exact p-value with ties" warning.
pval_df <- bind_rows(lapply(all_ptms, function(p_name) {
  sub <- plot_df %>% filter(ptm == p_name)
  test <- tryCatch(
    wilcox.test(additive ~ status, data = sub, alternative = "two.sided", exact = FALSE),
    error = function(e) NULL
  )
  pval <- if (is.null(test)) NA_real_ else test$p.value
  data.frame(
    ptm = p_name,
    p_label = if (is.na(pval)) "p = NA" else paste0("p = ", formatC(pval, format = "e", digits = 2)),
    stringsAsFactors = FALSE
  )
}))

p <- ggplot(plot_df, aes(x = status, y = additive, color = status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.5, na.rm = TRUE) +
  geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.6, na.rm = TRUE) +
  geom_text(data = pval_df, aes(x = 1.5, y = Inf, label = p_label),
            inherit.aes = FALSE, vjust = 1.4, size = 3) +
  facet_wrap(~ ptm, scales = "free_y", ncol = 4) +
  scale_color_manual(values = status_colors) +
  labs(x = NULL, y = "Metastatic fitness (additive)") +
  theme_cowplot() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90", color = NA)
  )

n_ptms <- length(all_ptms)
n_cols <- 4
n_rows <- ceiling(n_ptms / n_cols)

ggsave("mets_fitness_by_PTM_status.pdf", p,
       width = 2.5 * n_cols, height = 1.5 * n_rows, limitsize = FALSE)

cat("Wrote mets_fitness_by_PTM_status.pdf\n")
