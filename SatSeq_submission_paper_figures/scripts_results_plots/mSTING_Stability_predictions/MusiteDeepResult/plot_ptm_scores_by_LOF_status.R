rm(list=ls()) # Clear environment
setwd(getSrcDirectory(function(){})) # set working directory to where this script is

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

ptm_file <- "Prediction_results.txt"
lof_file <- "LOF_variants.csv"

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
# 2. LOF residues (already 1-378) vs. the rest of 1-378.
# ============================================================
lof_df <- read.csv(lof_file)

lof_status <- data.frame(position = 1:378) %>%
  mutate(status = ifelse(position %in% lof_df$residue, "LOF", "Non-LOF"))
cat("LOF sites:", sum(lof_status$status == "LOF"),
    "/ Non-LOF:", sum(lof_status$status == "Non-LOF"), "\n")

# ============================================================
# 3. Merge and plot: one facet per PTM, boxplot + points of PTM score
#    for LOF vs. non-LOF residues
# ============================================================
plot_df <- ptm_scores %>%
  left_join(lof_status, by = "position") %>%
  filter(!is.na(status)) %>%
  mutate(status = factor(status, levels = c("Non-LOF", "LOF")))

status_colors <- c("Non-LOF" = "grey40", "LOF" = "#d7301f")

# Two-sided Wilcoxon rank-sum test per PTM (LOF vs. Non-LOF). exact = FALSE
# avoids the "cannot compute exact p-value with ties" warning, which is
# expected here since most sites score exactly 0 for any given PTM.
pval_df <- bind_rows(lapply(all_ptms, function(p_name) {
  sub <- plot_df %>% filter(ptm == p_name)
  test <- tryCatch(
    wilcox.test(score ~ status, data = sub, alternative = "two.sided", exact = FALSE),
    error = function(e) NULL
  )
  pval <- if (is.null(test)) NA_real_ else test$p.value
  data.frame(
    ptm = p_name,
    p_label = if (is.na(pval)) "p = NA" else paste0("p = ", formatC(pval, format = "e", digits = 2)),
    stringsAsFactors = FALSE
  )
}))

p <- ggplot(plot_df, aes(x = status, y = score, color = status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.5, na.rm = TRUE) +
  geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.6, na.rm = TRUE) +
  geom_text(data = pval_df, aes(x = 1.5, y = Inf, label = p_label),
            inherit.aes = FALSE, vjust = 1.4, size = 3) +
  facet_wrap(~ ptm, scales = "free_y", ncol = 4) +
  scale_color_manual(values = status_colors) +
  labs(x = NULL, y = "PTM score") +
  theme_cowplot() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90", color = NA)
  )

n_ptms <- length(all_ptms)
n_cols <- 4
n_rows <- ceiling(n_ptms / n_cols)

ggsave("ptm_scores_by_LOF_status.pdf", p,
       width = 2.5 * n_cols, height = 1.5 * n_rows, limitsize = FALSE)

cat("Wrote ptm_scores_by_LOF_status.pdf\n")
