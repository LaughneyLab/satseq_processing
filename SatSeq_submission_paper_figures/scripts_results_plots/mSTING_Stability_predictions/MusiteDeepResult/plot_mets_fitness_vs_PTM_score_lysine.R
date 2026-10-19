rm(list=ls()) # Clear environment
setwd(getSrcDirectory(function(){})) # set working directory to where this script is

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(cowplot)

ptm_file     <- "Prediction_results.txt"
fitness_file <- "inferred_weights_Primary_vs_Met_4T1.csv"

# ============================================================
# 1. Parse the musite.net PTM prediction output (same approach as the
#    other PTM scripts): a column-header line, a FASTA description line
#    (starts with ">"), then one tab-delimited row per predicted site:
#      ID  Position  Residue  PTMscores  Cutoff=0.5
#    PTMscores is a ";"-separated list of "PTMname:score" for EVERY PTM
#    musite scored at that site (not just the >=0.5 cutoff column).
# ============================================================
raw_lines <- readLines(ptm_file)
raw_lines <- raw_lines[nzchar(raw_lines)]        # drop blank lines
raw_lines <- raw_lines[-1]                       # drop the column-header line
raw_lines <- raw_lines[!grepl("^>", raw_lines)]   # drop the FASTA description line

fields    <- strsplit(raw_lines, "\t")
positions <- as.integer(vapply(fields, `[`, character(1), 2))
residues  <- vapply(fields, `[`, character(1), 3)
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

# ============================================================
# 2. Restrict to lysine (K) sites only. Only PTM types that actually
#    have entries at K sites in the file end up as facets (in practice:
#    Ubiquitination, SUMOylation, N6-acetyllysine, Methyllysine,
#    Hydroxylysine) - other PTM types (e.g. phosphorylation) never
#    apply to lysine and would just be all-zero, so they're dropped
#    rather than plotted as a degenerate facet.
# ============================================================
lysine_positions <- sort(unique(positions[residues == "K"]))
cat("Lysine sites found:", length(lysine_positions), "\n")

lysine_ptm_long <- ptm_long %>% filter(position %in% lysine_positions)
lysine_ptms <- sort(unique(lysine_ptm_long$ptm))
cat("PTM types found at lysines:", paste(lysine_ptms, collapse = ", "), "\n")

# Every lysine site x every lysine-associated PTM type, missing = 0
ptm_scores_lys <- expand.grid(position = lysine_positions, ptm = lysine_ptms,
                               stringsAsFactors = FALSE) %>%
  left_join(lysine_ptm_long, by = c("position", "ptm")) %>%
  mutate(score = coalesce(score, 0))

# ============================================================
# 3. Metastatic fitness (additive), 0-377 index -> 1-378 position.
# ============================================================
fitness <- read.csv(fitness_file)
names(fitness)[1] <- "position0"
fitness <- fitness %>%
  mutate(position = position0 + 1) %>%
  select(position, additive)

# ============================================================
# 4. This is the reverse axis layout from the earlier PTM-vs-mets
#    boxplot: now PTM score is continuous on the x-axis (not
#    binarized), metastatic additive fitness is on the y-axis, and
#    each lysine site is one point, labeled by its WT+position
#    (e.g. "K337") with ggrepel so labels don't overlap.
# ============================================================
plot_df <- ptm_scores_lys %>%
  left_join(fitness, by = "position") %>%
  filter(!is.na(additive)) %>%
  mutate(site_label = paste0("K", position))

# stat_correlation()'s built-in combined label (via use_label()) is meant
# to be *parsed* as a plotmath expression - asking for output.type="text"
# does not reformat it, so it renders literally as e.g.
# `R = 0.86*", "*P < "0.001"` instead of clean text. Pulling the raw
# numeric `cor`/`p.value` (output.type="numeric") and building the label
# ourselves avoids that entirely and gives R (not R^2) plus a plain p-value.
p <- ggplot(plot_df, aes(x = score, y = additive)) +
  geom_point(color = "grey30", alpha = 0.8, size = 1.8, na.rm = TRUE) +
  stat_poly_line(method = "lm", color = "#2166ac", fill = "#2166ac",
                 alpha = 0.15, linewidth = 0.7, na.rm = TRUE) +
  stat_correlation(
    method = "pearson",
    mapping = aes(label = after_stat(
      paste0("R = ", round(cor, 2), ", p = ", formatC(p.value, format = "e", digits = 2))
    )),
    output.type = "numeric",
    label.x = "right", label.y = "bottom", size = 3, na.rm = TRUE
  ) +
  geom_text_repel(aes(label = site_label), size = 2.8, max.overlaps = Inf,
                   segment.size = 0.2, segment.color = "grey50") +
  facet_wrap(~ ptm, scales = "free_x", ncol = 3) +
  labs(x = "PTM score", y = "Metastatic fitness") +
  theme_cowplot() +
  theme(strip.background = element_rect(fill = "grey90", color = NA))

n_facets <- length(lysine_ptms)
n_cols <- 3
n_rows <- ceiling(n_facets / n_cols)

ggsave("mets_fitness_vs_PTM_score_lysine.pdf", p,
       width = 3.5 * n_cols, height = 2.5 * n_rows, limitsize = FALSE)

cat("Wrote mets_fitness_vs_PTM_score_lysine.pdf\n")
