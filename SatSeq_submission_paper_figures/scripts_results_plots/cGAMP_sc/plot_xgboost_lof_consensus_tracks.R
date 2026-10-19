rm(list=ls()) # Clear environment
setwd(getSrcDirectory(function(){})) # set working directory to where this script is

library(dplyr)
library(tidyr)
library(pROC)
library(ggplot2)
library(svglite)

lof_df <- read.csv("LOF_variants.csv")

# Trapezoidal-rule AUC, used below for a quick AUPRC readout
auc_trapz <- function(x, y) {
  ord <- order(x)
  x <- x[ord]; y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# ============================================================
# 1. For one cell-line's XGBoost output: build the same predictor as
#    analyze_xgboost_scores.R (frac_0 = count_label_0 / total), attach
#    the ground-truth is_lof from LOF_variants.csv, and pick a threshold
#    off the PRECISION-RECALL curve rather than the ROC curve.
#
#    LOF is rare (<10% of residues), so ROC/Youden thresholding (which
#    weighs specificity equally with sensitivity) tends to pick a cutoff
#    that looks great on the majority class but is too lenient for a rare
#    positive class. Instead we scan every threshold on the PR curve and
#    keep the one that maximizes F1 = 2*precision*recall/(precision+recall),
#    which is the standard "optimal threshold" choice for imbalanced/rare-
#    positive problems.
# ============================================================
classify_dataset <- function(file, dataset_label) {
  probs <- read.csv(file)
  probs$frac_0 <- probs$count_label_0 / (probs$count_label_1 + probs$count_label_0)

  probs <- probs %>%
    left_join(lof_df, by = "variant") %>%
    mutate(is_lof = ifelse(is.na(group), 0, 1))

  roc_obj <- roc(response = probs$is_lof, predictor = probs$frac_0, quiet = TRUE)

  pr <- coords(roc_obj, "all", ret = c("threshold", "precision", "recall"),
               transpose = FALSE)
  pr <- pr[is.finite(pr$threshold), ]  # drop the -Inf/Inf boundary rows
  auprc <- auc_trapz(pr$recall, pr$precision)
  cat(sprintf("[%s] AUC-ROC = %.4f, AUPRC = %.4f (n = %d variants, %d true LOF, %.1f%% positive)\n",
              dataset_label, as.numeric(auc(roc_obj)), auprc,
              nrow(probs), sum(probs$is_lof), 100 * mean(probs$is_lof)))

  pr$f1 <- 2 * pr$precision * pr$recall / (pr$precision + pr$recall)
  best <- pr[which.max(pr$f1), , drop = FALSE]
  cat(sprintf("[%s] optimal threshold (frac_0, max F1) = %.4f  (precision = %.3f, recall = %.3f, F1 = %.3f)\n",
              dataset_label, best$threshold, best$precision, best$recall, best$f1))

  # roc_obj$direction: "<" means controls (is_lof=0) < cases (is_lof=1), i.e. higher
  # frac_0 -> more likely LOF; ">" is the reverse. Applying whichever direction pROC
  # found keeps the binary call consistent with how precision/recall were computed.
  if (roc_obj$direction == "<") {
    predicted_lof <- probs$frac_0 > best$threshold
  } else {
    predicted_lof <- probs$frac_0 < best$threshold
  }

  # Position = leading WT-residue letter(s) stripped from the variant string
  # (e.g. "S365" -> 365). Non-positional variants (dCTT, RRDD, ...) -> NA and
  # are kept (not dropped) - they get plotted as extra columns after 378.
  position <- suppressWarnings(as.numeric(gsub("^[A-Za-z]+", "", probs$variant)))

  data.frame(variant = probs$variant, position = position, predicted_lof = predicted_lof,
             n_cells = probs$count_label_1 + probs$count_label_0,
             stringsAsFactors = FALSE)
}

lof_4T1   <- classify_dataset("XGBoost_output_4T1_cGAMP.csv",   "4T1")
lof_EO771 <- classify_dataset("XGBoost_output_EO771_cGAMP.csv", "EO771")

# ============================================================
# 2. One row per STING "column": positions 1-378 in order, followed by
#    any non-positional variants seen in either XGBoost file (e.g. dCTT,
#    RRDD), each column carrying all 5 binary calls plus a consensus call
#    (TRUE where XGBoost calls LOF in BOTH 4T1 and EO771, AND the residue
#    is annotated LOF in at least one of IFN/NFkB/UPR).
# ============================================================
main_df <- data.frame(x = 1:378, label = as.character(1:378), stringsAsFactors = FALSE)

extra_labels <- unique(c(
  lof_4T1$variant[is.na(lof_4T1$position)],
  lof_EO771$variant[is.na(lof_EO771$position)]
))
# Put known special variants in a familiar order; anything unrecognized is
# appended after, in whatever order it was first encountered.
preferred_order <- c("dCTT", "RRDD", "STOP")
ord <- match(extra_labels, preferred_order)
extra_labels <- extra_labels[order(ifelse(is.na(ord), Inf, ord))]

# Leave a visual gap between position 378 and the non-positional columns,
# and space the extra columns themselves out a bit so their (rotated) text
# labels don't collide
extra_gap   <- 8
extra_step  <- 4
extra_x_start <- 378 + extra_gap
extra_df <- data.frame(
  x = seq(extra_x_start, by = extra_step, length.out = length(extra_labels)),
  label = extra_labels,
  stringsAsFactors = FALSE
)

full_df <- bind_rows(main_df, extra_df)

# LOF pathway flags: match by residue number for the 1-378 part, and by
# variant name for anything else (covers dCTT/RRDD if they were ever
# annotated in LOF_variants.csv)
lof_lookup <- bind_rows(
  lof_df %>% transmute(label = as.character(residue), group),
  lof_df %>% transmute(label = variant, group)
) %>% distinct(label, .keep_all = TRUE)

# dCTT (deleted C-terminal tail) and RRDD are known to be LOF across all 3
# pathways, but aren't annotated in LOF_variants.csv - hardcode them here
# rather than relying on the (empty) lookup above.
manual_triple_lof <- c("dCTT", "RRDD")

full_df <- full_df %>%
  left_join(lof_lookup, by = "label") %>%
  mutate(group = coalesce(group, "")) %>%
  mutate(
    UPR  = grepl("UPR",  group) | label %in% manual_triple_lof,
    IFN  = grepl("IFNa", group) | label %in% manual_triple_lof,
    NFkB = grepl("NFkB", group) | label %in% manual_triple_lof
  ) %>%
  select(-group)

# XGBoost calls: match each dataset's variant/position directly to `label`
xg_4T1 <- lof_4T1 %>%
  transmute(label = ifelse(is.na(position), variant, as.character(position)),
            XGBoost_4T1 = predicted_lof, Cells_4T1 = n_cells)
xg_EO771 <- lof_EO771 %>%
  transmute(label = ifelse(is.na(position), variant, as.character(position)),
            XGBoost_EO771 = predicted_lof, Cells_EO771 = n_cells)

full_df <- full_df %>%
  left_join(xg_4T1, by = "label") %>%
  left_join(xg_EO771, by = "label") %>%
  mutate(
    XGBoost_4T1   = coalesce(XGBoost_4T1,   FALSE),
    XGBoost_EO771 = coalesce(XGBoost_EO771, FALSE),
    # Consensus: XGBoost calls LOF in BOTH cell lines, AND the residue is
    # annotated LOF in at least one of the 3 pathways (not all 3)
    consensus     = XGBoost_4T1 & XGBoost_EO771 & (UPR | IFN | NFkB)
    # Cells_4T1 / Cells_EO771 (count_label_1 + count_label_0) are left as NA
    # where a position wasn't in that dataset, so no bar is drawn there.
  ) %>%
  arrange(x)

cat("Columns where XGBoost agrees in both cell lines AND at least one pathway LOF (consensus):",
    paste(full_df$label[full_df$consensus], collapse = ", "), "\n")

# ============================================================
# 3. Plot: 6 binary tracks (tiles) + 2 cell-count tracks (bars), plain
#    ggplot2 (no heatmap package) - positions 1-378 followed by dCTT/RRDD
#    (or whatever non-positional variants were found). The two geoms
#    (tile for binary calls, col for cell counts) are stacked via
#    facet_grid on a shared 8-level row factor, so everything stays on
#    one aligned x-axis.
# ============================================================
track_display <- c(
  UPR = "UPR LOF", IFN = "IFN LOF", NFkB = "NFkB LOF",
  XGBoost_4T1 = "XGBoost 4T1", XGBoost_EO771 = "XGBoost EO771",
  consensus = "Consensus"
)
# top-to-bottom order in the plot - each XGBoost LOF-call track is
# immediately followed by the cell-count track backing it
track_order <- c("UPR LOF", "IFN LOF", "NFkB LOF",
                  "XGBoost 4T1", "Cells 4T1",
                  "XGBoost EO771", "Cells EO771",
                  "Consensus")

binary_df <- full_df %>%
  select(x, label, UPR, IFN, NFkB, XGBoost_4T1, XGBoost_EO771, consensus) %>%
  pivot_longer(cols = c(UPR, IFN, NFkB, XGBoost_4T1, XGBoost_EO771, consensus),
               names_to = "track", values_to = "value") %>%
  mutate(
    track = factor(track_display[track], levels = track_order),
    fill_key = paste0(track, "_", value)
  )

count_df <- full_df %>%
  select(x, label, Cells_4T1, Cells_EO771) %>%
  pivot_longer(cols = c(Cells_4T1, Cells_EO771), names_to = "track", values_to = "n_cells") %>%
  filter(!is.na(n_cells)) %>%  # no bar where that position wasn't in the dataset
  mutate(
    track = recode(track, Cells_4T1 = "Cells 4T1", Cells_EO771 = "Cells EO771"),
    track = factor(track, levels = track_order)
  )

fill_values <- c(
  "UPR LOF_TRUE" = "#35b6e3",  "UPR LOF_FALSE" = "white",
  "IFN LOF_TRUE" = "#e3312a",  "IFN LOF_FALSE" = "white",
  "NFkB LOF_TRUE" = "#7ac943", "NFkB LOF_FALSE" = "white",
  "XGBoost 4T1_TRUE" = "#2166ac",   "XGBoost 4T1_FALSE" = "white",
  "XGBoost EO771_TRUE" = "#2166ac", "XGBoost EO771_FALSE" = "white",
  "Consensus_TRUE" = "black",  "Consensus_FALSE" = "white"
)
cell_bar_color <- "grey35"

label_every <- 25
main_breaks <- c(1, seq(label_every, 378, by = label_every))
x_breaks <- c(main_breaks, extra_df$x)
x_labels <- c(as.character(main_breaks), extra_df$label)

separator_x <- 378 + extra_gap / 2

p <- ggplot() +
  geom_tile(data = binary_df, aes(x = x, y = 1, fill = fill_key),
            width = 1, height = 0.8, color = NA) +
  geom_col(data = count_df, aes(x = x, y = n_cells),
           width = 1, fill = cell_bar_color) +
  geom_vline(xintercept = separator_x, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  scale_fill_manual(values = fill_values, guide = "none") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels,
                      expand = expansion(add = 1)) +
  facet_grid(rows = vars(track), scales = "free_y", switch = "y") +
  labs(x = "STING position", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(color = "grey50"),
    panel.spacing.y = unit(0.8, "mm"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 10),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 75)
  )

ggsave("xgboost_LOF_consensus_tracks.svg",
       p, width = 14, height = 4.5, device = svglite::svglite)

cat("Wrote xgboost_LOF_consensus_tracks.svg\n")
