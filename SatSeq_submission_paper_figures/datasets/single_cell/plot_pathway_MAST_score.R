rm(list=ls()) # Clear environment

library(dplyr)
library(ggplot2)
library(cowplot)

# ============================================================
# Config
# ============================================================
gmt_file <- "mouse_hallmarks.gmt"
mast_dir <- "MAST_output_lung_mets"   # output of run_MAST_variant_vs_REF_lung_mets.R

celltype_levels <- c("4T1", "AM", "T cells")
variant_levels  <- c("STOP", "P13", "K337")

# Same hex codes as sc_analysis_lung_mets.py's variant_colors dict
variant_colors <- c(
  "P13"  = "#FBB040",
  "K337" = "#E43B3F",
  "STOP" = "#B4B4B4"
)
variant_labels <- c(
  "STOP" = "STOP vs REF",
  "P13"  = "P13 vs REF",    # rename to "P13A vs REF" if you know the exact
  "K337" = "K337 vs REF"    # substituted residue, e.g. "K337R vs REF"
)

# All 6 shown in every celltype panel, in this order (per your call — no
# per-celltype filtering, no significance testing/asterisks this round)
pathway_display_names <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE" = "Interferon Alpha",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" = "Interferon Gamma",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"   = "TNFA via NFkB",
  "HALLMARK_MYC_TARGETS_V1"            = "MYC Targets v1",
  "HALLMARK_MYC_TARGETS_V2"            = "MYC Targets v2",
  "HALLMARK_IL2_STAT5_SIGNALING"       = "IL2 STAT5"
)
pathway_order <- unname(pathway_display_names)

# ============================================================
# 1. Read the GMT file, keep only the 6 pathways of interest
# ============================================================
read_gmt <- function(file) {
  lines <- readLines(file)
  lines <- lines[nzchar(lines)]
  split_lines <- strsplit(lines, "\t")
  genesets <- lapply(split_lines, function(x) x[-c(1, 2)])
  names(genesets) <- vapply(split_lines, `[`, character(1), 1)
  genesets
}

all_pathways <- read_gmt(gmt_file)
missing_pw <- setdiff(names(pathway_display_names), names(all_pathways))
if (length(missing_pw) > 0) {
  stop("Pathway(s) not found in ", gmt_file, ": ", paste(missing_pw, collapse = ", "))
}
pathways <- all_pathways[names(pathway_display_names)]

# ============================================================
# 2. Read every MAST fcHurdle CSV and combine
# ============================================================
mast_files <- list.files(mast_dir, pattern = "^fcHurdle_.*\\.csv$",
                          recursive = TRUE, full.names = TRUE)
if (length(mast_files) == 0) {
  stop("No fcHurdle_*.csv files found under ", mast_dir)
}

mast_all <- bind_rows(lapply(mast_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  # `variant`/`celltype` columns are written by run_MAST_variant_vs_REF_lung_mets.R;
  # fall back to parsing the filename if an edited/older CSV is missing them
  if (!"variant" %in% names(df) || !"celltype" %in% names(df)) {
    base  <- sub("\\.csv$", "", basename(f))
    parts <- strsplit(sub("^fcHurdle_", "", base), "_")[[1]]
    df$variant  <- parts[length(parts)]
    df$celltype <- paste(parts[-length(parts)], collapse = "_")
  }
  df
}))

# Per-gene MAST score: -log10(FDR) * sign(coef) — same convention as your
# earlier heatmap script's `score <- -log10(fdr) * sign(coef)`
mast_all <- mast_all %>%
  mutate(score = -log10(fdr + 1e-300) * sign(coef))

# ============================================================
# 3. Average score per pathway / variant / celltype, over genes that were
#    actually tested in that MAST run (i.e. present in that fcHurdle file)
# ============================================================
pathway_scores <- bind_rows(lapply(names(pathways), function(pid) {
  genes <- pathways[[pid]]
  mast_all %>%
    filter(primerid %in% genes) %>%
    group_by(celltype, variant) %>%
    summarise(mean_score = mean(score, na.rm = TRUE),
              n_genes = n(), .groups = "drop") %>%
    mutate(pathway = pathway_display_names[[pid]])
}))

pathway_scores <- pathway_scores %>%
  filter(variant %in% names(variant_colors)) %>%
  mutate(
    celltype = factor(celltype, levels = celltype_levels),
    variant  = factor(variant, levels = variant_levels),
    pathway  = factor(pathway, levels = pathway_order)
  )

# ============================================================
# 4. Plot
# ============================================================
p <- ggplot(pathway_scores, aes(x = pathway, y = mean_score, fill = variant)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, na.rm = TRUE) +
  facet_wrap(~ celltype, nrow = 1) +
  scale_fill_manual(values = variant_colors, labels = variant_labels, name = NULL) +
  labs(x = NULL, y = "Pathway MAST Score") +
  theme_cowplot() +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = c(0.70, 0.3),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("pathway_MAST_score_by_variant_celltype.pdf", p,
 width = 6, height = 4, device = "pdf")
