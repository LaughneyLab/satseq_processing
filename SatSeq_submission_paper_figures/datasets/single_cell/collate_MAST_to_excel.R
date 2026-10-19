rm(list = ls())

# ============================================================
#  collate_MAST_to_excel.R
#
#  Run after all 33 SLURM array tasks have completed.
#  Reads CSVs from MAST_output/raw/ and writes one xlsx per
#  variant comparison (3 files per cell type):
#
#    <CellType>_P13_vs_REF.xlsx
#    <CellType>_K337_vs_REF.xlsx
#    <CellType>_STOP_vs_REF.xlsx
#
#  Each file has 2 columns:
#    gene | MAST_score_<variant>_vs_REF
#
#  All scores = tanh(k * raw), bounded in (-1, 1).
#    sign: raw = -log10(FDR) * sign(logFC)
#    NES:  raw = -log10(FDR) * logFC
#
#  Usage:
#    Rscript --vanilla collate_MAST_to_excel.R [raw_dir] [out_dir]
# ============================================================

suppressPackageStartupMessages(library(writexl))  # install.packages("writexl")

args    <- commandArgs(trailingOnly = TRUE)
RAW_DIR <- if (length(args) >= 1) args[1] else "MAST_output/raw"
OUT_DIR <- if (length(args) >= 2) args[2] else "MAST_output/excel"

if (!dir.exists(RAW_DIR)) stop("Raw CSV directory not found: ", RAW_DIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VARIANTS <- c("P13", "K337", "STOP")

# All scores use tanh(k * raw) to compress into (-1, 1).
#
# sign scores: raw = -log10(FDR) * sign(logFC)
#   k=0.01 -> variance-matched to gene_test example file (var~0.002, SD~0.045)
#             Solved numerically: brentq on var(tanh(k*raw)) == mean(var(gene_test cols))
#             Range across gene_test columns: k=0.009 to k=0.012
#   k=0.1  -> 10x broader; useful if cell type has weaker signal
#
# NES scores: raw = -log10(FDR) * logFC
#   logFC amplifies the raw value, so k must be smaller:
#   k=0.3 -> highlights only the very top hits
#   k=1.0 -> spreads the distribution more

tanh_score <- function(fdr, coef, type = c("sign", "NES"), k) {
  type <- match.arg(type)
  raw  <- if (type == "sign") -log10(pmax(fdr, 1e-300)) * sign(coef)
           else               -log10(pmax(fdr, 1e-300)) * coef
  tanh(k * raw)
}

# ---- Discover CSVs -----------------------------------------------------------
csv_files <- list.files(RAW_DIR, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("No CSV files found in: ", RAW_DIR)
cat(sprintf("Found %d CSV file(s)\n\n", length(csv_files)))

parse_filename <- function(f) {
  base <- sub("\\.csv$", "", basename(f))
  for (var in VARIANTS) {
    suffix <- paste0("_", var, "_vs_REF")
    if (endsWith(base, suffix))
      return(list(ct_safe = sub(paste0(suffix, "$"), "", base), variant = var))
  }
  NULL
}

parsed <- Filter(Negate(is.null), lapply(csv_files, function(f) {
  p <- parse_filename(f); if (!is.null(p)) c(p, list(path = f))
}))

ct_safes <- unique(sapply(parsed, `[[`, "ct_safe"))
cat(sprintf("Cell types detected: %d\n", length(ct_safes)))
cat(paste(" -", ct_safes, collapse = "\n"), "\n\n")

# ---- One xlsx per cell type per variant --------------------------------------
for (ct in ct_safes) {
  for (var in VARIANTS) {

    row <- Filter(function(p) p$ct_safe == ct && p$variant == var, parsed)
    if (length(row) == 0) {
      cat(sprintf("  WARNING: no CSV for %s %s – skipping\n", ct, var))
      next
    }

    dt <- tryCatch(read.csv(row[[1]]$path, stringsAsFactors = FALSE),
                   error = function(e) { cat("  ERROR:", e$message, "\n"); NULL })
    if (is.null(dt) || nrow(dt) == 0) {
      cat(sprintf("  WARNING: empty CSV for %s %s – skipping\n", ct, var))
      next
    }

    missing <- setdiff(c("primerid", "fdr", "coef"), names(dt))
    if (length(missing) > 0) {
      cat(sprintf("  WARNING: missing cols (%s) for %s %s – skipping\n",
                  paste(missing, collapse = ", "), ct, var))
      next
    }

    score_col <- paste0("MAST_score_", var, "_vs_REF")
    df <- data.frame(
      gene             = dt$primerid,
      score            = tanh_score(dt$fdr, dt$coef, "sign", k = 0.01),
      stringsAsFactors = FALSE
    )
    names(df)[2] <- score_col

    out_xlsx <- file.path(OUT_DIR, paste0(ct, "_", var, "_vs_REF.xlsx"))
    write_xlsx(df, out_xlsx)
    cat(sprintf("  Saved: %s  (%d genes)\n", out_xlsx, nrow(df)))
  }
  cat("\n")
}

cat("=== Done ===\n")
cat("Files are in:", OUT_DIR, "\n")
