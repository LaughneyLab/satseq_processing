rm(list=ls()) # Clear environment
setwd(getSrcDirectory(function(){})) # Set working directory to script location

library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(svglite)
library(ggplot2)
library(cowplot)
library(patchwork)

# ==============================================================================
# 1. PARAMETERS & INPUT CONFIGURATION
# ==============================================================================
# Normalization options: "none", "zscore_global", "zscore_subtract_t0", "zscore_subtract_ref"
norm_method <- "zscore_subtract_stop"

target_genes <- read.csv('STING_pathway_cGAMP_genes.csv')[,2]

rep1_tsv <- "7J9GL7_results/7J9GL7-expression-matrix.tsv"
rep1_mapping_dir <- "7J9GL7_results/7J9GL7_mapping-stats"

rep2_tsvs <- c("3F5XJK_results/3F5XJK-expression-matrix.tsv", "Q57ZYK_results/Q57ZYK-expression-matrix.tsv")
rep2_mapping_dirs <- c("3F5XJK_results/3F5XJK_mapping-stats", "Q57ZYK_results/Q57ZYK_mapping-stats")

safe_scale <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x))) else return(as.vector(scale(x)))
}

# ==============================================================================
# 2. EXTRACT AND PARSE REPLICATE 1 VALUES
# ==============================================================================
files_rep1 <- list.files(rep1_mapping_dir, pattern = "\\.tsv$", full.names = FALSE)

metadata_rep1 <- data.frame(filename = files_rep1, stringsAsFactors = FALSE) %>%
  mutate(
    sample_id = str_extract(filename, "^[A-Z0-9]+_\\d+"),
    Variant = str_split_i(filename, "_", 4),
    Time_raw = str_split_i(str_remove(filename, "\\.tsv$"), "_", 5)
  ) %>%
  filter(Time_raw != "chronic") %>% 
  mutate(
    Time = paste0(Time_raw, "h"),
    raw_colname = paste0("X", sample_id, "_cpm")
  ) %>%
  select(raw_colname, Variant, Time)

df_rep1 <- read.delim(rep1_tsv, sep = "\t",
                      stringsAsFactors = FALSE, check.names = TRUE)

rep1_processed <- df_rep1 %>% 
  filter(tolower(gene_name) %in% tolower(target_genes)) %>%
  group_by(gene_name) %>%
  summarise(across(ends_with("_cpm"), \(x) sum(x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(cols = ends_with("_cpm"), names_to = "raw_colname", values_to = "CPM") %>%
  inner_join(metadata_rep1, by = "raw_colname") %>%
  mutate(log_val = log10(1 + CPM))

# ==============================================================================
# 3. EXTRACT AND PARSE REPLICATE 2 VALUES
# ==============================================================================
files_rep2 <- unlist(lapply(rep2_mapping_dirs, function(d) list.files(d, pattern = "\\.tsv$", full.names = FALSE)))

metadata_rep2 <- data.frame(filename = files_rep2, stringsAsFactors = FALSE) %>%
  mutate(
    # Captures the single-index pattern: e.g., "3F5XJK_1" or "Q57ZYK_10"
    sample_id = str_extract(filename, "^[A-Z0-9]+_\\d+"),
    Variant = str_split_i(filename, "_", 3),
    # Captures trailing numeric digits or explicit text strings safely
    Time_raw = str_extract(str_remove(filename, "\\.tsv$"), "(chronic|\\d+)$")
  ) %>%
  filter(Time_raw != "chronic") %>% 
  mutate(
    Time = paste0(Time_raw, "h"),
    prefix = ifelse(str_detect(sample_id, "^[0-9]"), "X", ""),
    raw_colname = paste0(prefix, sample_id, "_cpm")
  ) %>%
  select(raw_colname, Variant, Time)

rep2_processed_list <- lapply(rep2_tsvs, function(p) {
  df_block <- read.delim(p, sep = "\t", stringsAsFactors = FALSE, check.names = TRUE)
  cpm_cols <- grep("_cpm$", colnames(df_block), value = TRUE)
  if(length(cpm_cols) == 0) return(NULL)
  
  df_block %>%
    filter(tolower(gene_name) %in% tolower(target_genes)) %>%
    group_by(gene_name) %>%
    summarise(across(all_of(cpm_cols), \(x) sum(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(cols = all_of(cpm_cols), names_to = "raw_colname", values_to = "CPM") %>%
    inner_join(metadata_rep2, by = "raw_colname") %>%
    mutate(log_val = log10(1 + CPM))
})

rep2_processed <- bind_rows(rep2_processed_list[!sapply(rep2_processed_list, is.null)])

# ==============================================================================
# 4. NORMALIZATION PROCESS (Controlled by norm_method string)
# ==============================================================================
apply_normalization <- function(df) {
  if (norm_method == "none") {
    df <- df %>% mutate(norm_val = log_val)
  } else if (norm_method == "zscore_global") {
    df <- df %>% 
      group_by(gene_name) %>% 
      mutate(norm_val = safe_scale(log_val)) %>% 
      ungroup()
  } else if (norm_method == "zscore_subtract_t0") {
    df <- df %>% 
      group_by(gene_name) %>% mutate(z_val = safe_scale(log_val)) %>% 
      group_by(gene_name, Variant) %>% 
      mutate(
        t0_val = mean(z_val[Time == "0h"], na.rm = TRUE),
        norm_val = z_val - ifelse(is.nan(t0_val) | is.na(t0_val), 0, t0_val)
      ) %>% ungroup()
  } else if (norm_method == "zscore_subtract_ref") {
    df <- df %>% 
      group_by(gene_name) %>% mutate(z_val = safe_scale(log_val)) %>% 
      group_by(gene_name, Time) %>% 
      mutate(
        ref_val = mean(z_val[Variant == "REF"], na.rm = TRUE),
        norm_val = z_val - ifelse(is.nan(ref_val) | is.na(ref_val), 0, ref_val)
      ) %>% ungroup()
  } else if (norm_method == "zscore_subtract_stop") {
    df <- df %>% 
      group_by(gene_name) %>% mutate(z_val = safe_scale(log_val)) %>% 
      group_by(gene_name, Time) %>% 
      mutate(
        ref_val = mean(z_val[Variant == "STOP"], na.rm = TRUE),
        norm_val = z_val - ifelse(is.nan(ref_val) | is.na(ref_val), 0, ref_val)
      ) %>% ungroup()
  }
  return(df)
}

rep1_norm <- apply_normalization(rep1_processed)
rep2_norm <- apply_normalization(rep2_processed)

# ==============================================================================
# 5. MERGE AND COMPUTE BIOLOGICAL AVERAGE
# ==============================================================================
combined_data <- bind_rows(
  rep1_norm %>% mutate(Replicate = "Rep1"),
  rep2_norm %>% mutate(Replicate = "Rep2")
) %>%
  group_by(gene_name, Variant, Time) %>%
  summarise(final_val = mean(norm_val, na.rm = TRUE), .groups = "drop") %>%
  mutate(sample_group_id = paste(Variant, Time, sep = "_"))

# ==============================================================================
# 6. DESIGN METADATA ANNOTATIONS AND COLUMNS ORDERING
# ==============================================================================
time_order <- c("0h", "2h", "6h", "12h")
custom_colors <- c("REF" = "#4E4E4E", "STOP" = "#B4B4B4",
                   "P13A" = "#FBB040", "K337R" = "#E43B3F")

metadata_final <- combined_data %>%
  select(sample_group_id, Variant, Time) %>%
  distinct() %>%
  mutate(
    Variant = factor(Variant, levels = names(custom_colors)),
    Time = factor(Time, levels = time_order)
  ) %>%
  arrange(Variant, Time)

wide_df <- combined_data %>% 
  select(gene_name, sample_group_id, final_val) %>% 
  pivot_wider(names_from = sample_group_id, values_from = final_val)

gene_names_vector <- wide_df$gene_name
wide_df$gene_name <- NULL

mat_raw <- as.matrix(wide_df)
rownames(mat_raw) <- gene_names_vector
class(mat_raw) <- "numeric"

ordered_groups <- metadata_final$sample_group_id
mat_ordered <- mat_raw[, ordered_groups, drop = FALSE]

annotation_base <- as.data.frame(metadata_final)
rownames(annotation_base) <- annotation_base$sample_group_id
annotation_base$sample_group_id <- NULL

variant_vector <- as.character(annotation_base$Variant)
gap_indices <- which(variant_vector[-1] != variant_vector[-length(variant_vector)])

time_colors <- c("0h" = "#E4D00A", "2h" = "#A1D99B", "6h" = "#74C476", "12h" = "#238443")
annotation_colors <- list(Variant = custom_colors, Time = time_colors)

heatmap_colors <- colorRampPalette(rev(brewer.pal(7, "RdBu")))(50)

if (norm_method == "none") {
  palette_breaks <- seq(0, as.numeric(quantile(abs(mat_ordered), 0.95,
                                               na.rm = TRUE)), length.out = 51)
} else {
  max_abs <- as.numeric(quantile(abs(mat_ordered), 0.95, na.rm = TRUE))
  if(max_abs == 0 || is.na(max_abs)) max_abs <- 1
  palette_breaks <- seq(-max_abs, max_abs, length.out = 51)
}

# ==============================================================================
# 7. GENERATE HEATMAP PLOT
# ==============================================================================
p_heatmap <- pheatmap(
  mat_ordered, 
  annotation_col     = annotation_base, 
  annotation_colors  = annotation_colors, 
  color              = heatmap_colors, 
  breaks             = palette_breaks, 
  cluster_rows       = TRUE, 
  cluster_cols       = FALSE, 
  gaps_col           = gap_indices, 
  scale              = "none", 
  show_colnames      = FALSE, 
  show_rownames      = TRUE, 
  border_color       = "white"
)


svglite("rnaseq_heatmap_STOP_subtracted.svg", width = 6, height = 15)
print(p_heatmap)
dev.off()


plot_variant_expression_over_time <- function(combined_data, gene_lst, ylabel) {
  
  # 1. Filter combined_data for the target gene list (case-insensitive)
  df_filtered <- combined_data %>%
    filter(tolower(gene_name) %in% tolower(gene_lst))
  
  if (nrow(df_filtered) == 0) {
    stop("None of the specified genes were found in combined_data.")
  }
  
  # 2. Compute summary metrics across the specified gene list subset per condition
  df_summary <- df_filtered %>% 
    group_by(Variant, Time) %>% 
    summarise(
      mean_val = mean(final_val, na.rm = TRUE),
      se_val   = sd(final_val, na.rm = TRUE) / sqrt(n()), 
      .groups  = "drop"
    )
  
  # 3. Enforce chronological ordering of the Time factor
  df_summary$Time <- factor(df_summary$Time, levels = c("0h", "2h", "6h", "12h"))
  
  # 4. Define requested custom color palette
  variant_colors <- c(REF = "#4E4E4E", STOP = "#B4B4B4", P13A = "#FBB040", K337R = "#E43B3F")
  
  # 5. Generate the line plot with points and error bars
  p <- ggplot(df_summary, aes(x = Time, y = mean_val, group = Variant, color = Variant)) + 
    geom_line(linewidth = 1) + 
    geom_point(size = 3) + 
    geom_linerange(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), linewidth = 0.8) + 
    scale_color_manual(values = variant_colors) + 
    labs(
      y = ylabel,
      x = "Time"
    ) + 
    theme_cowplot() + 
    theme(legend.position = "none")
  
  return(p)
}




# Read Hallmark genes and intersect with STING genes
hallmark_df <- read.csv('../../cGAMP_sc/mouse_hallmarks_cGAMP_extended.csv')
STING_genes <- read.csv('STING_pathway_cGAMP_genes.csv')[,2]

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

p_ifn <- plot_variant_expression_over_time(combined_data,ifna_genes,'Mean IFN gene TPM')

p_nfkb <- plot_variant_expression_over_time(combined_data,nfkb_genes,'Mean NFkB gene TPM')

p_upr <- plot_variant_expression_over_time(combined_data,er_genes,'Mean UPR gene TPM')

p <- p_ifn | p_nfkb | p_upr
svglite("rnaseq_pathway_score_REF_subtracted.svg", width = 12, height = 5)
print(p)
dev.off() 

plot_variant_expression_over_time_pval <- function(rep1_data, rep2_data, gene_lst, ylabel) {
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(tidyr)
  
  # 1. Combine raw replicate data sets internally
  combined_raw <- bind_rows(
    rep1_data %>% mutate(Replicate = "Rep1"),
    rep2_data %>% mutate(Replicate = "Rep2")
  )
  
  # 2. Filter for target gene(s) across all 4 time points
  df_filtered <- combined_raw %>% 
    filter(tolower(gene_name) %in% tolower(gene_lst)) %>% 
    filter(Time %in% c("0h", "2h", "6h", "12h"))
  
  if (nrow(df_filtered) == 0) {
    stop("None of the specified genes were found in the provided data sets.")
  }
  
  # 3. Compute summary metrics (mean and SE)
  df_summary <- df_filtered %>% 
    group_by(Variant, Time) %>% 
    summarise(
      mean_val = mean(norm_val, na.rm = TRUE),
      se_val = sd(norm_val, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Map Time to explicit sequential coordinates to normalize spacing gaps between panels
  df_summary <- df_summary %>%
    mutate(X_Coord = case_when(
      Time == "0h"  ~ 1,
      Time == "2h"  ~ 2,
      Time == "6h"  ~ 3,
      Time == "12h" ~ 4
    ))
  
  # 4. STATS: Run independent two-sample t-tests per time point (Skipped P13A vs K337R)
  time_points <- c("0h", "2h", "6h", "12h")
  stats_list <- list()
  
  for (tp in time_points) {
    sub_df <- df_filtered %>% filter(Time == tp)
    
    v_ref   <- sub_df$norm_val[sub_df$Variant == "REF" & !is.na(sub_df$norm_val)]
    v_p13a  <- sub_df$norm_val[sub_df$Variant == "P13A" & !is.na(sub_df$norm_val)]
    v_k337r <- sub_df$norm_val[sub_df$Variant == "K337R" & !is.na(sub_df$norm_val)]
    
    p_p13a_ref  <- NA
    p_k337r_ref <- NA
    
    if (length(v_p13a) >= 2 && length(v_ref) >= 2) {
      p_p13a_ref <- tryCatch(t.test(v_p13a, v_ref, var.equal = TRUE)$p.value, error = function(e) NA)
    }
    if (length(v_k337r) >= 2 && length(v_ref) >= 2) {
      p_k337r_ref <- tryCatch(t.test(v_k337r, v_ref, var.equal = TRUE)$p.value, error = function(e) NA)
    }
    
    stats_list[[tp]] <- data.frame(
      Time = tp,
      p_p13a_ref = p_p13a_ref,
      p_k337r_ref = p_k337r_ref,
      stringsAsFactors = FALSE
    )
  }
  stats_per_time <- bind_rows(stats_list)
  
  fmt_p <- function(p) {
    if (is.na(p)) return("p = NA")
    ifelse(p < 0.001, "p < 0.001", paste0("p = ", round(p, 4)))
  }
  
  # 5. Generate base line plot using sequential indexing scale
  df_plot <- df_summary %>% filter(Variant != "STOP")
  variant_colors <- c(REF = "#4E4E4E", P13A = "#FBB040", K337R = "#E43B3F")
  
  p <- ggplot(df_plot, aes(x = X_Coord, y = mean_val, group = Variant, color = Variant)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_linerange(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), linewidth = 0.8) +
    scale_color_manual(values = variant_colors) +
    labs(y = ylabel, x = "Time") +
    theme_cowplot() +
    theme(legend.position = "none")
  
  # 6. Draw side brackets right next to curves with maximized spacing
  time_coords <- c(1, 2, 3, 4)
  
  for (idx in seq_along(time_coords)) {
    xc <- time_coords[idx]
    ts <- time_points[idx]
    
    tp_means <- df_plot %>% filter(X_Coord == xc)
    y_ref   <- tp_means$mean_val[tp_means$Variant == "REF"]
    y_p13a  <- tp_means$mean_val[tp_means$Variant == "P13A"]
    y_k337r <- tp_means$mean_val[tp_means$Variant == "K337R"]
    
    tp_stats <- stats_per_time %>% filter(Time == ts)
    
    # Widened spacing between the remaining 2 bracket columns
    x1 <- xc + 0.20
    x2 <- xc + 0.55
    
    tick_w <- 0.04
    
    # Bracket 1: P13A vs REF
    if (length(y_p13a) > 0 && length(y_ref) > 0) {
      p <- p + 
        annotate("segment", x = x1, xend = x1, y = y_p13a, yend = y_ref, color = "black", linewidth = 0.4) +
        annotate("segment", x = x1 - tick_w, xend = x1, y = y_p13a, yend = y_p13a, color = "black") +
        annotate("segment", x = x1 - tick_w, xend = x1, y = y_ref, yend = y_ref, color = "black") +
        annotate("text", x = x1 + 0.12, y = (y_p13a + y_ref) / 2, label = fmt_p(tp_stats$p_p13a_ref), 
                 color = "black", angle = 90, size = 4.0, vjust = 0.5)
    }
    
    # Bracket 2: K337R vs REF
    if (length(y_k337r) > 0 && length(y_ref) > 0) {
      p <- p + 
        annotate("segment", x = x2, xend = x2, y = y_k337r, yend = y_ref, color = "black", linewidth = 0.4) +
        annotate("segment", x = x2 - tick_w, xend = x2, y = y_k337r, yend = y_k337r, color = "black") +
        annotate("segment", x = x2 - tick_w, xend = x2, y = y_ref, yend = y_ref, color = "black") +
        annotate("text", x = x2 + 0.12, y = (y_k337r + y_ref) / 2, label = fmt_p(tp_stats$p_k337r_ref), 
                 color = "black", angle = 90, size = 4.0, vjust = 0.5)
    }
  }
  
  # 7. Apply categorical labels to coordinates and adjust right-side limit for 2 columns
  p <- p + scale_x_continuous(
    breaks = c(1, 2, 3, 4), 
    labels = c("0h", "2h", "6h", "12h"),
    limits = c(0.7, 5.0)
  )
  
  return(p)
}


p_ifn <- plot_variant_expression_over_time_pval(rep1_norm, rep2_norm,ifna_genes,
                                                'LogFC IFN gene TPM')

p_nfkb <- plot_variant_expression_over_time_pval(rep1_norm, rep2_norm,nfkb_genes,
                                                 'LogFC NFkB gene TPM')

p_upr <- plot_variant_expression_over_time_pval(rep1_norm, rep2_norm,er_genes,
                                                'LogFC UPR gene TPM')

p <- p_ifn | p_nfkb | p_upr
svglite("rnaseq_pathway_score_REF_subtracted_pval.svg", width = 12, height = 5)
print(p)
dev.off() 

# Save plot
#ggsave(plot=p, "rnaseq_pathway_score_REF_subtracted.pdf",
#       width=12, height=5)
  
# # Hotspot module genes
# hotspot_genes <- read.csv('../tumor_sc_invivo/hotspot_gene_modules.csv')
# 
# prolif_genes_hotspot <- na.omit(hotspot_genes$Gene[hotspot_genes$Module == 1])
# ifn_genes_hotspot    <- na.omit(hotspot_genes$Gene[hotspot_genes$Module == 2])
# nfkb_genes_hotspot   <- na.omit(hotspot_genes$Gene[hotspot_genes$Module == 3])
# 
# p_ifn_hotspot <- plot_variant_expression_over_time(ifn_genes_hotspot,
#                                                    'Median IFN gene TPM')
# 
# p_nfkb_hotspot <- plot_variant_expression_over_time(nfkb_genes_hotspot,
#                                                     'Median NFkB gene TPM')
# 
# p_prolif_hotspot <- plot_variant_expression_over_time(prolif_genes_hotspot,
#                                                       'Median Proliferation gene TPM')
# 
# p <- p_prolif_hotspot | p_ifn_hotspot | p_nfkb_hotspot
# svglite("rnaseq_hotspot_score_median.svg", width = 12, height = 5)
# print(p)
# dev.off() 
# 
# 
# gene_list <- c(
#   'Htra1', 'Rnd3', 'Ptgs2', 'Tiparp', 'Fam107b', 'Wdr36', 'Rel', 'Pcna', 'Arih1', 'Nfkb1', 
#   'Asap1', 'Fn1', 'Psat1', 'Plk4', 'Msi2', 'Diaph3', 'Psip1', 'Ankrd1', 'Tra2a', 'Elf1', 
#   'Ehd4', 'Samhd1', 'Slfn9', 'Dtx3l', 'Rnf213', 'Ifit3', 'Stat1', 'Ifi204', 'Ifrd1', 'H3f3b', 
#   'Ltbp1', 'Cdt1', '1110038B12Rik', 'Gas5', 'Zfas1', 'Skil', 'Snhg1', 'Cdh1', 'Rsbn1', 'Sox9', 
#   'Rhob', 'Zfp36l2', 'Jun', 'Fmnl2', 'Kmt2c', 'Trps1', 'Mllt3', 'Arl6ip1', 'Pcyt1b', 'Atxn7l1', 
#   'Brip1', 'Diaph1', 'Diaph2', 'Rlf', 'Adamts6', 'Malat1', 'Dst', 'Chd2', 'Rbms1', 'Mast4', 
#   'Lyn', 'Atxn1', 'Myo10', 'Vegfc', 'Cxcl10', 'Irf1', 'Atf3', 'Ifi202b', 'Ifi203', 'Mndal', 
#   'Stat2', 'Parp14', 'Parp9', 'Gbp7', 'Ddx60', 'Oasl2', 'Isg15', 'Ifitm3', 'Ifi44', 'Sp100', 
#   'Xaf1', 'Psmd8', 'Npc2', 'Itm2b', 'Psme1', 'Psme2', 'Ly6e', 'H2-K1', 'Tapbp', 'Psmb8', 
#   'H2-D1', 'Psmb10', 'Gbp2', 'Gbp2b', 'Trim30a', 'Slfn2', 'Nmi', 'Gbp3', 'Herc6', 'Ghr', 
#   'Hspa4l', 'Nlrc5', 'Sik3', 'Ier5', 'Tmsb4x', 'Glis3', '2410006H16Rik', 'Slc38a2', 'Atf4', 
#   'Akt3', 'Glrp1', 'Spag9', 'Bbx', 'Utrn', 'Myof', 'Zbtb38', 'Sqstm1', 'Sat1', 'Ppp1r15a', 
#   'Gadd45a', 'Pfkp', 'Nupr1', 'Cdk14', 'Rsrp1', 'Neat1', 'Arhgap6', 'Hivep2', 'Fnip1', 'Ubc', 
#   'Tanc2', 'Rnf19b', 'Abtb2', 'Nfat5', 'Flnb', 'Phf21a', 'Hk2', 'Wee1')
# 
# 
# # --- Generate and save the LANDSCAPE layout ---
# p_landscape <- plot_gene_pheatmap_landscape(tsv_path = tsv_matrix_path,
#                                             gene_list)
# svglite("rnaseq_heatmap_selected_genes.svg", width = 20, height = 6)
# print(p_landscape)
# dev.off()
# 
# p_ifn_hotspot <- plot_variant_expression_over_time(intersect(ifn_genes_hotspot,gene_list),
#                                                    'Median IFN gene TPM')
# 
# p_nfkb_hotspot <- plot_variant_expression_over_time(intersect(nfkb_genes_hotspot,gene_list),
#                                                     'Median NFkB gene TPM')
# 
# p_prolif_hotspot <- plot_variant_expression_over_time(intersect(prolif_genes_hotspot,gene_list),
#                                                       'Median Proliferation gene TPM')
# 
# p <- p_prolif_hotspot | p_ifn_hotspot | p_nfkb_hotspot
# svglite("rnaseq_hotspot_score_median_selected_genes.svg", width = 12, height = 5)
# print(p)
# dev.off()
# 
# 



  