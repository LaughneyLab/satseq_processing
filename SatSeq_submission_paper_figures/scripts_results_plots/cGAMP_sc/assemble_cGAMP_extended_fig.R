rm(list=ls()) # Clear environment

#################
## libraries ####
#################
library(ggplot2)
library(patchwork)
library(magick)
library(grid)
library(ComplexHeatmap)
library(svglite)


### Gather and process data ####
setwd(getSrcDirectory(function(){})) # set working directory to where this script is

# =========================================================================
# 1. Load Figures & Base Data
# =========================================================================
p_scPC_4T1     <- wrap_plots(readRDS("PC1-10_pathways_heatmap_4T1.rds"))
p_scPC_EO771   <- wrap_plots(readRDS("PC1-10_pathways_heatmap_EO771.rds"))
p_elbow        <- wrap_plots(readRDS("MAST_cumsum_IFN_colored_4T1_EO771.rds"))
p_corr_pathway <- wrap_plots(readRDS("sc_scores_cell_lines_FACS_colored.rds"))
p_sc_FACS      <- wrap_plots(readRDS("sc_vs_FACS.rds"))

p_umap_4T1   <- image_read("umaps_4T1_cGAMP_LOF_colored.png") 
p_umap_4T1   <- wrap_elements(rasterGrob(p_umap_4T1)) + theme(plot.margin = margin(0,0,0,0,"pt"))

p_umap_EO771 <- image_read("umaps_EO771_cGAMP_LOF_colored.png") 
p_umap_EO771 <- wrap_elements(rasterGrob(p_umap_EO771)) + theme(plot.margin = margin(0,0,0,0,"pt"))

# Load raw heatmap data lists
ht_bundle_4T1   <- readRDS("cGAMP_gene_variant_MAST_heatmap_4T1_full.rds")
ht_bundle_EO771 <- readRDS("cGAMP_gene_variant_MAST_heatmap_EO771_full.rds")

# =========================================================================
# 2. Build Scalable Viewport Grobs (Fixed Viewport Dimensions)
# =========================================================================

# Turn off the RStudio warning message globally
ht_opt$message = FALSE

# CRITICAL FIX: Open a temporary null device file to give ComplexHeatmap
# real, non-zero physical layout dimensions to safely compute widths/heights.
pdf(tf <- tempfile(), width = 25, height = 15)

# ---- Process 4T1 Heatmap Grob ----
grid::grid.newpage()
draw(ht_bundle_4T1$heatmap, 
     annotation_legend_list = list(ht_bundle_4T1$legend), 
     heatmap_legend_side = "right", 
     annotation_legend_side = "bottom")
grob_4T1 <- grid::grid.grab()
p_ht_4T1 <- wrap_elements(grob_4T1)

# ---- Process EO771 Heatmap Grob ----
grid::grid.newpage()
draw(ht_bundle_EO771$heatmap, 
     annotation_legend_list = list(ht_bundle_EO771$legend), 
     heatmap_legend_side = "right", 
     annotation_legend_side = "bottom")
grob_EO771 <- grid::grid.grab()
p_ht_EO771 <- wrap_elements(grob_EO771)

# Close the temporary file device and clean it up immediately
dev.off()
if(file.exists(tf)) file.remove(tf)

# =========================================================================
# 3. Assemble and Print Layout to Canvas
# =========================================================================

# LEFT COLUMN: Stacks the 5 original layout rows vertically
p_left_column <- (p_scPC_4T1 | p_scPC_EO771) / 
  free(p_elbow) / 
  free(p_corr_pathway) / 
  free(p_sc_FACS) /
  free(p_umap_4T1 | p_umap_EO771) +
  plot_layout(heights = c(4, 1.7, 2, 2, 2))

# RIGHT COLUMN: Stacks heatmaps using clean vectorized grob wrappers
p_right_column <- p_ht_4T1 / p_ht_EO771 + 
  plot_layout(heights = c(1.3, 1)) 

# Complete compilation grid side-by-side
p <- (p_left_column | p_right_column) + plot_layout(widths = c(0.75, 1.25)) 

# Coordinates span from 0 (bottom/left) to 1 (top/right) of the entire 35x23 PDF.
# Tweak these decimals directly to shift labels anywhere across the canvas.
xa <- 0.0; ya <- 1
xb <- 0.0; yb <- 0.65
xc <- 0.0; yc <- 0.51
xd <- 0.0; yd <- 0.34
xe <- 0.0; ye <- 0.18


xf <- 0.43; yf <- 1
xg <- 0.43; yg <- 0.45


# Construct individual custom floating label elements
labels_grob_layer <- grobTree(
  textGrob("a", x = unit(xa, "npc"), y = unit(ya, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 30)),
  textGrob("b", x = unit(xb, "npc"), y = unit(yb, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 30)),
  textGrob("c", x = unit(xc, "npc"), y = unit(yc, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 30)),
  textGrob("d", x = unit(xd, "npc"), y = unit(yd, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 30)),
  textGrob("e", x = unit(xe, "npc"), y = unit(ye, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 30)),
  textGrob("f", x = unit(xf, "npc"), y = unit(yf, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 30)),
  textGrob("g", x = unit(xg, "npc"), y = unit(yg, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 30)))

p_final_grob <- grobTree(
  patchwork::patchworkGrob(p),
  labels_grob_layer
)

# =========================================================================
# 5. Export Final Canvas Layout Vector Graphic
# =========================================================================
pdf('cGAMP_sc_extended_fig.pdf', width = 35, height = 23)
grid::grid.draw(p_final_grob)
dev.off()

svglite::svglite('cGAMP_sc_extended_fig.svg', width = 35, height = 23)
grid::grid.draw(p_final_grob)
dev.off()
