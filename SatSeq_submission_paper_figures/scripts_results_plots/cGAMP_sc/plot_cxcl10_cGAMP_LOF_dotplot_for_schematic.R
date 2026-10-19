rm(list=ls())#Clear environment

library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(svglite)
library(cowplot)
library(ggrepel)
library(scales)


setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

# empty df
ranked_variants = data.frame()

cell_line <- "4T1"
  
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

# Take only cxcl10 data and sort variants
# 2. Extract Cxcl10 data and prepare for plotting
# Assuming all_data_combined is: rownames = genes, colnames = variants
cxcl10_data <- all_data_combined["Cxcl10", ]
plot_df <- data.frame(
  Variant = colnames(all_data_combined),
  Expression = as.numeric(cxcl10_data)
)

# 3. Calculate mean and standard deviation (SD)
mean_val <- mean(plot_df$Expression)
sd_val <- sd(plot_df$Expression)

# 4. Define gain (>2SD) and loss (< -2SD) of function
plot_df <- plot_df %>%
  mutate(
    Type = case_when(
      Expression > ( 1 * sd_val) ~ "Gain of function",
      Expression < (- 1 * sd_val) ~ "Loss of function",
      TRUE ~ "Unchanged"
    ),
    # Only label gain/loss for clarity
    Label = ifelse(Type != "Unchanged", Variant, "")
  ) %>%
 arrange(Expression) %>%
  mutate(
    Rank = row_number(),
    Label = ifelse(Type != "Unchanged", Variant, "")
  ) %>%
  mutate(
    Position = as.numeric(gsub("\\D", "", Variant)))

# 2. Create the ranked dot plot
p <- ggplot(plot_df, aes(x = Rank, y = Expression)) +
  # Main points
  geom_point(aes(color = Type), size = 1) +
  # Customize colors
  scale_color_manual(values = c("Gain of function" = "#E41A1C", 
                                "Loss of function" = "#377EB8", 
                                "Unchanged" = "black")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7393B3") + 
  # Add labels for the outliers
  geom_text_repel(aes(label = Label), 
                  size = 2.5,
                  box.padding = 0.35,
                  min.segment.length = 0.08,
                  max.overlaps = 20) +
  theme_cowplot() +
  # Apply pseudo_log here (this is what brings the point back for geom_text_repel)
  scale_y_continuous(trans = pseudo_log_trans(sigma = 10, base = 10),
                     breaks = scales::breaks_pretty(n = 4)) +
  labs(x = "Variant rank",
       y = "Differential expression",
       color = NULL)

ggsave(plot=p, "cGAMP_4T1_Cxcl10_variant_rankplot.svg", width=6, height=3)

# Create position ordered plot
p <- ggplot(plot_df, aes(x = Position, y = Expression)) +
  # Main points
  geom_point(aes(color = Type), size = 1) +
  geom_line() +
  # Customize colors
  scale_color_manual(values = c("Gain of function" = "#E41A1C", 
                                "Loss of function" = "#377EB8", 
                                "Unchanged" = "black")) +
  # Add labels for the outliers
  geom_text_repel(aes(label = Label), 
                  size = 3,
                  box.padding = 0.35,
                  min.segment.length = 0.08,
                  max.overlaps = 20) +
  theme_cowplot() +
  ylim(-50,0.1) +
  # # Apply pseudo_log here (this is what brings the point back for geom_text_repel)
  # scale_y_continuous(trans = pseudo_log_trans(sigma = 10, base = 10),
  #                    breaks = scales::breaks_pretty(n = 4)) +
  labs(x = "Mutated residue position",
       y = "Differential expression",
       color = NULL)

ggsave(plot=p, "cGAMP_4T1_Cxcl10_variant_pos_ordered.svg", width=9, height=2.5)

