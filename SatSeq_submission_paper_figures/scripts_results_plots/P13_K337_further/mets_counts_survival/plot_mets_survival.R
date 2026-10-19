# Load libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(survival)
library(survminer)
library(readxl)
library(svglite)
library(stringr)



# Read Excel file (update path + sheet name)
df <- read_excel("Met counts Matt.xlsx",
                 sheet = "prometastatic variants")

# Wide → long
df_long <- df %>%
  pivot_longer(
    cols = everything(),
    names_to = "variant",
    values_to = "count"
  ) %>%
  filter(!is.na(count))

# Set variant order
variant_levels <- c("REF", "STOP", "P13X", "K337X")
df_long <- df_long %>%
  mutate(variant = factor(variant, levels = variant_levels))

# Compute medians
medians <- df_long %>%
  group_by(variant) %>%
  summarise(med = median(count, na.rm = TRUE))

# Get y-position for significance bar (a bit above max)
y_max <- max(df_long$count, na.rm = TRUE)
y_sig <- 140

# Plot
p <- ggplot(df_long, aes(x = variant, y = count, color = variant)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.5) +
  
  # Median bars (colored like points)
  geom_crossbar(
    data = medians,
    aes(x = variant, y = med, ymin = med, ymax = med, color = variant),
    width = 0.5,
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c("REF" = "#66a61e",
                                "STOP" = "#a6761d",
                                "P13X" = "#e7298a",
                                "K337X" = "#542788")) +
  # Significance line between REF (1) and STOP (2)
  annotate("segment", x = 1, xend = 2, y = y_sig, yend = y_sig, color = "black") +
  annotate("text", x = 1.5, y = y_sig * 1.01, label = "**", size = 7) +
  theme_cowplot() +
  ylim(NA,250) +
  labs(x = "Variant", y = "Surface lung mets count") +
  theme(legend.position = "none")

# ggsave(plot=p, 'mets_counts.pdf',
#        width=3.5, height=5)

svglite("mets_counts.svg", width=3.5, height=5)
print(p)
dev.off()

# Read Excel file (update path + sheet name)
df <- read_excel("Met counts Matt.xlsx",
                 sheet = "variant pilot")

# Wide → long
df_long <- df[c("REF", "STOP", "Y104X",	"L138X",	"L329X")] %>%
  pivot_longer(
    cols = everything(),
    names_to = "variant",
    values_to = "count"
  ) %>%
  filter(!is.na(count))

# Set variant order
variant_levels <- c("REF", "STOP", "Y104X",	"L138X",	"L329X")
df_long <- df_long %>%
  mutate(variant = factor(variant, levels = variant_levels))

# Compute medians
medians <- df_long %>%
  group_by(variant) %>%
  summarise(med = median(count, na.rm = TRUE))

# Get y-position for significance bar (a bit above max)
y_max <- max(df_long$count, na.rm = TRUE)
y_sig <- 100

# Plot
p <- ggplot(df_long, aes(x = variant, y = count, color = variant)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  
  # Median bars (colored like points)
  geom_crossbar(
    data = medians,
    aes(x = variant, y = med, ymin = med, ymax = med, color = variant),
    width = 0.5,
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c("REF" = "#66a61e",
                                "STOP" = "#a6761d",
                                "Y104X" = "#fc8d62" ,
                                "L138X" = "#8da0cb",
                                "L329X" = "#A9A9A9")) +
  # Significance line between REF (1) and STOP (2)
  annotate("segment", x = 1, xend = 5, y = y_sig, yend = y_sig, color = "black") +
  annotate("text", x = 3, y = y_sig * 1.01, label = "*", size = 7) +
  theme_cowplot() +
  ylim(NA,110) +
  labs(x = "Variant", y = "Surface lung mets count") +
  theme(legend.position = "none")

# ggsave(plot=p, 'mets_counts.pdf',
#        width=3.5, height=5)

svglite("mets_counts_antimets.svg", width=4, height=5)
print(p)
dev.off()

## Do survival curves
survival_curve_variant <- data.frame()

death_df <- read.csv('lung_mets_points.csv')

for (variant in c('P13A','K337R','REF','STOP')){
  
  # Subset for Variant
  df_variant <- death_df[death_df$variant==variant,]
  df_variant$Status <- 1
  
  # Fit survival model
  surv_obj <- Surv(time = df_variant$Days, event = df_variant$Status)
  fit <- survfit(surv_obj ~ 1, data = df_variant)
  fit$variant <- variant
  
  # add to data
  survival_curve_variant <- rbind(survival_curve_variant,
                                  data.frame(time = fit$time,
                                             surv = fit$surv,
                                             variant = fit$variant)
                                  )
}

# plot survival curves
# 1. Ensure each variant starts at time 0 with survival 1
plot_df <- survival_curve_variant %>%
  group_by(variant) %>%
  group_modify(~ bind_rows(tibble(time = 0, surv = 1), .x))

# 2. Create the plot
p <- ggplot(plot_df, aes(x = time, y = surv, color = variant)) +
  geom_step(linewidth = 1) +
  scale_color_manual(values = c("REF" = "#66a61e",
                                "STOP" = "#a6761d",
                                "P13A" = "#e7298a",
                                "K337R" = "#542788")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0.02)) + 
  xlim(NA,50) +
  theme_cowplot() +
  labs(x = "Days", 
       y = "Survival probability",
       color = "Variant")

svglite("survival.svg", width=6, height=3)
print(p)
dev.off()

# Create P13X mets plot for different AAs
# 1. Load and Clean Data
df <- read.csv("Mammary_Fat_Pad_Primary_vs_Met_boot_4T1.csv")

plot_data <- df %>%
  # Filter for P13X variants
  filter(grepl("^P13[A-Z]$", variant)) %>%
  # Remove "P13" to keep only the substituted amino acid (e.g., "A")
  mutate(amino_acid = str_replace(variant, "P13", ""))

plot_data <- plot_data[plot_data$variant!='P13X',]

# 2. Calculate summary statistics for ordering
df_summary <- plot_data %>%
  group_by(amino_acid) %>%
  summarize(
    mean_fit = mean(observed_fitness, na.rm = TRUE),
    se_fit = sd(observed_fitness, na.rm = TRUE) / sqrt(n())
  ) %>%
  arrange(mean_fit)

# Reorder factor levels so bars increase to the right
plot_data$amino_acid <- factor(plot_data$amino_acid, levels = df_summary$amino_acid)
df_summary$amino_acid <- factor(df_summary$amino_acid, levels = df_summary$amino_acid)

# 3. Create the Plot
p <- ggplot() +
  # Bar layer: color = NA removes the borders
  geom_bar(data = df_summary, 
           aes(x = amino_acid, y = mean_fit, fill = as.numeric(amino_acid)), 
           stat = "identity", color = NA, show.legend = FALSE) +
  # Error bars
  geom_errorbar(data = df_summary, 
                aes(x = amino_acid, ymin = mean_fit - se_fit, ymax = mean_fit + se_fit), 
                width = 0.2, color = "black") +
  # Individual data points
  geom_jitter(data = plot_data, 
              aes(x = amino_acid, y = observed_fitness), 
              color = "black", alpha = 0.5, width = 0.15) +
  # Gradient from light to dark #e7298a
  scale_fill_gradient(low = "#fce7f3", high = "#e7298a") +
  theme_cowplot() +
  geom_hline(yintercept = 0) +
  labs(x = "P13 Substitution", 
       y = "Metastatic Fitness")

svglite("P13_variant_fitness.svg", width=6, height=3)
print(p)
dev.off()
