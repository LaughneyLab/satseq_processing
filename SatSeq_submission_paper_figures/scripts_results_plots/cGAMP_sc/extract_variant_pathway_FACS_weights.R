rm(list=ls())#Clear environment

#################
## libraries ####
#################
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(stringr)
library(tidyverse)

setwd(getSrcDirectory(function(){})[1]) # set working directory to where this script is

## Plot FDR and NES match across cell lines and FACS-----
# Load mets fitness weights
# STING sequence
wt_seq = paste0("MPYSNLHPAIPRPRGHRSKYVALIFLVASLMILWVAKDPPNHTLKYLALHLASHELGLLL",
                "KNLCCLAEELCHVQSRYQGSYWKAVRACLGCPIHCMAMILLSSYFYFLQNTADIYLSWMF",
                "GLLVLYKSLSMLLGLQSLTPAEVSAVCEEKKLNVAHGLAWSYYIGYLRLILPGLQARIRM",
                "FNQLHNNMLSGAGSRRLYILFPLDCGVPDNLSVVDPNIRFRDMLPQQNIDRAGIKNRVYS",
                "NSVYEILENGQPAGVCILEYATPLQTLFAMSQDAKAGFSREDRLEQAKLFCRTLEEILED",
                "VPESRNNCRLIVYQEPTDGNSFSLSQEVLRHIRQEEKEEVTMNAPMTSVAPPPSVLSQEP",
                "RLLISGMDQPLPLRTDLI")

# split the string into a vector so we can index it by position
seq_vector <- strsplit(wt_seq, "")[[1]]

# Read inferred additive weights for variants
weights_df <- data.frame()
readouts <- c('ATF3', 'Dual_IFN', 'NFkB', 'STING', 'Proliferation')

for (readout in readouts){
  
  # Read weight intercept and average over experiments
  weight_intercept <- mean(read.csv(paste0("../FACS_position_subs_modeling/inferred_intercept_experiment_FACS_4T1_", 
                                           readout, ".csv"))[, 2], 
                           na.rm = TRUE)
  # Read weights
  weights <- read.csv(paste0('../FACS_position_subs_modeling/inferred_weights_FACS_4T1_',readout,'.csv'))
  
  colnames(weights) <- c('position', 'additive','interaction')
  
  # Add intercept
  weights$additive <- weights$additive
  
  weights <- weights %>%
    pivot_longer(!position, names_to = 'effect', values_to = 'weight')
  
  weights$position <- weights$position +1
  
  weights$readout <- readout
  
  # Add annotation for top and bottom 10
  weights <- weights %>%
    mutate(
      label = paste0(seq_vector[position], position))
  
  
  weights_df <- rbind(weights_df,weights)
  
}

weights_df <- weights_df %>%
  filter(effect == "additive") %>%
  # Combine position and label to create the 'variant' ID (e.g., A1)
  mutate(variant = label) %>%
  # Select only the necessary columns for the final shape
  select(variant, readout, weight) %>%
  # Pivot the readout names into separate columns
  pivot_wider(
    names_from = readout, 
    values_from = weight,
    names_glue = "{readout}_expression"
  )

# For STOP, RRDD and dCTT, read and average over original data
# Helper to define sheet names and regex patterns
get_config <- function(readout) {
  if (readout == 'Proliferation') {
    list(up = '4T1_Initial_vs_MitogenHigh_p20', dn = NULL, 
         p_up = "4T1_p20_(\\d+)_", p_dn = NULL)
  } else if (readout == 'ATF3') {
    list(up = '4T1_ATF3_FACS_All_vs_Up', dn = '4T1_ATF3_FACS_All_vs_Down', 
         p_up = "ATF3_R(\\d+)_Up", p_dn = "ATF3_R(\\d+)_Dn")
  } else if (readout == 'NFkB') {
    list(up = '4T1_NFkB_FACS_All_vs_Up', dn = '4T1_NFkB_FACS_All_vs_Down', 
         p_up = "NFkB_UP(\\d+)", p_dn = "NFkB_DN(\\d+)")
  } else if (readout == 'Dual_IFN') {
    list(up = '4T1_Dual_IFN_FACS_All_vs_Up', dn = '4T1_Dual_IFN_FACS_All_vs_Down', 
         p_up = "R(\\d+)_Up", p_dn = "R(\\d+)_Down")
  } else { # STING
    list(up = '4T1_STING_FACS_All_vs_Up', dn = '4T1_STING_FACS_All_vs_Down', 
         p_up = "Up_(\\d+)", p_dn = "Down_(\\d+)")
  }
}

file_path <- '../../datasets/FACS/FACS_ATF3_IFN_NFkB_STING_Proliferation_4T1_EO771_all.xlsx'

results_list <- lapply(readouts, function(name) {
  cfg <- get_config(name)
  df_up <- read_excel(file_path, sheet = cfg$up)
  
  if (name == "Proliferation") {
    # For Proliferation, use observed_fitness directly from df_up
    df_processed <- df_up %>%
      rename(expression = observed_fitness)
  } else {
    # For others, merge and subtract
    df_down <- read_excel(file_path, sheet = cfg$dn)
    
    # Extract reps (default to 1 if NA)
    df_up$rep <- as.integer(str_extract(df_up$experiment, cfg$p_up, group = 1))
    df_down$rep <- as.integer(str_extract(df_down$experiment, cfg$p_dn, group = 1))
    df_up$rep[is.na(df_up$rep)] <- 1
    df_down$rep[is.na(df_down$rep)] <- 1
    
    df_processed <- inner_join(df_up, df_down, by = c("variant", "rep", "is_control"), suffix = c("_Up", "_Down")) %>%
      mutate(expression = observed_fitness_Up - observed_fitness_Down)
  }
  
  # Common cleaning, scaling, and WT subtraction
  df_final <- df_processed %>%
    mutate(variant = str_replace_all(variant, c("REF"= "X379X", "STOP" = "X380X",
                                                "ΔCTT" = "X381X", "RRDD" = "X382X"))) %>%
    extract(variant, into = c('wt', 'pos', 'mut'), regex = '^([a-zA-Z]+)([0-9.]+)([a-zA-Z]+)$', remove = FALSE) %>%
    mutate(variant = str_c(wt, pos)) %>%
    group_by(variant) %>%
    summarize(expression = mean(expression, na.rm = TRUE)) %>%
    #mutate(expression = as.numeric(scale(expression))) %>%
    drop_na() %>%
    mutate(variant = str_replace_all(variant, c("X380"= "STOP",
                                                "X381" = "dCTT",
                                                "X382" = "RRDD")))
  
  # Subtract WT (X379) value
  wt_val <- df_final$expression[df_final$variant == 'X379']
  df_final$expression <- df_final$expression - wt_val
  
  # Rename column to {readout}_expression
  colnames(df_final)[colnames(df_final) == "expression"] <- paste0(name, "_expression")
  return(df_final)
})

# Final merge
weights_df_original <- reduce(results_list, full_join, by = "variant")

# Add STOP and dCTT to weights df
weights_df <- rbind(weights_df,
                    weights_df_original[
                      weights_df_original$variant %in% c("STOP","dCTT","RRDD"),]) %>%
  rename(prolif_expression = Proliferation_expression) %>%
  # Scale FACS data to SDs
  mutate(across(ends_with("_expression"), ~ . / sd(., na.rm = TRUE)))

# Save data
write.csv(weights_df,'FACS_weights_all.csv')
