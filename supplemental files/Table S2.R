###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-25

# Supplemental Table S2

# ---------------------------------------
# Packages
# ---------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)
library(showtext)
library(lme4)
library(lmerTest)
library(performance)
library(lubridate)
library(gridExtra)
library(cowplot)
library(flextable)
library(performance)
library(rcompanion)
library(lme4)
library(lmerTest)
library(nlme)
library(MASS)
library(glmmTMB)

# --------------------------------------
# PLOT THEMES
# --------------------------------------
source("seq diversity - plot themes.R")

# ---------------------------------------
# FUNCTIONS
# ---------------------------------------
source("seq diversity - functions.R")

# ---------------------------------------
# Load data
# ---------------------------------------

# combined data files for each geography
load(file = "data/combined_data.Rdata")

# --------------------------------------
# Table S2  spearman correlation values
# --------------------------------------


# correlation for each genome region predictor for pi and H
row_1 <- cor_ci_function(
  dataframe = dat_state,
  Region = "Genome-wide",
  predictor_name = "genomewide_pi_state_3w",
  predictor_group = "Pi"
)

row_2 <- cor_ci_function(
  dataframe = dat_state,
  Region = "NSPs 5 and 6",
  predictor_name = "orf_pi_state_3w",
  predictor_group = "Pi"
)

row_3 <- cor_ci_function(
  dataframe = dat_state,
  Region = "2' O-Mtase",
  predictor_name = "cov_mt_2_pi_state_3w",
  predictor_group = "Pi"
)

row_4 <- cor_ci_function(
  dataframe = dat_state,
  Region = "Spike Protein",
  predictor_name = "spike_pi_state_3w",
  predictor_group = "Pi"
)

row_5 <- cor_ci_function(
  dataframe = dat_state,
  Region = "S1 NTD",
  predictor_name = "ntd_pi_state_3w",
  predictor_group = "Pi"
)

row_6 <- cor_ci_function(
  dataframe = dat_state,
  Region = "S1 RBD",
  predictor_name = "rbd_pi_state_3w",
  predictor_group = "Pi"
)

# repeat for h values
row_7 <- cor_ci_function(
  dataframe = dat_state,
  Region = "Genome-wide",
  predictor_name = "genomewide_h_state_3w",
  predictor_group = "H"
)

row_8 <- cor_ci_function(
  dataframe = dat_state,
  Region = "NSPs 5 and 6",
  predictor_name = "orf_h_state_3w",
  predictor_group = "H"
)

row_9 <- cor_ci_function(
  dataframe = dat_state,
  Region = "2' O-Mtase",
  predictor_name = "cov_mt_2_h_state_3w",
  predictor_group = "H"
)

row_10 <- cor_ci_function(
  dataframe = dat_state,
  Region = "Spike Protein",
  predictor_name = "spike_h_state_3w",
  predictor_group = "H"
)

row_11 <- cor_ci_function(
  dataframe = dat_state,
  Region = "S1 NTD",
  predictor_name = "ntd_h_state_3w",
  predictor_group = "H"
)

row_12 <- cor_ci_function(
  dataframe = dat_state,
  Region = "S1 RBD",
  predictor_name = "rbd_h_state_3w",
  predictor_group = "H"
)

# correlation for variant count
row_13 <- cor_ci_function(
  dataframe = dat_state,
  Region = "",
  predictor_name = "n_variants_no_thresh_3w_mean",
  predictor_group = "Freyja variant count"
)

# combine into one table
cor_table <- bind_rows(
  row_1, row_2, row_3, row_4, row_5, row_6,
  row_7, row_8, row_9, row_10, row_11, row_12,
  row_13
)

# make it a flextable and export
# save table as flextable with rounded values
cor_ft <- table_as_flex_function(
  dataframe = cor_table,
  title = "Table 2: Spearman correlations for diversity")

# save
save_as_docx(FitFlextableToPage(cor_ft), 
             path = paste("Supplemental files/",
                          "Table S2.docx",
                          sep = "")
)
