###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-26

# Supplemental Table S9

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
library(lmtest)

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

# Table S9 descriptives for conc data
# summarize by pcr lab
conc_df_summary <- dat_sewershed %>%
  ungroup() %>%
  filter(!is.na(pcr_lab)) %>%
  group_by(pcr_lab) %>%
  summarize(`Number of sampling sites` = length(unique(facility_id)),
            `Number of samples` = sum(samples, na.rm  = TRUE),
            `First collection week` = min(week, na.rm = TRUE),
            `Last collection week` = max(week, na.rm= TRUE),
            `Mean concentration` = signif(mean(mean_sars2_conc, na.rm = TRUE),
                                          2),
            `Min concentration` = signif(min(mean_sars2_conc, na.rm = TRUE),
                                         2),
            `Median concentration` = signif(median(mean_sars2_conc, 
                                                   na.rm = TRUE),
                                            2),
            `Max concentration` = signif(max(mean_sars2_conc, na.rm = TRUE), 2),
            `SD concentration` = signif(sd(mean_sars2_conc, na.rm = TRUE),
                                        2)
            
  ) %>%
  ungroup() %>%
  mutate(pcr_lab = case_when(
    pcr_lab == "NYC" ~ "NYC DOH",
    pcr_lab == "genesee_orleans_health" ~ "GO Health",
    pcr_lab == "quadrant" ~ "Quadrant",
    pcr_lab == "suny_buffalo" ~ "SUNY Buffalo",
    pcr_lab == "suny_stony_brook" ~ "SUNY Stony Brook",
    pcr_lab == "wadsworth" ~ "Wadsworth Center"
  )) %>%
  rename(`PCR Lab` = pcr_lab) %>%
  mutate(`Min concentration` = ifelse(
    `Min concentration` == 1, 0, `Min concentration`
  ))

# make it flextable
conc_ft <- table_as_flex_function(
  dataframe = conc_df_summary,
  title = "Table S9: Descriptive statistics for concentration data")

# save
save_as_docx(FitFlextableToPage(conc_ft), 
             path = paste0("Supplemental files",
                           "/Table S9- concentration descriptive statistics.docx")
)
