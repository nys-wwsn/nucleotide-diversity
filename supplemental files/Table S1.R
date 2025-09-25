###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-25

# Supplemental Table 1

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

# max and min weeks for the study
min_week <- "2023-01-01"
max_week <- "2025-04-20"

# mean depth per sample - data prior to July 2024
depth_1 <- read.csv("data/other data/depth_stats_20250424.txt") %>%
  tidyr::separate(
    sample, into = c("date", "cdc_id"), sep = 8
  ) %>%
  mutate(date = lubridate::ymd(date),
         mean_depth = mean) %>%
  select(date, cdc_id, mean_depth)

# mean depth post july 2024
depth_2 <- read.csv("data/other data/genomwide_depth_2025-07-02.csv") %>%
  mutate(date = ymd(date)) %>%
  select(date, cdc_id, mean_depth) %>%
  tidyr::separate(
    cdc_id, into = c("cdc_id", "drop"), sep = 12
  )

# bind together, then merge to the dat_sewershed object
depth <- bind_rows(depth_1, depth_2) %>%
  group_by(cdc_id, week = floor_date(date, unit = "weeks", week_start = 7)) %>%
  summarize(mean_depth = mean(mean_depth, na.rm =  TRUE)
  ) %>%
  ungroup() %>%
  rename(facility_id = cdc_id)

dat_sewershed <- left_join(dat_sewershed, depth, by = c("facility_id", "week"))

# --------------------------------------
# Table S1 - descriptive stats for diveristy data
# --------------------------------------

# merge lab data
# summarize by pcr lab
diversity_df_summary <- dat_sewershed %>%
  group_by(pcr_lab) %>%
  filter(!is.na(genomewide_pi)) %>%
  summarize(`Number of sampling sites` = length(unique(facility_id)),
            `Number of samples` = sum(samples, na.rm = TRUE),
            `First collection week` = min(week, na.rm = TRUE),
            `Last collection week` = max(week, na.rm= TRUE),
            `Mean Pi` = signif(mean(genomewide_pi, na.rm = TRUE),2),
            `Min Pi` = signif(min(genomewide_pi, na.rm = TRUE), 2),
            `Median Pi` = signif(median(genomewide_pi, na.rm = TRUE), 2),
            `Max Pi` = signif(max(genomewide_pi, na.rm = TRUE), 2),
            `SD Pi` = signif(sd(genomewide_pi, na.rm = TRUE), 2)
            
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
  rename(`PCR Lab` = pcr_lab) 


# save table as flextable with rounded values
pi_ft <- table_as_flex_function(
  dataframe = diversity_df_summary,
  title = "Table: Descriptive statistics for genome-wide Pi")

# save
save_as_docx(FitFlextableToPage(pi_ft), 
             path = paste("Supplemental files/",
                          "Table S1.docx",
                          sep = "")
)
