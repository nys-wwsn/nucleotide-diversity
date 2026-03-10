###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-26

# Supplemental Table S5

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


# --------------------------------------------------------------------------
# Table S5 - region model
# --------------------------------------------------------------------------


# REGIONAL MODEL - TIME SERIES CORRECTION
#  samples by region
samples_region <- dat_sewershed %>%
  group_by(region, week) %>%
  summarize(seq_samples = sum(samples, na.rm = TRUE),
            conc_samples = sum(conc_samples, na.rm = TRUE)) %>%
  ungroup()

dat_region <- left_join(dat_region, samples_region,
                        by = c("region", "week"))

# make week a factor
dat_region$date_Factor <- factor(dat_region$week,
                                 levels = seq.Date(
                                   from = min(as_date(dat_region$week)), 
                                   to = max(as_date(dat_region$week)),
                                   by = 1)
)
dat_region$date_Factor <- cut(dat_region$week, breaks = "week")

# glmmtmb with negative binomial and time series correction
model_region_nb_time <- glmmTMB(case_incidence ~ 
                                  + scale(ntd_pi_region_3w)
                                + scale(log(mean_sars2_conc_region_3w))
                                + scale(depth_region_3w)
                                + scale(seq_samples)
                                + scale(conc_samples)
                                + scale(mean_coverage)
                                +ar1(date_Factor + 0|region),
                                family = nbinom2()
                                ,
                                data = dat_region,
                                control=glmmTMBControl(optimizer=optim, 
                                                       optArgs=list(method="BFGS"))
)

summary(model_region_nb_time)
AIC(model_region_nb_time)

region_table <- model_summary_function(
  dataframe = dat_region,
  outcome = "case_incidence",
  glmer_option = "glmmtmb",
  model = model_region_nb_time
) %>%
  mutate(group = "S1 NTD Pi") %>%
  rename(p_value = `Pr(>|z|)`)

#  save table
# edit estimate column
region_table$variable<- c(
  "Intercept",
  "S1 NTD Pi",
  "SARS-CoV-2 concentration",
  "Mean sample depth of read",
  "Number of samples sequenced",
  "Number of samples collected",
  "Mean genome coverage",
  "n",
  "Efron R2",
  "AIC"
) 

# round p value
region_table$p_value <- round_2(as.numeric(region_table$p_value))

# make it an ft table
t <- table_as_flex_function(dataframe = region_table,
                            title = "Tale: Region generalized liner model results for S1 NTD")

t <- set_header_labels(t,
                       values = list(
                         group = "Model",
                         variable = "Variable/Metric",
                         Est = "Estimate",
                         `Std. Error` = "Standard Error (SE)",
                         `t value2` = "t value",
                         `Pr(>|t|)` = "P value"
                       ))

t

# save
save_as_docx(FitFlextableToPage(t), 
             path = paste("Figures/",
                          "Table S5.docx",
                          sep = "")
)