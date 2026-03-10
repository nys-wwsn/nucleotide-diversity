###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-25

# Supplemental Table S4

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
# Table S4 - county model
# --------------------------------------------------------------------------

# COUNTY MODEL - TIME SERIES CORRECTION
# number of samples by county
samples_county <- dat_sewershed %>%
  group_by(county, week) %>%
  summarize(seq_samples = sum(samples, na.rm = TRUE),
            conc_samples = sum(conc_samples, na.rm = TRUE)) %>%
  ungroup()

dat_county <- left_join(dat_county, samples_county,
                        by = c("county", "week"))

# glmmtbm
# ar1 adjustment
dat_county$date_Factor <- factor(dat_county$week,
                                 levels = seq.Date(
                                   from = min(as_date(dat_county$week)), 
                                   to = max(as_date(dat_county$week)),
                                   by = 1)
)
dat_county$date_Factor <- cut(dat_county$week, breaks = "week")


# gaussian with random effect only
# negative binomial time series glmm
model_county_nb_time <- glmmTMB(case_incidence ~ 
                                  + scale(ntd_pi_county_3w)
                                + scale(log(mean_sars2_conc_county_3w))
                                + scale(depth_county_3w)
                                + scale(seq_samples)
                                + scale(conc_samples)
                                + scale(mean_coverage)
                                +ar1(date_Factor + 0|county),
                                family = nbinom2()
                                ,
                                data = dat_county,
                                control=glmmTMBControl(optimizer=optim, 
                                                       optArgs=list(method="BFGS"))
)

summary(model_county_nb_time)
AIC(model_county_nb_time)

# table

county_table <-  model_summary_function(
  dataframe = dat_county,
  outcome = "case_incidence",
  model = model_county_nb_time,
  glmer_option = "glmmtmb"
)%>%
  mutate(group = "S1 NTD model") %>%
  rename(p_value = `Pr(>|z|)`)

#
#  save table
# edit estimate column
county_table$variable<- c(
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
county_table$p_value <- round_2(as.numeric(county_table$p_value))

# make it an ft table
t <- table_as_flex_function(dataframe = county_table,
                            title = "Tale: County generalized liner model results for S1 NTD")

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
                          "Table S4.docx",
                          sep = "")
)
