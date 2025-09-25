###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-25

# Main analysis - Figure 4

# --------------------------------------
# PACKAGES
# --------------------------------------

library(docstring)
library(showtext)
library(ggplot2)
library(sf)
library(zoo)
library(stringr)
library(lubridate)
library(flextable)
library(tidyr)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(confintr)
library(MetBrewer)
library(broom)
library(cowplot)
library(performance)
library(rcompanion)

# Load functions
source("seq diversity - functions.R")

# Load plot themes
source("seq diversity - plot themes.R")

# -----------------------------------------------------------------------------

# --------------------------------------
# DATA LOAD
# --------------------------------------

# combined data files for each geography
load(file = "data/combined_data.Rdata")

# max and min weeks for the study
min_week <- "2023-01-01"
max_week <- "2025-04-20"

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Figure 4 - correlation scatterplots for all 3 measures plus concentration
# -----------------------------------------------------------------------------

# Fisher test for different correlations
# compare correlation coefficients to see if they are significantly different

# resources
# https://rpubs.com/JimGrange/comparing_correlations
# https://www.statisticssolutions.com/comparing-correlation-coefficients/#:~:text=Researchers%20recommend%20this%20comparison%20when,transformation%20to%20compare%20z%20scores.
# https://search.r-project.org/CRAN/refmans/DescTools/html/FisherZ.html

# compute Rho for each correlation
cor_ntd <- cor(dat_state$case_incidence, dat_state$ntd_pi_state_3w,
               use = "na.or.complete",
               method = "spearman")

cor_conc <- cor(dat_state$case_incidence, dat_state$mean_sars2_conc_state_3w,
                use = "na.or.complete",
                method = "spearman")

# package
library(DescTools)

# Fisher z score for the correlations
z1 <- FisherZ(cor_ntd)
z2 <- FisherZ(cor_conc)

# number of observations
n1 <- 119
n2 <- 119

# compare z scores to get Fisher's Z
(z1 - z2)/sqrt((1/(n1-3)) + (1/(n2-3)))

# try replacing R with Rho


# ntd plot
pi_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = ntd_pi_state_3w,
             y = case_incidence))+
  geom_point(color = pal[6])+
  stat_cor(method = "spearman",
           cor.coef.name = "rho")+
  theme_dth_1+
  labs(
    title = "S1 NTD region",
    x = expression(paste( pi[ww])),
    y = "Case incidence"
  )

# shannon diversity correlation plot
h_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = ntd_h_state_3w,
             y = case_incidence))+
  geom_point(color = pal[5])+
  stat_cor(method = "spearman",
           cor.coef.name = "rho")+
  theme_dth_1+
  labs(
    title = "S1 NTD region",
    x = expression(paste( "H"[ww], sep = "")),
    y = "Case incidence"
  )

# freyja variant count
var_count_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = n_variants_no_thresh_3w_mean,
             y = case_incidence))+
  geom_point(color = pal[1])+
  stat_cor(method = "spearman",
           cor.coef.name = "rho")+
  theme_dth_1+
  labs(
    title = "Freyja variant count",
    x = "Variant count",
    y = "Case incidence"
  )

# concentration
conc_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = mean_sars2_conc_state_3w,
             y = case_incidence))+
  geom_point(color = pal[2])+
  stat_cor(method = "spearman",
           cor.coef.name = "rho")+
  theme_dth_1+
  labs(
    title = "Concentration",
    x = expression(paste("Conc, copies/", "m","L", sep = "")),
    y = "Case incidence"
  )

# title
title <- ggdraw() + draw_label("Correlation with case incidence per 100,000", 
                               fontface='bold')

# panel
panel <-
  plot_grid(
    pi_cor_plot,
    h_cor_plot,
    var_count_cor_plot,
    conc_cor_plot,
    nrow = 2,
    align = "v",
    axis = "l",
    labels = c("A", "B", "C", "D")
  )

# add title to panel
case_panel <- plot_grid(title, panel, ncol=1, rel_heights=c(0.1, 1))

# add hospitalization correlations as panels E through H

# ntd plot
pi_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = ntd_pi_state_3w,
             y = hosp_incidence))+
  geom_point(color = pal[6])+
  stat_cor(method = "spearman",
           cor.coef.name = "rho")+
  theme_dth_1+
  labs(
    title = "S1 NTD region",
    x = expression(paste( pi[ww])),
    y = "Hospitalization incidence"
  )

# shannon diversity correlation plot
h_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = ntd_h_state_3w,
             y = hosp_incidence))+
  geom_point(color = pal[5])+
  stat_cor(method = "spearman",
           cor.coef.name = "rho")+
  theme_dth_1+
  labs(
    title = "S1 NTD region",
    x = expression(paste( "H"[ww], sep = "")),
    y = "Hospitalization incidence"
  )

# freyja variant count
var_count_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = n_variants_no_thresh_3w_mean,
             y = hosp_incidence))+
  geom_point(color = pal[1])+
  stat_cor(method = "spearman",
           cor.coef.name = "rho")+
  theme_dth_1+
  labs(
    title = "Freyja variant count",
    x = "Variant count",
    y = "Hospitalization incidence"
  )

# concentration
conc_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = mean_sars2_conc_state_3w,
             y = hosp_incidence))+
  geom_point(color = pal[2])+
  stat_cor(method = "spearman",
           cor.coef.name = "rho")+
  theme_dth_1+
  labs(
    title = "Concentration",
    x = expression(paste("Conc, copies/", "m","L", sep = "")),
    y = "Hospitalization incidence"
  )

# title
title <- ggdraw() + draw_label(
  "Correlation with hospitalization incidence per 100,000", 
  fontface='bold')

# panel
panel <-
  plot_grid(
    pi_cor_plot,
    h_cor_plot,
    var_count_cor_plot,
    conc_cor_plot,
    nrow = 2,
    align = "v",
    axis = "l",
    labels = c("E", "F", "G", "H")
  )

# add title to panel
hosp_panel <- plot_grid(title, panel, ncol=1, rel_heights=c(0.1, 1))

# combine panels
plot_grid(
  case_panel,
  hosp_panel,
  nrow = 2
)

# save
png("Figures/Figure 4.png",
    units = "in",
    width = 8.5, height = 11,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(
  case_panel,
  hosp_panel,
  nrow = 2
)
showtext_end()
dev.off()
