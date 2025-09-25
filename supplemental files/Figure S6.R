###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-25

# Supplemental Figures S6

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

# ----------------------------------------------------
# Figure S6 - lag and lead for all genome regions
# ----------------------------------------------------

# correllation with hospital admissions

pi_lag <- lag_lead_function_hosp(dataframe = dat_state,
                                 columns = c("genomewide_pi_state_3w",
                                             "spike_pi_state_3w",
                                             "orf_pi_state_3w",
                                             "ntd_pi_state_3w",
                                             "rbd_pi_state_3w",
                                             "cov_mt_2_pi_state_3w"),
                                 lags = seq(1:10))


# plot
lag_plot <- 
  ggplot(data = pi_lag,
         aes(x = lag, y = spearman_cor, color = name))+
  geom_point(shape = 1, size = 2)+
  geom_point(shape = 17,
             data = pi_lag %>%
               group_by(name) %>%
               filter(spearman_cor == max(spearman_cor)),
             size = 3
  )+
  geom_line()+
  ylim(0,1)+
  theme_dth_1+
  theme(legend.position = "bottom")+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genome-wide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = seq(-5, 5, by = 1),
                     limits = c(-5,5))+
  labs(x = "Lead time (weeks)",
       y = "Spearman correlation\ncoefficient",
       title= "Hospital incidence ~ Pi (wastewater) for\nstatewide data")+
  annotate("text",
           x = -3, y = 0.09,
           label = "signal lags incidence")+
  annotate("text",
           x = 3, y = 0.09,
           label = "signal leads incidence")

# repeat for shannon

lag_h <- lag_lead_function_hosp(dataframe = dat_state,
                                columns = c("genomewide_h_state_3w",
                                            "spike_h_state_3w",
                                            "orf_h_state_3w",
                                            "ntd_h_state_3w",
                                            "rbd_h_state_3w",
                                            "cov_mt_2_h_state_3w"),
                                lags = seq(1:10))

lag_plot_2 <- ggplot(data = lag_h ,
                     aes(x = lag, y = pearson_cor, color = name))+
  geom_point(shape = 1, size = 2)+
  geom_point(shape = 17,
             data = lag_h %>%
               group_by(name) %>%
               filter(pearson_cor == max(pearson_cor)),
             size = 3
  )+
  geom_line()+
  ylim(0,1)+
  theme_dth_1+
  theme(legend.position = "bottom")+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genome-wide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = seq(-5, 5, by = 1),
                     limits = c(-5,5))+
  labs(x = "Lead time (weeks)",
       y = "Spearman correlation\ncoefficient",
       title= "Hospital incidence ~ Shannon H (wastewater) for\nstatewide data")+
  annotate("text",
           x = -3, y = 0.09,
           label = "signal lags incidence")+
  annotate("text",
           x = 3, y = 0.09,
           label = "signal leads incidence")

# add  cases
pi_lag_2 <- lag_lead_function(dataframe = dat_state,
                              columns = c("genomewide_pi_state_3w",
                                          "spike_pi_state_3w",
                                          "orf_pi_state_3w",
                                          "ntd_pi_state_3w",
                                          "rbd_pi_state_3w",
                                          "cov_mt_2_pi_state_3w"),
                              lags = seq(1:10))


# plot
lag_plot_case <- 
  ggplot(data = pi_lag_2,
         aes(x = lag, y = spearman_cor, color = name))+
  geom_point(shape = 1, size = 2)+
  geom_point(shape = 17,
             data = pi_lag_2 %>%
               group_by(name) %>%
               filter(spearman_cor == max(spearman_cor)),
             size = 3
  )+
  geom_line()+
  ylim(0,1)+
  theme_dth_1+
  theme(legend.position = "bottom")+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genome-wide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = seq(-5, 5, by = 1),
                     limits = c(-5,5))+
  labs(x = "Lead time (weeks)",
       y = "Spearman correlation\ncoefficient",
       title= "Case incidence ~ Pi (wastewater) for\nstatewide data")+
  annotate("text",
           x = -3, y = 0.09,
           label = "signal lags incidence")+
  annotate("text",
           x = 3, y = 0.09,
           label = "signal leads incidence")


# repeat for shannon

lag_h_2 <- lag_lead_function(dataframe = dat_state,
                             columns = c("genomewide_h_state_3w",
                                         "spike_h_state_3w",
                                         "orf_h_state_3w",
                                         "ntd_h_state_3w",
                                         "rbd_h_state_3w",
                                         "cov_mt_2_h_state_3w"),
                             lags = seq(1:10))

lag_plot_2_case <- ggplot(data = lag_h_2 ,
                          aes(x = lag, y = pearson_cor, color = name))+
  geom_point(shape = 1, size = 2)+
  geom_point(shape = 17,
             data = lag_h_2 %>%
               group_by(name) %>%
               filter(pearson_cor == max(pearson_cor)),
             size = 3
  )+
  geom_line()+
  ylim(0,1)+
  theme_dth_1+
  theme(legend.position = "bottom")+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genome-wide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = seq(-5, 5, by = 1),
                     limits = c(-5,5))+
  labs(x = "Lead time (weeks)",
       y = "Spearman correlation\ncoefficient",
       title= "Case incidence ~ Shannon H (wastewater) for\nstatewide data")+
  annotate("text",
           x = -3, y = 0.09,
           label = "signal lags incidence")+
  annotate("text",
           x = 3, y = 0.09,
           label = "signal leads incidence")

# Add freyja variant counts and concentration data

# case
lag_freyja_case <- lag_lead_function(dataframe = dat_state,
                                     columns = c("n_variants_no_thresh_3w_mean",
                                                 "mean_sars2_conc_state_3w"),
                                     lags = seq(1:10))

lag_freyja_case <- ggplot(data = lag_freyja_case ,
                          aes(x = lag, y = pearson_cor, color = name))+
  geom_point(shape = 1, size = 2)+
  geom_point(shape = 17,
             data = lag_h_2 %>%
               group_by(name) %>%
               filter(pearson_cor == max(pearson_cor)),
             size = 3
  )+
  geom_line()+
  ylim(0,1)+
  theme_dth_1+
  theme(legend.position = "bottom")+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genome-wide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = seq(-5, 5, by = 1),
                     limits = c(-5,5))+
  labs(x = "Lead time (weeks)",
       y = "Spearman correlation\ncoefficient",
       title= "Case incidence ~ Shannon H (wastewater) for\nstatewide data")+
  annotate("text",
           x = -3, y = 0.09,
           label = "signal lags incidence")+
  annotate("text",
           x = 3, y = 0.09,
           label = "signal leads incidence")
# hosp

# extract legend function
mylegend <- g_legend(lag_plot)



# save
png("Supplemental files/Figure S6.png",
    units = "in",
    width = 10, height = 11,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
ggarrange(ggarrange(lag_plot_case+ theme(legend.position = "none"),
                    lag_plot_2_case + theme(legend.position = "none"),
                    nrow = 1,
                    labels = c("A", "B")),
          ggarrange(lag_plot+ theme(legend.position = "none"),
                    lag_plot_2 + theme(legend.position = "none"),
                    nrow = 1,
                    labels = c("C", "D")),
          mylegend, 
          nrow=3,
          heights=c(5,5, 1))
dev.off()
