###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-25

# Supplemental Figure S4

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


# -----------------------------------------------
# Figure S4 - pi and genome regions time series
# -----------------------------------------------

# ------------------------------------------------
# CORRELATION BETWEEN CLINICAL AND GENOME REGIONS
# ------------------------------------------------

# Statewide time series figure, all genome regions, panel for NYC
dat_state_long <- dat_state %>%
  pivot_longer(cols = c(
    genomewide_pi_state_3w,
    ntd_pi_state_3w,
    cov_mt_2_pi_state_3w,
    orf_pi_state_3w,
    spike_pi_state_3w,
    rbd_pi_state_3w
  ))%>%
  mutate(name_factor = factor(name,
                              levels = c("genomewide_pi_state_3w",
                                         "spike_pi_state_3w",
                                         "orf_pi_state_3w",
                                         "ntd_pi_state_3w",
                                         "rbd_pi_state_3w",
                                         "cov_mt_2_pi_state_3w"),
                              labels = c("Genome",
                                         "Spike",
                                         "ORF NSPs 5 and 6",
                                         "S1 NTD",
                                         "S1 RBD",
                                         "CoV-MT-2"))
  )

pal <- met.brewer(name = "Troy", n = 8)
pal

time_series_1 <- 
  ggplot(data = dat_state_long) +
  geom_bar(aes(x = week, y = case_incidence,
               fill = "grey60"),
           position = "dodge", stat = "identity")+
  geom_line(aes(x = week, y = value * 10000,
                color = name),
            lwd= 1.15)+
  facet_wrap(~name_factor, nrow = 6)+
  theme_dth_1+
  theme(legend.position = "none",
        legend.background = element_rect(color = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0)
  )+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genome-wide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  scale_x_date(date_labels = "%b %y",
               date_breaks = "1 month",
               limits = c(ymd("2023-01-01"),
                          ymd("2025-04-20"))
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    labels = c("Cases\nper 100k"),
                    name = ""
  )+
  labs(x = "",
       title = "Genome region Pi values")+
  scale_y_continuous(
    "COVID-19 cases per 100k", 
    sec.axis = sec_axis(~ ./ 10000, name = expression("Pi"[ww]))
  )

time_series_1

# figure with correlations by genome region
cor_plot <- dat_state %>%
  pivot_longer(cols = c("genomewide_pi_state_3w",
                        "spike_pi_state_3w",
                        "orf_pi_state_3w",
                        "ntd_pi_state_3w",
                        "rbd_pi_state_3w",
                        "cov_mt_2_pi_state_3w")
  ) %>%
  mutate(name_factor = factor(name,
                              levels = c("genomewide_pi_state_3w",
                                         "spike_pi_state_3w",
                                         "orf_pi_state_3w",
                                         "ntd_pi_state_3w",
                                         "rbd_pi_state_3w",
                                         "cov_mt_2_pi_state_3w"),
                              labels = c("Genome",
                                         "Spike",
                                         "ORF NSPs 5 and 6",
                                         "S1 NTD",
                                         "S1 RBD",
                                         "CoV-MT-2"))
  ) %>%
  ggplot()+
  geom_point(aes(x = value, y = case_incidence, color = name))+
  facet_wrap(~name_factor, nrow = 6)+
  stat_cor(aes(x = value, y = case_incidence),
           method = "spearman",
           size = 4)+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genome-wide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  theme_dth_1 +
  theme(legend.position = "none")+
  labs(title = "Spearman correlation:\nPi ~ case incidence",
       x = "Pi",
       y = "")

cor_plot


cowplot::plot_grid(time_series_1, cor_plot, nrow = 1,
                   rel_widths = c(2,1),
                   align = 'h', axis = 'tb')

# combined plot for pi
# save
png("Supplemental files/Figure S4.png",
    units = "in",
    width = 8.5, height = 11,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
cowplot::plot_grid(time_series_1, cor_plot, nrow = 1,
                   rel_widths = c(2,1),
                   align = 'h', axis = 'tb',
                   labels = c("A", "B"))
showtext_end()
dev.off()

