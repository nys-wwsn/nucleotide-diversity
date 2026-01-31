###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-26

# Supplemental Figure S7

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

# --------------------------------------------------------------------------
# Figure S7 - selective sweep figure showin the virus takeover with the 5% 
# threshold for lineages
# --------------------------------------------------------------------------

# JN.1 was first detected in statewide wasteater samples October 29, 2023
# JN.1 was no longer dominant after April 14, 2024

dat_state$jn1 <- ifelse(
  dat_state$week >= "2023-10-29" & dat_state$week <= "2024-04-14",
  "JN.1", "Other Variants"
)

# correlation with freyja lineages
freyja_cor_thresh <- 
  ggplot(data = dat_state,
         aes(x = n_lineages_5_3w_mean,
             y = case_incidence)
  )+
  geom_point(color = "grey60", size = 3, alpha = 0.7,
             show.legend = FALSE)+
  geom_point(data = dat_state %>% filter(jn1 == "JN.1"),
             aes(color = jn1), size =3, alpha = 1)+
  stat_cor(method = "spearman"
  )+
  theme_dth_1+
  labs(x = "Mean number of Freyja lineages",
       y = "COVID-19 cases per 100k")+
  scale_color_manual(values = c("JN.1" = "gold"),
                     name = "",
                     label = "JN.1 lineage dominant")+
  theme(legend.position = "bottom",
        legend.background = element_blank())
freyja_cor_thresh

# time series plot
var_plot_thresh <- 
  ggplot()+
  geom_bar(data = dat_state,
           aes(x = week,
               y = case_incidence,
               fill = "grey60"),
           position = "stack",
           stat = "identity")+
  geom_line(data = dat_state,
            aes(x = week, y= n_lineages_5_3w_mean*15,
                color = "darkblue"),
            linewidth = 1
  )+
  theme_dth_1+
  theme(legend.position = "bottom",
        legend.background = element_blank())+
  labs(
    x = "")+
  scale_y_continuous(
    "COVID-19 cases per 100k", 
    sec.axis = sec_axis(~ ./ 15, name = "Lineage count")
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    label = c("Cases"),
                    name = "")+
  scale_color_manual(values = c("darkblue" = "darkblue"),
                     label = "Freyja lineages",
                     name = "")+
  
  geom_vline(xintercept = as.Date("2023-10-29"),
             color = "black",
             linetype = "dashed")+
  geom_vline(xintercept = as.Date("2024-04-14"),
             color = "black",
             linetype = "dashed")+
  annotate("label",
           label = "JN.1\nvariant surge",
           x = as.Date("2024-01-31"),
           y = 100,
           fill = "gold",
           alpha =0.5
  )

var_plot_thresh

# make it a panel
plot_grid(
  var_plot_thresh,
  freyja_cor_thresh,
  nrow = 1,
  rel_widths = c(2,1),
  align = 'h', axis = 'tb',
  labels = c("A", "B")
)

# save
png("Figures/Figure S7.png",
    units = "in",
    width = 13, height = 6,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(
  var_plot_thresh,
  freyja_cor_thresh,
  nrow = 1,
  rel_widths = c(2,1.5),
  align = 'h', axis = 'tb',
  labels = c("A", "B")
)
showtext_end()
dev.off()