###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2026-01-16

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


# ------------------------------------------------
# Figure 3 - correlation with increased aggregation
# ------------------------------------------------
# case incidence over time
case_plot <- 
  ggplot(data = dat_state)+
  geom_bar(aes(x = week,
               y = case_incidence,
               fill = "grey80"),
           position = "dodge",
           stat = "identity")+
  geom_line(
    aes(x = week,
        y = mean_sars2_conc_state_3w/12,
        color = "#ef8a47"),
    lwd = 1
  )+
  geom_line(
    aes(x = week,
        y = n_variants_no_thresh_3w_mean,
        color = "#e76254"),
    lwd = 1
  )+
  geom_line(
    aes(x = week,
        y = ntd_h_state_3w*10000,
        color = "#72bcd5"),
    lwd = 1
  )+
  geom_line(
    aes(x = week,
        y = ntd_pi_state_3w*10000,
        color = "#376795"),
    lwd = 1
  )+
  
  theme_dth_1+
  labs(
    title = "Case incidence (per 100,000)",
    x = "",
    y = "Cases"
  )+
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(size = 0.5))+
  scale_color_manual(
    values = c("#e76254" ="#e76254",
               "#72bcd5" = "#72bcd5",
               "#376795" = "#376795",
               "#ef8a47" = "#ef8a47"),
    labels = c(
      expression("S1 NTD "~ pi[ww]),
      expression("S1 NTD H"[ww]),
      "Freyja variant count",
      "Concentration"),
    
    name = ""
  )+
  scale_fill_manual(values = c("grey80" = "grey80"),
                    name = "",
                    label = "Cases")+
  scale_y_continuous(
    "COVID-19 cases per 100k", 
    sec.axis = sec_axis(~ ., name = expression(paste("Variant count/",
                                                     ~Pi[ww], "/", "H"[ww],  "/",
                                                     "Conc.",))
    )
  )+
  guides(color=guide_legend(nrow=2,byrow=TRUE))
case_plot

# hosp incidence
hosp_plot <- 
  ggplot(data = dat_state)+
  geom_bar(aes(x = week,
               y = hosp_incidence,
               fill = "grey80"),
           position = "dodge",
           stat = "identity")+
  geom_bar(aes(fill = "grey80",
               x = week,
               y = 0),
           position = "dodge",
           stat = "identity")+
  geom_line(
    aes(x = week,
        y = mean_sars2_conc_state_3w/100,
        color = "#ef8a47"),
    lwd = 1
  )+
  geom_line(
    aes(x = week,
        y = n_variants_no_thresh_3w_mean/10,
        color = "#e76254"),
    lwd = 1
  )+
  geom_line(
    aes(x = week,
        y = ntd_h_state_3w*1000,
        color = "#72bcd5"),
    lwd = 1
  )+
  geom_line(
    aes(x = week,
        y = ntd_pi_state_3w*1000,
        color = "#376795"),
    lwd = 1
  )+
  theme_dth_1+
  labs(
    title = "Hospitalizations (per 100,000)",
    x = ""
  )+
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(size = 0.5))+
  scale_color_manual(
    values = c("#e76254" ="#e76254",
               "#72bcd5" = "#72bcd5",
               "#376795" = "#376795",
               "#ef8a47" = "#ef8a47"),
    labels = c(
      expression("S1 NTD "~ pi[ww]),
      expression("S1 NTD H"[ww]),
      "Freyja variant count",
      "Concentration"),
    
    name = ""
  )+
  scale_fill_manual(values = c("grey80" = "grey80"),
                    name = "",
                    label = "Cases/Hospitalizations")+
  scale_y_continuous(
    "COVID-19 hospitalizations per 100k", 
    sec.axis = sec_axis(~ .*10, name = expression(paste("Variant count/",
                                                        ~Pi[ww], "/", "H"[ww], "/",
                                                        "Conc.",))
    )
  )+
  guides(color=guide_legend(nrow=2,byrow=TRUE))
hosp_plot

mylegend <- g_legend(hosp_plot)

plot_grid(case_plot + theme(legend.position = "none"),
          hosp_plot + theme(legend.position = "none"),
          mylegend,
          nrow = 3,
          align = "v",
          axis = "l",
          labels = c("A", "B"),
          rel_heights = c(3,3,1))


################################################################################

# Add the correlation plots

# correlation coefficient and p value
cor_1 <- cor.test(dat_state$ntd_pi_state_3w,
                  dat_state$case_incidence,
                  method = "spearman")
cor_1 <- round_2(cor_1$estimate)
signif_2(cor_1$p.value)

# ntd plot
pi_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = ntd_pi_state_3w,
             y = case_incidence))+
  geom_point(color = pal[6])+
  theme_dth_1+
  labs(
    title = "S1 NTD region",
    x = expression(paste( pi[ww])),
    y = "Case incidence",
    subtitle = expression(paste(rho, " = ", "0.92", ", p value < 0.001"))
  )+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_line(size = 0.5),
    plot.subtitle = element_text(hjust = 0.5))


# correlation coefficient and p value
cor_1 <- cor.test(dat_state$ntd_h_state_3w,
                  dat_state$case_incidence,
                  method = "spearman")
round_2(cor_1$estimate)
signif_2(cor_1$p.value)

# shannon diversity correlation plot
h_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = ntd_h_state_3w,
             y = case_incidence))+
  geom_point(color = pal[5])+
  theme_dth_1+
  labs(
    title = "S1 NTD region",
    x = expression(paste( "H"[ww], sep = "")),
    y = "Case incidence",
    subtitle = expression(paste(rho, " = ", "0.91", ", p value < 0.001"))
  )+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_line(size = 0.5),
    plot.subtitle = element_text(hjust = 0.5))

# correlation coefficient and p value
cor_1 <- cor.test(dat_state$n_variants_no_thresh_3w_mean,
                  dat_state$case_incidence,
                  method = "spearman")
round_2(cor_1$estimate)
signif_2(cor_1$p.value)

# freyja variant count
var_count_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = n_variants_no_thresh_3w_mean,
             y = case_incidence))+
  geom_point(color = pal[1])+
  theme_dth_1+
  labs(
    title = "Freyja variant count",
    x = "Variant count",
    y = "Case incidence",
    subtitle = expression(paste(rho, " = ", "0.85", ", p value < 0.001"))
  )+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_line(size = 0.5),
    plot.subtitle = element_text(hjust = 0.5))

# correlation coefficient and p value
cor_1 <- cor.test(dat_state$mean_sars2_conc_state_3w,
                  dat_state$case_incidence,
                  method = "spearman")
round_2(cor_1$estimate)
signif_2(cor_1$p.value)

# concentration
conc_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = mean_sars2_conc_state_3w,
             y = case_incidence))+
  geom_point(color = pal[2])+
  theme_dth_1+
  labs(
    title = "Concentration",
    x = expression(paste("Conc, copies/", "m","L", sep = "")),
    y = "Case incidence",
    subtitle = expression(paste(rho, " = ", "0.73", ", p value < 0.001"))
  )+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_line(size = 0.5),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))

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
    labels = c("B", "C", "D", "E")
  )

# add title to panel
case_panel <- plot_grid(title, panel, ncol=1, rel_heights=c(0.1, 1))

# add hospitalization correlations as panels E through H

# correlation coefficient and p value
cor_1 <- cor.test(dat_state$ntd_pi_state_3w,
                  dat_state$hosp_incidence,
                  method = "spearman")
round_2(cor_1$estimate)
signif_2(cor_1$p.value)

# ntd plot
pi_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = ntd_pi_state_3w,
             y = hosp_incidence))+
  geom_point(color = pal[6])+
  theme_dth_1+
  labs(
    title = "S1 NTD region",
    x = expression(paste( pi[ww])),
    y = "Hospitalization incidence",
    subtitle = expression(paste(rho, " = ", "0.88", ", p value < 0.001"))
  )+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_line(size = 0.5),
    plot.subtitle = element_text(hjust = 0.5))

# correlation coefficient and p value
cor_1 <- cor.test(dat_state$ntd_h_state_3w,
                  dat_state$hosp_incidence,
                  method = "spearman")
round_2(cor_1$estimate)
signif_2(cor_1$p.value)

# shannon diversity correlation plot
h_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = ntd_h_state_3w,
             y = hosp_incidence))+
  geom_point(color = pal[5])+
  theme_dth_1+
  labs(
    title = "S1 NTD region",
    x = expression(paste( "H"[ww], sep = "")),
    y = "Hospitalization incidence",
    subtitle = expression(paste(rho, " = ", "0.87", ", p value < 0.001"))
  )+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_line(size = 0.5),
    plot.subtitle = element_text(hjust = 0.5))

# correlation coefficient and p value
cor_1 <- cor.test(dat_state$n_variants_no_thresh_3w_mean,
                  dat_state$hosp_incidence,
                  method = "spearman")
round_2(cor_1$estimate)
signif_2(cor_1$p.value)

# freyja variant count
var_count_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = n_variants_no_thresh_3w_mean,
             y = hosp_incidence))+
  geom_point(color = pal[1])+
  theme_dth_1+
  labs(
    title = "Freyja variant count",
    x = "Variant count",
    y = "Hospitalization incidence",
    subtitle = expression(paste(rho, " = ", "0.81", ", p value < 0.001"))
  )+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_line(size = 0.5),
    plot.subtitle = element_text(hjust = 0.5))

# correlation coefficient and p value
cor_1 <- cor.test(dat_state$mean_sars2_conc_state_3w,
                  dat_state$hosp_incidence,
                  method = "spearman")
round_2(cor_1$estimate)
signif_2(cor_1$p.value)

# concentration
conc_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = mean_sars2_conc_state_3w,
             y = hosp_incidence))+
  geom_point(color = pal[2])+
  theme_dth_1+
  labs(
    title = "Concentration",
    x = expression(paste("Conc, copies/", "m","L", sep = "")),
    y = "Hospitalization incidence",
    subtitle = expression(paste(rho, " = ", "0.72", ", p value < 0.001"))
  )+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_line(size = 0.5),
    plot.subtitle = element_text(hjust = 0.5))

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
    labels = c("G", "H", "I", "J")
  )

# add title to panel
hosp_panel <- plot_grid(title, panel, ncol=1, rel_heights=c(0.1, 1))

# combine panels
plot_grid(
  case_panel,
  hosp_panel,
  nrow = 2
)

# combine them both
plot_grid(case_plot + theme(legend.position = "none"),
          case_panel,
          hosp_plot + theme(legend.position = "none"),
          hosp_panel,
          mylegend,
          nrow = 3,
          align = "v",
          axis = "l",
          labels = c("A", "", "F", ""),
          rel_heights = c(3,3,1),
          rel_widths = c(2,1.5))


# save
png("Figures/Figure 4.png",
    units = "in",
    width = 14, height = 12,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(case_plot + theme(legend.position = "none"),
          case_panel,
          hosp_plot + theme(legend.position = "none"),
          hosp_panel,
          mylegend,
          nrow = 3,
          align = "v",
          axis = "l",
          labels = c("A", "", "F", ""),
          rel_heights = c(3,3,1),
          rel_widths = c(2,1.5))
showtext_end()
dev.off()
