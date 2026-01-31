###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-25

# Main analysis - Figure 2

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

# bird island
bird_pi <- readRDS("data/other data/erie_NY029028410A.rds")%>%
  mutate(week = floor_date(ymd(date), unit = "weeks", week_start = 7)
  )
# -----------------------------------------------------------------------------
# ------------------------------------------------
# Figure 2 - Single sample figures
# ------------------------------------------------


# make a tiled heatmap of diversity by base over time
h_heatmap <- ggplot()+
  geom_tile(data = bird_pi %>%
              filter(ymd(date) >= "2023-11-01"&
                       ymd(date) <= "2024-02-01"),
            aes(x = center,
                y = floor_date(ymd(date), 
                               unit = "weeks", 
                               week_start = 7),
                fill = avg_h))+
  theme_dth_1+
  theme(legend.position = "top",
        plot.margin=unit(c(0.9,1,-0.19,1), "cm"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_gradientn(colors=rev(met.brewer("Hiroshige")),
                       name = expression("Mean H"[ww]),
                       limits=c(0,0.023),
                       breaks = c(0,0.01,0.02))+
  labs(x = "",
       y = "")+
  scale_y_date(date_labels = "%b %Y",
               date_breaks = "1 month")+
  scale_x_continuous(expand = c(0,0))

h_time <- ggplot()+
  geom_line(data = bird_pi %>%
              filter(ymd(date) >= "2023-11-01"&
                       ymd(date) <= "2024-02-01"),
            aes(x = center,
                y = avg_h,
                group = date,
                color = "grey")
  )+
  geom_line(data = bird_pi %>%
              filter(ymd(date) >= "2023-11-01"&
                       ymd(date) <= "2024-02-01") %>%
              group_by(center) %>%
              mutate(mean_h = mean(avg_h, na.rm = TRUE)
              ),
            aes(x = center, y = mean_h, color = "black"),
            lwd = 1
  )+
  ylim(0,0.03)+
  theme_dth_1+
  labs(x = "",
       y = expression("Mean H"[ww]))+
  theme(plot.margin=unit(c(0.9,1,-0.19,1), "cm"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top")+
  scale_x_continuous(expand = c(0,0))+
  scale_color_manual(values = c("black" = "black",
                                "grey" = "grey"),
                     labels = c(expression(paste("Mean H"[ww])),
                                expression(paste("Sample H"[ww]))),
                     name = "Values"
                     )

pi_time <- ggplot()+
  geom_line(data = bird_pi %>%
              filter(ymd(date) >= "2023-11-01"&
                       ymd(date) <= "2024-02-01"),
            aes(x = center,
                y = avg_pi,
                group = date,
                color = "grey")
  )+
  geom_line(data = bird_pi %>%
              filter(ymd(date) >= "2023-11-01"&
                       ymd(date) <= "2024-02-01") %>%
              group_by(center) %>%
              mutate(mean_pi = mean(avg_pi, na.rm = TRUE)
              ),
            aes(x = center, y = mean_pi,
                color = "black"),
            lwd = 1
  )+
  ylim(0,0.015)+
  theme_dth_1+
  labs(x = "",
       y = expression(paste("Mean ",pi[ww])))+
  theme(plot.margin=unit(c(0.9,1,-0.19,1), "cm"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top")+
  scale_x_continuous(expand = c(0,0))+
  scale_color_manual(values = c("black" = "black",
                                "grey" = "grey"),
                     labels = c(expression(paste("Mean ",Pi[ww])),
                                expression(paste("Sample ",Pi[ww]))),
                     name = "Values"
  )

pi_heatmap <- ggplot()+
  geom_tile(data = bird_pi %>%
              filter(ymd(date) >= "2023-11-01"&
                       ymd(date) <= "2024-02-01"),
            aes(x = center,
                y = floor_date(ymd(date), 
                               unit = "weeks", 
                               week_start = 7),
                fill = avg_pi))+
  theme_dth_1+
  theme(legend.position = "top",
        plot.margin=unit(c(0.9,1,-0.19,1), "cm"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_gradientn(colors=rev(met.brewer("Hiroshige")),
                       name = expression(paste("Mean ",pi[ww])),
                       limits=c(0,0.015),
                       breaks = c(0,0.007,0.014)
  )+
  labs(x = "",
       y = "")+
  scale_y_date(date_labels = "%b %Y",
               date_breaks = "1 month")+
  scale_x_continuous(expand = c(0,0))


# base pair plots
# plot for sample data with sequencing locations
plot_bp_labels <-ggplot()+
  annotate("rect",
           xmin = 0, xmax = 30000,
           ymin = 0, ymax = 0.5,
           alpha = 0.5,
           fill = "white",
           color = "black")+
  annotate("text",
           label = "ORF region",
           x = 10000,
           y = 0.25
  )+
  annotate("rect",
           xmin = 10063, xmax = 11842,
           ymin = 0.5, ymax = 1,
           fill = pal[5])+
  annotate("text",
           label = "NSPs 5+6",
           x = 6500,
           y = 0.75,
           size = 4,
           color = pal[5]
  )+
  annotate("rect",
           xmin = 20661, xmax = 21549,
           ymin = 0.5, ymax = 1,
           fill = pal[1])+
  annotate("text",
           label = "2'O-Mtase ",
           x = 17500,
           y = 0.75,
           size = 3.5,
           color = pal[1]
  )+
  annotate("rect",
           xmin = 21598, xmax = 22474,
           ymin = 0.5, ymax = 1,
           fill = pal[6])+
  annotate("text",
           label = "NTD",
           x = 24500,
           y = 0.89,
           color = pal[6],
           size = 3
  )+
  annotate("rect",
           xmin = 22516, xmax = 23185,
           ymin = 0.5, ymax = 1,
           fill = pal[2])+
  annotate("text",
           label = "RBD",
           x = 24500,
           y = 0.65,
           color = pal[2],
           size = 3
  )+
  annotate("rect",
           xmin = 21563, xmax = 25384,
           ymin = 0, ymax = 0.5,
           fill = "grey")+
  annotate("text",
           label = "Spike",
           x = 23250,
           y = 0.25,
           color = "black",
           size = 4
  )+
  labs(x = "Base pair position",
       y = "Genome")+
  theme_dth_1+
  theme(axis.text.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin=unit(c(-0.19,1,0.5,1), "cm"),
        plot.background = element_blank(),
        panel.border = element_blank())+
  scale_x_continuous(expand = c(0,0))
plot_bp_labels


# save plots
png("Figures/Figure 2.png",
    units = "in",
    width = 9.5, height = 11,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
cowplot::set_null_device("png")
plot_grid(
  pi_heatmap, h_heatmap,
  plot_bp_labels, plot_bp_labels,
  pi_time, h_time, 
  plot_bp_labels, plot_bp_labels,
  
  nrow = 4,
  align = "v",
  axis = "l",
  rel_heights = c(5,1,5,1),
  labels = c("A", "B", "", "", "C", "D", "", "")
)
dev.off()