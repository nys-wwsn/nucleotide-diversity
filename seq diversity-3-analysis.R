###########
###########
#### GENOMEWIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-08-28

# Main analysis - figures and tables for sequencing diversity paper

# This paper uses several functions and the docstring package to document
# the function usage. You can view help pages for each function by writing in 
# the console ?function_name

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

# load county spatial data
ny_counties <- st_read("data/Counties-Shoreline/Counties_Shoreline.shp") %>%
  mutate(county = NAME)
ny_counties$county[ny_counties$county == "St Lawrence"] <- "St Lawrence"

# -----------------------------------------------------------------------------


# --------------------------------------------------------------------

# --------------
# Figure 1
# --------------
# MAP OF THE SAMPLING SITES AND POSITIVE DETECTION RATE

# summarize by pcr lab
detection_rate <- dat_sewershed %>%
  filter(week >= ymd(min_week) & week <= ymd(max_week)) %>%
  filter(!is.na(genomewide_pi)) %>%
  mutate(sars_positive = ifelse(mean_sars2_conc > 1, 1, 0)) %>%
  group_by(facility_id, wwtp_latitude, wwtp_longitude,capacity_mgd
  ) %>%
  summarize(`Number of positive samples` = sum(sars_positive, na.rm = TRUE),
            `Number of samples` = n()
  ) %>%
  ungroup() %>%
  mutate(`Detection rate` = `Number of positive samples`/`Number of samples`) 

# plot by site

# make the facility data into points
detection_rate <- detection_rate %>%
  st_as_sf(coords =  c("wwtp_longitude", "wwtp_latitude"),
           crs = "+proj=longlat +zone=18 +datum=WGS84") %>% 
  st_transform(st_crs(ny_counties))


# map
map_dr <- 
  ggplot()+
  geom_sf(data = ny_counties,
          #aes(fill = pcr_lab),
          fill = "white",
          color = "grey65",
          lwd = 1
          
  )+
  geom_sf(data = detection_rate,
          aes(color = as.numeric(`Detection rate`),
              size = capacity_mgd),
          alpha = 0.7)+
  geom_sf(data = detection_rate,
          aes(
            size = capacity_mgd),
          shape = 1)+
  theme_dth_maps +
  scale_color_gradientn(colors = rev(MetBrewer::met.brewer("Hiroshige")),
                        name = "Detection\nrate")+
  theme(legend.position.inside = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.direction = "horizontal")+
  labs(size = "WWTP size\n(mgd)",
       title = "SARS-CoV-2 detection rate in wastewater")+
  theme(legend.position = c(0.6, 0.3),
        legend.background = element_blank(),
        legend.box.just = "left")+
  scale_size(range = c(0, 20))

map_dr

# sample collected over time and samples sequenced stacked bar
samples_time <- dat_sewershed %>%
  filter(week >= ymd(min_week) & week <= ymd(max_week)) %>%
  filter(!is.na(mean_sars2_conc)) %>%
  group_by(week) %>%
  summarise(samples_collected = sum(conc_samples, na.rm = TRUE)
  ) %>%
  ungroup() 

# sequenced samples over time
seq_time <- dat_sewershed  %>%
  filter(week >= ymd(min_week) & week <= ymd(max_week)) %>%
  filter(!is.na(genomewide_pi)) %>%
  group_by(week) %>%
  summarise(samples_sequenced = sum(samples, na.rm = TRUE)
  ) %>%
  ungroup() 

# join
samples_time <- left_join(samples_time, seq_time, by = c("week"))
samples_time$samples_no_seq <- 
  samples_time$samples_collected - samples_time$samples_sequenced

samples_time <- samples_time %>%
  pivot_longer(cols = c(samples_sequenced, samples_no_seq))

pal_1 <- met.brewer(name = "Hiroshige", n = 6)

# combine into one panel figure
samples_plot <- 
  ggplot(data = samples_time) +
  geom_bar(position = "stack",
           stat = "identity",
           aes(x = week, y = value,
               fill = name))+
  theme_dth_1+
  theme(legend.position = "bottom")+
  labs(x = "",
       y = "Samples collected",
       fill = "",
       title = "Samples collected over time")+
  scale_fill_manual(values = c("samples_no_seq" = pal_1[1],
                               "samples_sequenced" = pal_1[6]),
                    labels = c("Samples not sequenced",
                               "Samples sequenced"))

samples_plot 



# panel
ggarrange(map_dr, samples_plot , nrow = 1, labels = c("A","B"))

# save map
png("figures and tables/Figure 1 -Detection rate map and samples over time.png",
    units = "in",
    width = 12, height = 6,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
ggarrange(map_dr, samples_plot , nrow = 1, labels = c("A","B"))
showtext_end()
dev.off()


# ------------------------------------------------
# Figure 2 - Single sample figures
# ------------------------------------------------

# bird island
bird_pi <- readRDS("data/erie_NY029028410A.rds")%>%
  mutate(week = floor_date(ymd(date), unit = "weeks", week_start = 7)
  )

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
                group = date),
            color = "grey"
  )+
  geom_line(data = bird_pi %>%
              filter(ymd(date) >= "2023-11-01"&
                       ymd(date) <= "2024-02-01") %>%
              group_by(center) %>%
              mutate(mean_h = mean(avg_h, na.rm = TRUE)
              ),
            aes(x = center, y = mean_h),
            color = "black",
            lwd = 1
  )+
  ylim(0,0.03)+
  theme_dth_1+
  labs(x = "",
       y = expression("Mean H"[ww]))+
  theme(plot.margin=unit(c(0.9,1,-0.19,1), "cm"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_x_continuous(expand = c(0,0))

pi_time <- ggplot()+
  geom_line(data = bird_pi %>%
              filter(ymd(date) >= "2023-11-01"&
                       ymd(date) <= "2024-02-01"),
            aes(x = center,
                y = avg_pi,
                group = date),
            color = "grey"
  )+
  geom_line(data = bird_pi %>%
              filter(ymd(date) >= "2023-11-01"&
                       ymd(date) <= "2024-02-01") %>%
              group_by(center) %>%
              mutate(mean_pi = mean(avg_pi, na.rm = TRUE)
              ),
            aes(x = center, y = mean_pi),
            color = "black",
            lwd = 1
  )+
  ylim(0,0.015)+
  theme_dth_1+
  labs(x = "",
       y = expression(paste("Mean ",pi[ww])))+
  theme(plot.margin=unit(c(0.9,1,-0.19,1), "cm"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_x_continuous(expand = c(0,0))

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
png("figures and tables/Figure 2 - pi and h per sample comparison.png",
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


# ------------------------------------------------
# Figure 3 - correlation with increased aggregation
# ------------------------------------------------
# case incidence over time
case_plot <- 
  ggplot(data = dat_state)+
  geom_bar(aes(x = week,
               y = case_incidence,
               fill = "grey60"),
           position = "dodge",
           stat = "identity")+
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
        legend.background = element_blank())+
  scale_color_manual(
    values = c("#e76254" ="#e76254",
               "#72bcd5" = "#72bcd5",
               "#376795" = "#376795"),
    labels = c(
      expression("S1 NTD "~ pi[ww]),
      expression("S1 NTD H"[ww]),
      "Freyja variant count"),
    
    name = ""
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    name = "",
                    label = "Cases")+
  scale_y_continuous(
    "COVID-19 cases per 100k", 
    sec.axis = sec_axis(~ ., name = expression(paste("Variant count/",
                                                     ~Pi[ww], "/", "H"[ww], ))
    )
  )
case_plot

# hosp incidence
hosp_plot <- 
  ggplot(data = dat_state)+
  geom_bar(aes(x = week,
               y = hosp_incidence,
               fill = "grey70"),
           position = "dodge",
           stat = "identity")+
  geom_bar(aes(fill = "grey60",
               x = week,
               y = 0),
           position = "dodge",
           stat = "identity")+
  geom_line(
    aes(x = week,
        y = n_variants_no_thresh_3w_mean/8,
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
        legend.background = element_blank())+
  scale_color_manual(
    values = c("#e76254" ="#e76254",
               "#72bcd5" = "#72bcd5",
               "#376795" = "#376795"),
    labels = c(
      expression("S1 NTD "~ pi[ww]),
      expression("S1 NTD H"[ww]),
      "Freyja variant count"),
    
    name = ""
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    name = "",
                    label = "Cases/Hospitalizations")+
  scale_y_continuous(
    "COVID-19 hospitalizations per 100k", 
    sec.axis = sec_axis(~ .*8, name = expression(paste("Variant count/",
                                                       ~Pi[ww], "/", "H"[ww], ))
    )
  )
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

# save
png("figures and tables/Figure 3 - time series statewide samples_v3.png",
    units = "in",
    width = 6, height = 8,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(case_plot + theme(legend.position = "none"),
          hosp_plot + theme(legend.position = "none"),
          mylegend,
          nrow = 3,
          align = "v",
          axis = "l",
          labels = c("A", "B"),
          rel_heights = c(3,3,1))
showtext_end()
dev.off()


# -----------------------------------------------------------------------------
# Figure 4 - correlation scatterplots for all 3 measures plus concentration
# -----------------------------------------------------------------------------

# ntd plot
pi_cor_plot <- 
  ggplot(data = dat_state,
         aes(x = ntd_pi_state_3w,
             y = case_incidence))+
  geom_point(color = pal[6])+
  stat_cor(method = "spearman")+
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
  stat_cor(method = "spearman")+
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
  stat_cor(method = "spearman")+
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
  stat_cor(method = "spearman")+
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
  stat_cor(method = "spearman")+
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
  stat_cor(method = "spearman")+
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
  stat_cor(method = "spearman")+
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
  stat_cor(method = "spearman")+
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
png("figures and tables/Figure 4 - correlation scatterplots.png",
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


# --------------------------------------------------------------------------
# Figure 5 - lag and lead times for all measures
# --------------------------------------------------------------------------

# case
lag_case <- lag_lead_function(dataframe = dat_state,
                              columns = c(# pi
                                "ntd_pi_state_3w",
                                # H
                                "ntd_h_state_3w",
                                # var count
                                "n_variants_no_thresh_3w_mean",
                                # conc
                                "mean_sars2_conc_state_3w"),
                              lags = seq(1:10))

lag_case_plot <- ggplot(data = lag_case ,
                        aes(x = lag, y = spearman_cor, color = name))+
  geom_point(shape = 1, size = 2)+
  geom_point(shape = 17,
             data = lag_case %>%
               group_by(name) %>%
               filter(spearman_cor == max(spearman_cor)),
             size = 3
  )+
  geom_line()+
  ylim(0,1)+
  theme_dth_1+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_line(color = "gray80"),
        panel.grid.minor.y = element_blank())+
  scale_color_manual(values = c(
    "ntd_pi_state_3w"=pal[6],
    "ntd_h_state_3w"=pal[5],
    "n_variants_no_thresh_3w_mean"=pal[1],
    "mean_sars2_conc_state_3w"=pal[2]
  ),
  name = "",
  labels = c(
    "ntd_pi_state_3w"=expression("S1 NTD "~ pi[ww]),
    "ntd_h_state_3w"=expression( "S1 NTD H"[ww]),
    "n_variants_no_thresh_3w_mean"="Freyja Variant counts",
    "mean_sars2_conc_state_3w"="Concentration"
  )
  )+
  scale_x_continuous(breaks = seq(-5, 5, by = 1),
                     limits = c(-5,5))+
  labs(x = "Lead time (weeks)",
       y = "spearman correlation\ncoefficient",
       title= "Case incidence ~ Wastewater measure for\nstatewide data")+
  annotate("text",
           x = -3, y = 0.09,
           label = "signal lags incidence")+
  annotate("text",
           x = 3, y = 0.09,
           label = "signal leads incidence")+
  guides(color=guide_legend(nrow=2,byrow=TRUE))
lag_case_plot

# hosp

lag_hosp <- lag_lead_function_hosp(dataframe = dat_state,
                                   columns = c(# pi
                                     "ntd_pi_state_3w",
                                     # H
                                     "ntd_h_state_3w",
                                     # var count
                                     "n_variants_no_thresh_3w_mean",
                                     # conc
                                     "mean_sars2_conc_state_3w"),
                                   lags = seq(1:10))

lag_hosp_plot <- ggplot(data = lag_hosp ,
                        aes(x = lag, y = spearman_cor, color = name))+
  geom_point(shape = 1, size = 2)+
  geom_point(shape = 17,
             data = lag_hosp %>%
               group_by(name) %>%
               filter(spearman_cor == max(spearman_cor)),
             size = 3
  )+
  geom_line()+
  ylim(0,1)+
  theme_dth_1+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_line(color = "gray80"),
        panel.grid.minor.y = element_blank())+
  scale_color_manual(values = c(
    "ntd_pi_state_3w"=pal[6],
    "ntd_h_state_3w"=pal[5],
    "n_variants_no_thresh_3w_mean"=pal[1],
    "mean_sars2_conc_state_3w"=pal[2]
  ),
  name = "",
  labels = c(
    "ntd_pi_state_3w"=expression("S1 NTD "~ pi[ww]),
    "ntd_h_state_3w"=expression( "S1 NTD H"[ww]),
    "n_variants_no_thresh_3w_mean"="Freyja Variant counts",
    "mean_sars2_conc_state_3w"="Concentration"
  )
  )+
  scale_x_continuous(breaks = seq(-5, 5, by = 1),
                     limits = c(-5,5))+
  labs(x = "Lead time (weeks)",
       y = "spearman correlation\ncoefficient",
       title= "Hospitalization incidence ~ Wastewater measure for\nstatewide data")+
  annotate("text",
           x = -3, y = 0.09,
           label = "signal lags incidence")+
  annotate("text",
           x = 3, y = 0.09,
           label = "signal leads incidence")


# extract legend function
mylegend <- g_legend(lag_case_plot)


# save
png("figures and tables/Figure 5 - lag plot.png",
    units = "in",
    width = 10, height = 6.5,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(
  plot_grid(lag_case_plot+ theme(legend.position = "none"),
            lag_hosp_plot + theme(legend.position = "none"),
            labels = c("A", "B"), nrow = 1, ncol = 2),
  plot_grid(NULL, mylegend, NULL, nrow = 1, rel_widths = c(1, 0.5, 1)),
  nrow = 2,
  rel_heights = c(5,1)
)
dev.off()

