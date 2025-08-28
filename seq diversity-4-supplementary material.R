###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-08-28

# Supplemental figures and tables for the sars-cov-2 sequencing diversity 
# paper

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
depth_1 <- read.csv("data/depths/depth_stats_20250424.txt") %>%
  tidyr::separate(
    sample, into = c("date", "cdc_id"), sep = 8
    ) %>%
  mutate(date = lubridate::ymd(date),
         mean_depth = mean) %>%
  select(date, cdc_id, mean_depth)

# mean depth post july 2024
depth_2 <- read.csv("data/depths/genomwide_depth_2025-07-02.csv") %>%
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
             path = paste("figures and tables/",
                          "Table S1 - genome-wide pi descriptive statistics.docx",
                          sep = "")
)



# ---------------------------------------
# Supplemental Figure S1 - depth and ct
# ---------------------------------------

# Depth with conc and Depth with ct - 4 figure panel with correlations and 
# time series.

depth_ct <- 
  ggplot(data = dat_sewershed,
         aes(x = mean_depth,
             y = mean_ct))+
  geom_point(alpha = 0.5, color = "tomato")+
  stat_cor(method = "spearman")+
  theme_dth_1+
  labs(x = "Mean depth",
       y = "Mean Ct",
       title = "Ct ~ depth")
depth_ct

# ct over time

ct_time <- 
  ggplot(data = dat_sewershed,
         aes(x = week,
             y = mean_ct))+
  geom_point(alpha = 0.5, color = "tomato")+
  theme_dth_1+
  labs(x = "",
       y = "Mean Ct",
       title = "Ct over time")
ct_time

# depth by concentration
depth_conc <- 
  ggplot(data = dat_sewershed,
         aes(x = mean_depth,
             y = mean_sars2_conc))+
  geom_point(alpha = 0.5, color = "orange")+
  stat_cor(method = "spearman")+
  theme_dth_1+
  labs(x = "Mean depth",
       y = "Conc, copies/mL",
       title = "Concentration ~ depth")
depth_conc

# ct over time

conc_time <- 
  ggplot(data = dat_sewershed,
         aes(x = week,
             y = mean_sars2_conc))+
  geom_point(alpha = 0.5, color = "orange")+
  theme_dth_1+
  labs(x = "",
       y = expression(paste("Conc, copies/mL",
                            )),
       title = "Concentration over time")
conc_time

plot_grid(depth_ct, depth_conc,
             ct_time, conc_time,
             nrow = 2,
             labels = c("A", "B", "C", "D"))

# save
png("figures and tables/Figure S1-depth and ct.png",
    units = "in",
    width = 12, height = 6,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(depth_ct, depth_conc,
          ct_time, conc_time,
          nrow = 2,
          labels = c("A", "B", "C", "D"))
showtext_end()
dev.off()

# --------------------------------------------------------------------
# Figure S2 - county level ntd concentration and case correlation
# --------------------------------------------------------------------

# ntd pi over time by county
ntd_pi_plot <-
  ggplot(data = dat_county)+
  geom_point(aes(x = week,
                 y= ntd_pi_county_3w),
             alpha = 0.6,
             color = "grey60")+
  theme_dth_1+
  labs(
    title = expression(paste( pi[ww])),
    x = "",
    y = expression(paste( pi[ww]))
  )

ntd_pi_plot

# ntd h over time by county
ntd_h_plot <-
  ggplot(data = dat_county)+
  geom_point(aes(x = week,
                 y= ntd_h_county_3w),
             alpha = 0.6,
             color = "grey60")+
  theme_dth_1+
  labs(
    title = expression(paste("H"[ww])),
    x = "",
    y = expression(paste("H"[ww]))
  )

# variant count over time by county
var_count_plot <- 
  ggplot(data = dat_county)+
  geom_point(aes(x = week,
                 y= n_variants_no_thresh_3w_mean,
                 fill = "grey60"),
             alpha = 0.6,
             color = "grey60")+
  theme_dth_1+
  labs(
    title = "Freyja variant counts",
    x = "",
    y = "Variant counts"
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    name = "",
                    labels = c("County weighted average"))+
  theme(legend.position = "bottom",
        legend.background = element_blank())
var_count_plot

# case incidence over time
case_plot <- 
  ggplot(data = dat_state)+
  geom_bar(aes(x = week,
               y = case_incidence),
           position = "dodge",
           stat = "identity")+
  theme_dth_1+
  labs(
    title = "Case incidence (per 100,000)",
    x = "",
    y = "Cases"
  )
mylegend <- g_legend(var_count_plot)


# put the plots in a panel
plot_grid(ntd_pi_plot, 
          ntd_h_plot, 
          var_count_plot+ 
            theme(legend.position = "none"),
          case_plot, 
          mylegend,
          nrow = 5,
          align = "v",
          axis = "l",
          labels = c("A", "B", "C", "D"),
          rel_heights = c(3,3,3,3,1))

# save
png("figures and tables/Figure S2 - time series county average.png",
    units = "in",
    width = 8.5, height = 11,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(ntd_pi_plot, 
          ntd_h_plot, 
          var_count_plot+ 
            theme(legend.position = "none"),
          case_plot, 
          mylegend,
          nrow = 5,
          align = "v",
          axis = "l",
          labels = c("A", "B", "C", "D"),
          rel_heights = c(3,3,3,3,1))
showtext_end()
dev.off()

# --------------------------------------------------------------------
# Figure S3 - region level ntd concentration and case correlation
# --------------------------------------------------------------------

# ntd pi over time by region
ntd_pi_plot <-
  ggplot(data = dat_region)+
  geom_point(aes(x = week,
                 y= ntd_pi_region_3w),
             alpha = 0.6,
             color = "grey60")+
  theme_dth_1+
  labs(
    title = expression(paste( pi[ww])),
    x = "",
    y = expression(paste( pi[ww]))
  )

ntd_pi_plot

# ntd h over time by region
ntd_h_plot <-
  ggplot(data = dat_region)+
  geom_point(aes(x = week,
                 y= ntd_h_region_3w),
             alpha = 0.6,
             color = "grey60")+
  theme_dth_1+
  labs(
    title = expression(paste("H"[ww])),
    x = "",
    y = expression(paste("H"[ww]))
  )

# variant count over time by region
var_count_plot <- 
  ggplot(data = dat_region)+
  geom_point(aes(x = week,
                 y= n_variants_no_thresh_3w_mean,
                 fill = "grey60"),
             alpha = 0.6,
             color = "grey60")+
  theme_dth_1+
  labs(
    title = "Freyja variant counts",
    x = "",
    y = "Variant counts"
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    name = "",
                    labels = c("Regional weighted average"))+
  theme(legend.position = "bottom",
        legend.background = element_blank())
var_count_plot

# case incidence over time
case_plot <- 
  ggplot(data = dat_state)+
  geom_bar(aes(x = week,
               y = case_incidence),
           position = "dodge",
           stat = "identity")+
  theme_dth_1+
  labs(
    title = "Case incidence (per 100,000)",
    x = "",
    y = "Cases"
  )
mylegend <- g_legend(var_count_plot)


# put the plots in a panel
plot_grid(ntd_pi_plot, 
          ntd_h_plot, 
          var_count_plot+ 
            theme(legend.position = "none"),
          case_plot, 
          mylegend,
          nrow = 5,
          align = "v",
          axis = "l",
          labels = c("A", "B", "C", "D"),
          rel_heights = c(3,3,3,3,1))

# save
png("figures and tables/Figure S3 - time series region average.png",
    units = "in",
    width = 8.5, height = 11,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(ntd_pi_plot, 
          ntd_h_plot, 
          var_count_plot+ 
            theme(legend.position = "none"),
          case_plot, 
          mylegend,
          nrow = 5,
          align = "v",
          axis = "l",
          labels = c("A", "B", "C", "D"),
          rel_heights = c(3,3,3,3,1))
showtext_end()
dev.off()

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
png("figures and tables/Figure S4 - pi time series and correlation.png",
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


# -----------------------------------------------
# Figure S5 - H and genome regions time series
# -----------------------------------------------
# SHANNON

# Statewide time series figure, all genome regions, panel for NYC
dat_state_long <- dat_state %>%
  pivot_longer(cols = c(
    genomewide_h_state_3w,
    ntd_h_state_3w,
    cov_mt_2_h_state_3w,
    orf_h_state_3w,
    spike_h_state_3w,
    rbd_h_state_3w
  ))%>%
  mutate(name_factor = factor(name,
                              levels = c("genomewide_h_state_3w",
                                         "spike_h_state_3w",
                                         "orf_h_state_3w",
                                         "ntd_h_state_3w",
                                         "rbd_h_state_3w",
                                         "cov_mt_2_h_state_3w"),
                              labels = c("Genome",
                                         "Spike",
                                         "ORF NSPs 5 and 6",
                                         "S1 NTD",
                                         "S1 RBD",
                                         "CoV-MT-2"))
  )

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
       title = "Genome region Shannon's H values")+
  scale_y_continuous(
    "COVID-19 cases per 100k", 
    sec.axis = sec_axis(~ ./ 10000, name = "Shannon's H")
  )

time_series_1

# figure with correlations by genome region
cor_plot <- dat_state %>%
  pivot_longer(cols = c("genomewide_h_state_3w",
                        "spike_h_state_3w",
                        "orf_h_state_3w",
                        "ntd_h_state_3w",
                        "rbd_h_state_3w",
                        "cov_mt_2_h_state_3w")
  ) %>%
  mutate(name_factor = factor(name,
                              levels = c("genomewide_h_state_3w",
                                         "spike_h_state_3w",
                                         "orf_h_state_3w",
                                         "ntd_h_state_3w",
                                         "rbd_h_state_3w",
                                         "cov_mt_2_h_state_3w"),
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
  labs(title = "Spearman correlation:\nShannon's H ~ case incidence",
       x = "Shannon's H",
       y = "")

cor_plot


cowplot::plot_grid(time_series_1, cor_plot, nrow = 1,
                   rel_widths = c(2,1),
                   align = 'h', axis = 'tb')

# combined plot for pi
# save
png("figures and tables/Figure S5 - H time series and correlation.png",
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

# --------------------------------------
# Table S2  spearman correlation values
# --------------------------------------


# correlation for each genome region predictor for pi and H
row_1 <- cor_ci_function(
  dataframe = dat_state,
  Region = "Genome-wide",
  predictor_name = "genomewide_pi_state_3w",
  predictor_group = "Pi"
)

row_2 <- cor_ci_function(
  dataframe = dat_state,
  Region = "NSPs 5 and 6",
  predictor_name = "orf_pi_state_3w",
  predictor_group = "Pi"
)

row_3 <- cor_ci_function(
  dataframe = dat_state,
  Region = "2' O-Mtase",
  predictor_name = "cov_mt_2_pi_state_3w",
  predictor_group = "Pi"
)

row_4 <- cor_ci_function(
  dataframe = dat_state,
  Region = "Spike Protein",
  predictor_name = "spike_pi_state_3w",
  predictor_group = "Pi"
)

row_5 <- cor_ci_function(
  dataframe = dat_state,
  Region = "S1 NTD",
  predictor_name = "ntd_pi_state_3w",
  predictor_group = "Pi"
)

row_6 <- cor_ci_function(
  dataframe = dat_state,
  Region = "S1 RBD",
  predictor_name = "rbd_pi_state_3w",
  predictor_group = "Pi"
)

# repeat for h values
row_7 <- cor_ci_function(
  dataframe = dat_state,
  Region = "Genome-wide",
  predictor_name = "genomewide_h_state_3w",
  predictor_group = "H"
)

row_8 <- cor_ci_function(
  dataframe = dat_state,
  Region = "NSPs 5 and 6",
  predictor_name = "orf_h_state_3w",
  predictor_group = "H"
)

row_9 <- cor_ci_function(
  dataframe = dat_state,
  Region = "2' O-Mtase",
  predictor_name = "cov_mt_2_h_state_3w",
  predictor_group = "H"
)

row_10 <- cor_ci_function(
  dataframe = dat_state,
  Region = "Spike Protein",
  predictor_name = "spike_h_state_3w",
  predictor_group = "H"
)

row_11 <- cor_ci_function(
  dataframe = dat_state,
  Region = "S1 NTD",
  predictor_name = "ntd_h_state_3w",
  predictor_group = "H"
)

row_12 <- cor_ci_function(
  dataframe = dat_state,
  Region = "S1 RBD",
  predictor_name = "rbd_h_state_3w",
  predictor_group = "H"
)

# correlation for variant count
row_13 <- cor_ci_function(
  dataframe = dat_state,
  Region = "",
  predictor_name = "n_variants_no_thresh_3w_mean",
  predictor_group = "Freyja variant count"
)

# combine into one table
cor_table <- bind_rows(
  row_1, row_2, row_3, row_4, row_5, row_6,
  row_7, row_8, row_9, row_10, row_11, row_12,
  row_13
)

# make it a flextable and export
# save table as flextable with rounded values
cor_ft <- table_as_flex_function(
  dataframe = cor_table,
  title = "Table 2: Spearman correlations for diversity")

# save
save_as_docx(FitFlextableToPage(cor_ft), 
             path = paste("figures and tables/",
                          "Table S2 - correlation table.docx",
                          sep = "")
)


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
png("figures and tables/Figure S6 - lag plot.png",
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



# ----------------------------------------------------------------
# Table S3 - GLM results for variance explained by S1 NTD and 
# cases/hospitalizations
# ----------------------------------------------------------------


# STATE MODEL 

# remove nas
s <- dat_state %>%
  filter(!is.na(case_incidence)) %>%
  filter(!is.na(ntd_pi_state_3w)) %>%
  filter(!is.na(mean_sars2_conc_state_3w)) %>%
  filter(!is.na(depth_state_3w)) %>%
  filter(!is.na(n_variants_no_thresh_3w_mean))

# model to fit
model_gaussian_pi <- glm(
  case_incidence ~ 
    scale(ntd_pi_state_3w),
  family = gaussian(),
  data = s
)

# create table
state_table_pi <- model_summary_function(
  dataframe = s,
  outcome = "case_incidence",
  model = model_gaussian_pi,
  glmer_option = "No"
) %>%
  mutate(group = "S1 NTD Pi")

# ntd pi and hospital admissions
# model to fit
model_gaussian_pi_2 <- glm(
  hosp_incidence ~ 
    scale(ntd_pi_state_3w),
  family = gaussian(),
  data = s
)

# create table
state_table_pi_2 <- model_summary_function(
  dataframe = s,
  outcome = "hosp_incidence",
  model = model_gaussian_pi_2,
  glmer_option = "No"
) %>%
  mutate(group = "S1 NTD Pi")

# ntd H
model_gaussian_h <- glm(
  case_incidence ~ 
    scale(ntd_h_state_3w),
  family = gaussian(),
  data = s
) 

# table
state_table_h <- model_summary_function(
  dataframe = s,
  outcome = "case_incidence",
  model = model_gaussian_h,
  glmer_option = "No"
)%>%
  mutate(group = "S1 NTD H")

# freyja variant counts
model_gaussian_freyja <- glm(
  case_incidence ~ 
    scale(n_variants_no_thresh_3w_mean),
  family = gaussian(),
  data = s
) 

# table
state_table_freyja <- model_summary_function(
  dataframe = s,
  outcome = "case_incidence",
  model = model_gaussian_freyja,
  glmer_option = "No"
)%>%
  mutate(group = "Freyja variant counts")

# conc only
model_gaussian_conc <- glm(
  case_incidence ~ 
    scale(mean_sars2_conc_state_3w),
  family = gaussian(),
  data = s
)

# table
state_table_conc <- model_summary_function(
  dataframe = s,
  outcome = "case_incidence",
  model = model_gaussian_conc,
  glmer_option = "No"
)%>%
  mutate(group = "Concentration only model")

# CREATE TABLES FOR PAPER

# state model table
table <- bind_rows(state_table_pi,
                   state_table_h,
                   state_table_freyja,
                   state_table_conc) %>%
  dplyr::select(group, everything())
table$`Pr(>|t|)` <- ifelse(table$`Pr(>|t|)` < 0.0001,
                           "<0.0001",
                           table$`Pr(>|t|)`)

# edit estimate column
table$variable<- c(
  "Intercept",
  "S1 NTD Pi",
  "n",
  "Efron R2",
  "AIC",
  "Intercept",
  "S1 NTD H",
  "n",
  "Efron R2",
  "AIC",
  "Intercept",
  "Freyja variant count",
  "n",
  "Efron R2",
  "AIC",
  "Intercept",
  "SARS-CoV-2 concentration",
  "n",
  "Efron R2",
  "AIC"
)


# make it an ft table
t <- table_as_flex_function(dataframe = table,
                            title = "Table S3: Generalized liner model results for S1 NTD")

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
             path = paste("figures and tables/",
                          "Table S3 - state model results.docx",sep = ""))


# -------------------------------------------------------
# Table S4 - county glmm
# -------------------------------------------------------

# county model
# try glmer with gaussian distribution
model_county_pi <- glmer(case_incidence ~ 
                        + scale(ntd_pi_county_3w)
                      +(1|county),
                      family = gaussian()
                      ,
                      data = dat_county)

# table
county_gaus_table_pi <- model_summary_function(
  dataframe = dat_county,
  outcome = "case_incidence",
  model = model_county_pi,
  glmer_option = "Yes"
)%>%
  mutate(group = "S1 NTD Pi")

# H as predictor for county data
model_county_h <- glmer(case_incidence ~ 
                           + scale(ntd_h_county_3w)
                         +(1|county),
                         family = gaussian()
                         ,
                         data = dat_county)

# table
county_gaus_table_h <- model_summary_function(
  dataframe = dat_county,
  outcome = "case_incidence",
  model = model_county_h,
  glmer_option = "Yes"
)%>%
  mutate(group = "S1 NTD H")

# Freyja variants as predictor for county data
model_county_freyja <- glmer(case_incidence ~ 
                           + scale(n_variants_no_thresh_3w_mean)
                         +(1|county),
                         family = gaussian()
                         ,
                         data = dat_county)

# table
county_gaus_table_freyja <- model_summary_function(
  dataframe = dat_county,
  outcome = "case_incidence",
  model = model_county_freyja,
  glmer_option = "Yes"
)%>%
  mutate(group = "Freyja variant counts")

# conc model
model_county_conc <- glmer(case_incidence ~ 
                           + scale(mean_sars2_conc_county_3w)
                         +(1|county),
                         family = gaussian()
                         ,
                         data = dat_county)

# table
county_gaus_table_conc <- model_summary_function(
  dataframe = dat_county,
  outcome = "case_incidence",
  model = model_county_conc,
  glmer_option = "Yes"
)%>%
  mutate(group = "Concentration only model")

# combine into one table
# state model table
table <- bind_rows(county_gaus_table_pi,
                   county_gaus_table_h,
                   county_gaus_table_freyja,
                   county_gaus_table_conc) %>%
  dplyr::select(group, everything())
table$`Pr(>|t|)` <- ifelse(table$`Pr(>|t|)` < 0.0001,
                           "<0.0001",
                           table$`Pr(>|t|)`)

# edit estimate column
table$variable<- c(
  "Intercept",
  "S1 NTD Pi",
  "n",
  "Efron R2",
  "AIC",
  "Intercept",
  "S1 NTD H",
  "n",
  "Efron R2",
  "AIC",
  "Intercept",
  "Freyja variant count",
  "n",
  "Efron R2",
  "AIC",
  "Intercept",
  "SARS-CoV-2 concentration",
  "n",
  "Efron R2",
  "AIC"
)


# make it an ft table
t <- table_as_flex_function(dataframe = table,
                            title = "Table S4: glmm for county data")

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
             path = paste("figures and tables/",
                          "Table S4 - glmm for county data.docx",sep = ""))


# -------------------------------------------------------
# Table S5 - regional glmm
# -------------------------------------------------------
# region model
# try glmer with gaussian distribution
model_region_pi <- glmer(case_incidence ~ 
                           + scale(ntd_pi_region_3w)
                         +(1|region),
                         family = gaussian()
                         ,
                         data = dat_region)

# table
region_gaus_table_pi <- model_summary_function(
  dataframe = dat_region,
  outcome = "case_incidence",
  model = model_region_pi,
  glmer_option = "Yes"
)%>%
  mutate(group = "S1 NTD Pi")

# H as predictor for region data
model_region_h <- glmer(case_incidence ~ 
                          + scale(ntd_h_region_3w)
                        +(1|region),
                        family = gaussian()
                        ,
                        data = dat_region)

# table
region_gaus_table_h <- model_summary_function(
  dataframe = dat_region,
  outcome = "case_incidence",
  model = model_region_h,
  glmer_option = "Yes"
)%>%
  mutate(group = "S1 NTD H")

# Freyja variants as predictor for region data
model_region_freyja <- glmer(case_incidence ~ 
                               + scale(n_variants_no_thresh_3w_mean)
                             +(1|region),
                             family = gaussian()
                             ,
                             data = dat_region)

# table
region_gaus_table_freyja <- model_summary_function(
  dataframe = dat_region,
  outcome = "case_incidence",
  model = model_region_freyja,
  glmer_option = "Yes"
)%>%
  mutate(group = "Freyja variant counts")

# conc model
model_region_conc <- glmer(case_incidence ~ 
                             + scale(mean_sars2_conc_region_3w)
                           +(1|region),
                           family = gaussian()
                           ,
                           data = dat_region)

# table
region_gaus_table_conc <- model_summary_function(
  dataframe = dat_region,
  outcome = "case_incidence",
  model = model_region_conc,
  glmer_option = "Yes"
)%>%
  mutate(group = "Concentration only model")

# combine into one table
# state model table
table <- bind_rows(region_gaus_table_pi,
                   region_gaus_table_h,
                   region_gaus_table_freyja,
                   region_gaus_table_conc) %>%
  dplyr::select(group, everything())
table$`Pr(>|t|)` <- ifelse(table$`Pr(>|t|)` < 0.0001,
                           "<0.0001",
                           table$`Pr(>|t|)`)

# edit estimate column
table$variable<- c(
  "Intercept",
  "S1 NTD Pi",
  "n",
  "Efron R2",
  "AIC",
  "Intercept",
  "S1 NTD H",
  "n",
  "Efron R2",
  "AIC",
  "Intercept",
  "Freyja variant count",
  "n",
  "Efron R2",
  "AIC",
  "Intercept",
  "SARS-CoV-2 concentration",
  "n",
  "Efron R2",
  "AIC"
)


# make it an ft table
t <- table_as_flex_function(dataframe = table,
                            title = "Table S5: glmm for region data")

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
             path = paste("figures and tables/",
                          "Table S5 - glmm for region data.docx",sep = ""))


# -------------------------------------------------------------------------
# Table S6 - Granger causality results
# -------------------------------------------------------------------------

# Granger causality tests for whether time series x predicts time series y
# Both time series are adjusted to be stationary.

# Remove NAs from the data and arrange the time series
s <- dat_state %>%
  filter(!is.na(case_incidence)) %>%
  filter(!is.na(ntd_pi_state_3w)) %>%
  filter(!is.na(mean_sars2_conc_state_3w)) %>%
  filter(!is.na(depth_state_3w)) %>%
  filter(!is.na(n_variants_no_thresh_3w_mean)) %>%
  arrange(week) 

# granger test for s1 ntd pi and cases
g_ntd_pi_cases <- 
  granger_casuality_function(
    dataframe = s,
    y = "case_incidence",
    x = "ntd_pi_state_3w",
    lag = 1,
    group = "S1 NTD Pi",
    outcome = "Case incidence"
  )

# granger test for s1 ntd h and cases
g_ntd_h_cases <- 
  granger_casuality_function(
    dataframe = s,
    y = "case_incidence",
    x = "ntd_h_state_3w",
    lag = 1,
    group = "S1 NTD H",
    outcome = "Case incidence"
  )

# granger test for freyja variant count and  cases
g_var_count_cases <- 
  granger_casuality_function(
    dataframe = s,
    y = "case_incidence",
    x = "n_variants_no_thresh_3w_mean",
    lag = 1,
    group = "Freyja variant count",
    outcome = "Case incidence"
  )

# granger test for s1 ntd pi and hospitalizations
g_ntd_pi_hosp <- 
  granger_casuality_function(
    dataframe = s,
    y = "hosp_incidence",
    x = "ntd_pi_state_3w",
    lag = 1,
    group = "S1 NTD Pi",
    outcome = "Hospitalization incidence"
  )

# granger test for s1 ntd h and hospitalizations
g_ntd_h_hosp <- 
  granger_casuality_function(
    dataframe = s,
    y = "hosp_incidence",
    x = "ntd_h_state_3w",
    lag = 1,
    group = "S1 NTD H",
    outcome = "Hospitalization incidence"
  )

# granger test for freyja variant count and hospitalizations
g_var_count_hosp <- 
  granger_casuality_function(
    dataframe = s,
    y = "hosp_incidence",
    x = "n_variants_no_thresh_3w_mean",
    lag = 1,
    group = "Freyja variant count",
    outcome = "Hospitalization incidence"
  )

# combine into one dataframe
granger_results <- bind_rows(
  g_ntd_pi_cases,
  g_ntd_h_cases,
  g_var_count_cases,
  g_ntd_pi_hosp,
  g_ntd_h_hosp,
  g_var_count_hosp
)

# round results
granger_results <- granger_results %>%
  mutate_if(., is.numeric, round_signifi_function)

# label p values
granger_results <- granger_results %>%
  mutate(`P value` = case_when(
    `P value` > 0.05 ~ as.character(`P value`),
    `P value` < 0.05 & `P value` > 0.01~ paste(as.character(`P value`), "*", 
                                               sep = ""),
    `P value` < 0.01 & `P value` > 0.001 ~ paste(as.character(`P value`), "**", 
                                                 sep = ""),
    `P value` < 0.001 ~ paste("<0.001***")
  ))

# edit column names
colnames(granger_results) <- c("Diversity measure (x)",
                               "Outcome variable (y)",
                               "Model",
                               "F statistic",
                               "P value")

# make it a flextable object
granger_ft <- table_as_flex_function(granger_results,
                                     title = "Table S6: Granger causality results")

granger_ft

# save
save_as_docx(FitFlextableToPage(granger_ft), 
             path = paste("figures and tables/",
                          "Table S6 - granger causality.docx",sep = ""))


# --------------------------------------------------------------------------
# Figure S8 - selective sweep figure showin the virus takeover with the 5% 
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
png("figures and tables/Figure S8 - Jn1 surge.png",
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

# --------------------------------------
# Figure S9 freyja variant counts with 5 % threshold
# --------------------------------------

# JN.1 was first detected in statewide wasteater samples October 29, 2023
# JN.1 was no longer dominant after April 14, 2024

dat_state$jn1 <- ifelse(
  dat_state$week >= "2023-10-29" & dat_state$week <= "2024-04-14",
  "JN.1", "Other Variants"
)

# correlation with freyja lineages
freyja_cor_thresh <- 
  ggplot(data = dat_state,
         aes(x = n_variants_5_3w_mean,
             y = case_incidence)
  )+
  geom_point(color = "grey60", size = 3, alpha = 0.7,
             show.legend = FALSE)+
  stat_cor(method = "spearman"
  )+
  theme_dth_1+
  labs(x = "Mean number of Freyja variants",
       y = "COVID-19 cases per 100k")+
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
            aes(x = week, y= n_variants_5_3w_mean*15,
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
    sec.axis = sec_axis(~ ./ 15, name = "Variant count")
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    label = c("Cases"),
                    name = "")+
  scale_color_manual(values = c("darkblue" = "darkblue"),
                     label = "Freyja lineages",
                     name = "")

var_plot_thresh


# save
png("figures and tables/Figure S9 -  threshold variant count.png",
    units = "in",
    width = 13, height = 6,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(
  plot_grid(
  var_plot_thresh+theme(legend.position = "none"),
  freyja_cor_thresh,
  nrow = 1,
  rel_widths = c(2,1.5),
  align = 'h', axis = 'bt',
  labels = c("A", "B")
  ), 
  g_legend(var_plot_thresh),
  nrow = 2,
  rel_heights = c(2,0.5)
  )
showtext_end()
dev.off()




# -----------------------
# Figure S10 - depth and infections
# -----------------------

p1 <- 
  ggplot(data = dat_state)+
  geom_bar(aes(x = week,
               y = case_incidence,
               fill = "grey60"),
           position = "dodge",
           stat = "identity")+
  geom_line(aes(x = week,
                y = depth_state_3w/10,
                color = "black"),
            lwd = 1)+
  theme_dth_1+
  scale_y_continuous(
    "COVID-19 cases per 100k", 
    sec.axis = sec_axis(~ .*10, name = "Mean depth")
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    name = "",
                    label = c("Cases"))+
  scale_color_manual(values = c("black" = "black"),
                     name = "",
                     label = "Depth")+
  theme(legend.position = "bottom",
        legend.background = element_blank())+
  labs(x = "")

# correlation plot
p2 <- 
  ggplot(data = dat_state,
         aes(x = depth_state_3w,
             y = case_incidence))+
  geom_point(size = 2, color = "black")+
  stat_cor()+
  theme_dth_1+
  labs(x = "Mean depth",
       y = "COVID-19 cases per 100k")



# save
png("figures and tables/Figure S7 -depth and cases over time.png",
    units = "in",
    width = 8.5, height = 6,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(
  plot_grid(p1 + theme(legend.position = "none"), p2, nrow = 1,
            align = "h",
            axis = "bt",
            rel_widths = c(2,1.5),
            labels = c("A", "B")),
  g_legend(p1),
  nrow = 2,
  rel_heights = c(2,0.5))
showtext_end()
dev.off()


# ------------------------------------------
# Figure S11 - depth correlation with pi
# ------------------------------------------

depth_cor <-
  ggplot(data = dat_sewershed,
       aes(x = mean_depth,
           y = ntd_pi))+
  geom_point(
             alpha = 0.5,
             color = "orange")+
  stat_cor(method = "spearman")+
  theme_dth_1+
  labs(x = "Depth of sample",
        y = expression(pi[ww]))

# save
png("figures and tables/Figure S11 -  depth correlation.png",
    units = "in",
    width = 6, height = 6,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
depth_cor
showtext_end()
dev.off()


# ---------------------------------------------------
# Table S7 - variance explained per sample by depth
# ---------------------------------------------------

# fit model predicting pi from depth
model <- lm(genomewide_pi ~ scale(mean_depth),
            data = dat_sewershed)

# model output
table <- 
  model_summary_function(
    dataframe = dat_sewershed,
    outcome = "genomewide_pi",
    model = model,
    glmer_option = "No"
  )

table$`Pr(>|t|)` <- ifelse(table$`Pr(>|t|)` < 0.0001,
                           "<0.0001",
                           table$`Pr(>|t|)`)

# edit estimate column
table$variable<- c(
  "Intercept",
  "Sample mean depth",
  "n",
  "Efron R2",
  "AIC"
)

table <- table %>%
  mutate_if(is.numeric, round_signifi_function)


# make it an ft table
depth_ft <- table_as_flex_function(table,
                                   title = "Table S7: Model for predicting Pi from sample depth")

depth_ft

# save
save_as_docx(FitFlextableToPage(depth_ft), 
             path = paste("figures and tables/",
                          "Table S7 - depth model.docx" ,sep = ""))


# ------------------------------------
# Figure S12 - random subsampling
# ------------------------------------


# Randomly subsample each sample to equal read depth and reassess correlation

# relationship prior to subsample
ggplot(data = dat_sewershed %>%
         filter(region != "New York City"))+
  geom_point(aes(x = week, y = ntd_pi, color = facility_id))+
  theme_dth_1+
  theme(legend.position = "none")+
  geom_smooth(aes(x = week, y = ntd_pi),
              method = "loess",
              span = 0.1,
              color = "darkblue",
              lwd = 1.5)+
  labs(x = "",
       y = "S1 NTD Pi wastewater",
       title = "All samples")+
  scale_color_manual(values = MetBrewer::met.brewer(name = "Austria",
                                                    n = 204))


# subsample for depth reading between 100 and 500
set.seed(23)

# Depth subsample - 3 week moving average for state values
# statewide average - generate example for the figure legend

cases <- dat_state %>%
  select(week, case_incidence)

# sample 3
sample_3_3w <- dat_sewershed  %>%
  filter(depth >= 500 & depth <= 1000) %>%
  sample_n(1000) %>%
  arrange(week) %>%
  group_by(week)%>%
  summarise(
    # statewide pi
    genomewide_pi_state_3w = weighted.mean(x = genomewide_pi_ma3, 
                                           w = population_served, 
                                           na.rm = TRUE),
    ntd_pi_state_3w = weighted.mean(x = ntd_pi_ma3, 
                                    w = population_served, 
                                    na.rm = TRUE)
  )%>%
  ungroup() %>%
  left_join(cases, by = c("week"))

# generate legend
p_time <- 
  ggplot(data = sample_3_3w)+
  geom_bar(data = dat_state %>%
             filter(group == "State"),
           aes(x = week,
               y = case_incidence,
               fill = "grey60"),
           position = "dodge",
           stat = "identity")+
  geom_line(aes(x = week,
                y = ntd_pi_state_3w*10000,
                color = "darkblue"),
            lwd = 1.5)+
  theme_dth_1+
  scale_color_manual(values = c("darkblue"),
                     labels = c("Pi"),
                     name = "")+
  scale_fill_manual(values = c("grey60"),
                    labels = "Cases",
                    name = "")+
  theme(legend.background = element_blank())+
  theme(legend.position = "bottom")

# extract legend function
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# extract legend
mylegend<-g_legend(p_time)



sample_plots_1 <- time_cor_plot_function(
  dataframe = dat_sewershed,
  depth_start = 100,
  depth_end = 400,
  title = "Depth read of sample 100 to 400\n1000 randomly selected samples"
)

sample_plots_2 <- time_cor_plot_function(
  dataframe = dat_sewershed,
  depth_start = 400,
  depth_end = 800,
  title = "Depth read of sample 400 to 800\n1000 randomly selected samples")

sample_plots_3 <- time_cor_plot_function(
  dataframe = dat_sewershed,
  depth_start = 800,
  depth_end = 1500,
  title = "Depth read of sample > 800\n1000 randomly selected samples")

sample_plots_original <- time_cor_plot_function(
  dataframe = dat_sewershed,
  depth_start = 0,
  depth_end = 45000,
  title = "Depth read of sample (all reads)\n1000 randomly selected samples")

ggarrange(ggarrange(sample_plots_original,
                    sample_plots_1,
                    sample_plots_2,
                    sample_plots_3,
                    nrow = 4),mylegend, nrow=2,heights=c(10, 1))


png("figures and tables/Figure S9 -depth subsample moving average.png",
    units = "in",
    width = 8, height = 12,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)

ggarrange(ggarrange(sample_plots_original,
                    sample_plots_1,
                    sample_plots_2,
                    sample_plots_3,
                    nrow = 4),mylegend, nrow=2,heights=c(10, 1))
showtext_end()
dev.off()

# --------------------------------------------------------------------------


# Table S 5 descriptives for conc data
# summarize by pcr lab
conc_df_summary <- dat_sewershed %>%
  ungroup() %>%
  filter(!is.na(pcr_lab)) %>%
  group_by(pcr_lab) %>%
  summarize(`Number of sampling sites` = length(unique(facility_id)),
            `Number of samples` = sum(samples, na.rm  = TRUE),
            `First collection week` = min(week, na.rm = TRUE),
            `Last collection week` = max(week, na.rm= TRUE),
            `Mean concentration` = signif(mean(mean_sars2_conc, na.rm = TRUE),
                                          2),
            `Min concentration` = signif(min(mean_sars2_conc, na.rm = TRUE),
                                         2),
            `Median concentration` = signif(median(mean_sars2_conc, 
                                                   na.rm = TRUE),
                                            2),
            `Max concentration` = signif(max(mean_sars2_conc, na.rm = TRUE), 2),
            `SD concentration` = signif(sd(mean_sars2_conc, na.rm = TRUE),
                                        2)
            
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
  rename(`PCR Lab` = pcr_lab) %>%
  mutate(`Min concentration` = ifelse(
    `Min concentration` == 1, 0, `Min concentration`
  ))

# make it flextable
conc_ft <- table_as_flex_function(
  dataframe = conc_df_summary,
  title = "Table S5: Descriptive statistics for concentration data")

# save
save_as_docx(FitFlextableToPage(conc_ft), 
             path = paste0("figures and tables",
                           "/Table S5 - concentration descriptive statistics.docx")
)

median(depth_1$mean_depth)

# -------------------------------------------------------------------
# Table S7 for case + hosp data
# -------------------------------------------------------------------

hosp_summary <- dat_county %>%
  ungroup() %>%
  filter(week >= "2023-01-01" & week <= "2025-04-20") %>%
  summarize(`Number of weeks with data` = length(unique(week)),
            `Mean` = mean(hospitalizations, na.rm = TRUE),
            `Min` = min(hospitalizations, na.rm = TRUE),
            `Median` = median(hospitalizations, na.rm = TRUE),
            `Max` = max(hospitalizations, na.rm = TRUE),
            `SD` = sd(hospitalizations, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(group = "Hospitalizations") %>%
  mutate_if(is.numeric, ~round(., 2))

# incidence
hosp_incidence_summary <- dat_county %>%
  ungroup() %>%
  filter(week >= "2023-01-01" & week <= "2025-04-20") %>%
  summarize(`Number of weeks with data` = length(unique(week)),
            `Mean` = mean(hosp_incidence, na.rm = TRUE),
            `Min` = min(hosp_incidence, na.rm = TRUE),
            `Median` = median(hosp_incidence, na.rm = TRUE),
            `Max` = max(hosp_incidence, na.rm = TRUE),
            `SD` = sd(hosp_incidence, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(group = "Hospital incidence (per 100k)") %>%
  mutate_if(is.numeric, ~round(., 2))

# cases
cases_summary <- dat_county %>%
  ungroup() %>%
  filter(week >= "2023-01-01" & week <= "2025-04-20") %>%
  summarize(`Number of weeks with data` = length(unique(week)),
            `Mean` = mean(cases, na.rm = TRUE),
            `Min` = min(cases, na.rm = TRUE),
            `Median` = median(cases, na.rm = TRUE),
            `Max` = max(cases, na.rm = TRUE),
            `SD` = sd(cases, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(group = "COVID-19 cases") %>%
  mutate_if(is.numeric, ~round(., 2))

# case incidence
case_incidence_summary <- dat_county %>%
  ungroup() %>%
  filter(week >= "2023-01-01" & week <= "2025-04-20") %>%
  summarize(`Number of weeks with data` = length(unique(week)),
            `Mean` = mean(case_incidence, na.rm = TRUE),
            `Min` = min(case_incidence, na.rm = TRUE),
            `Median` = median(case_incidence, na.rm = TRUE),
            `Max` = max(case_incidence, na.rm = TRUE),
            `SD` = sd(case_incidence, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(group = "COVID-19 case incidence (per 100k)") %>%
  mutate_if(is.numeric, ~round(., 2))

# combine
clinical_df_summary <- bind_rows(
  cases_summary,
  case_incidence_summary,
  hosp_summary,
  hosp_incidence_summary
) %>%
  select(group, everything())

# make it an ft table
# make it flextable
clinical_ft <- table_as_flex_function(
  dataframe = clinical_df_summary,
  title = "Table S7: Descriptive statistics for COVID-19 clinical data")

# save
save_as_docx(FitFlextableToPage(clinical_ft), 
             path = paste0("figures and tables/",
                           "Table S7- clinical descriptive statistics.docx")
)


# ---------------------------------------------
# Figure S14
# ---------------------------------------------
# Depth from freyja variants (alteration detected and depth recorded)
# compared to depth from the total sample
depth_plot_7 <- 
  ggplot(data = dat_sewershed ,
         aes(x = depth,
             y = mean_depth)
  )+
  geom_point(alpha = 0.5, color = "seagreen")+
  stat_cor(method = "spearman",
           digits = 4)+
  theme_dth_1+
  labs(x = "Depth - alterations only",
       y = "Depth - complete sample")

depth_plot_7

# save
png("figures and tables/Figure S14 -depth per sample and depth from alterations.png",
    units = "in",
    width = 5, height = 5,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
depth_plot_7
showtext_end()
dev.off()
