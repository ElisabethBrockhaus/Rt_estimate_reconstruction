library(readr)

setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid.csv")

# vary window size
windows <- c(1, 3, 7, 13)
if (exists("estimates_window")) rm(estimates_window)
for (window in windows){
  R_est <- estimate_RKI_R(incid, method = "EpiEstim",
                          window = window,
                          gt_type = "gamma",
                          gt_mean = 5,
                          gt_sd = 4,
                          delay = 0)
  names(R_est) <- c("date", window, paste0(window, ".lower"), paste0(window, ".upper"))
  if (!exists("estimates_window")){
    estimates_window <- R_est
  } else {
    estimates_window <- estimates_window %>% full_join(R_est, by = "date")
  }
}
plot_for_comparison(estimates_window, comp_methods = windows,
                    legend_name = "window size", filenames = "influence_window_m.pdf",
                    method = "EpiEstim", variation = "with different window sizes",
                    comparing_parameters = TRUE)

windows <- 1:7
if (exists("estimates_window")) rm(estimates_window)
for (window in windows){
  R_est <- estimate_RKI_R(incid, method = "EpiEstim",
                          window = window,
                          gt_type = "gamma",
                          gt_mean = 5,
                          gt_sd = 4,
                          delay = 0)
  names(R_est) <- c("date", window, paste0(window, ".lower"), paste0(window, ".upper"))
  if (!exists("estimates_window")){
    estimates_window <- R_est
  } else {
    estimates_window <- estimates_window %>% full_join(R_est, by = "date")
  }
}
plot_for_comparison(estimates_window, comp_methods = windows,
                    legend_name = "window size", filename = "influence_window_c.png",
                    method = "EpiEstim", variation = "with different window sizes",
                    comparing_parameters = TRUE)

# vary mean generation time
means <- c(3.4, 3.6, 4.0, 4.7, 4.8, 5.0, 5.6, 7.0)
if (exists("estimates_mean")) rm(estimates_mean)
for (mean in means){
  R_est <- estimate_RKI_R(incid, method = "EpiEstim",
                          window = 7,
                          gt_type = "gamma",
                          gt_mean = mean,
                          gt_sd = 4,
                          delay = 0)
  names(R_est) <- c("date", mean, paste0(mean, ".lower"), paste0(mean, ".upper"))
  if (!exists("estimates_mean")){
    estimates_mean <- R_est
  } else {
    estimates_mean <- estimates_mean %>% full_join(R_est, by = "date")
  }
}
plot_for_comparison(estimates_mean, comp_methods = means,
                    legend_name = "mean GT", filename = "influence_meanGT_m.png",
                    method = "EpiEstim", variation = "with different mean generation times",
                    comparing_parameters = TRUE)

means <- 1:10
if (exists("estimates_mean")) rm(estimates_mean)
for (mean in means){
  R_est <- estimate_RKI_R(incid, method = "EpiEstim",
                          window = 7,
                          gt_type = "gamma",
                          gt_mean = mean,
                          gt_sd = 4,
                          delay = 0)
  names(R_est) <- c("date", mean, paste0(mean, ".lower"), paste0(mean, ".upper"))
  if (!exists("estimates_mean")){
    estimates_mean <- R_est
  } else {
    estimates_mean <- estimates_mean %>% full_join(R_est, by = "date")
  }
}
plot_for_comparison(estimates_mean, comp_methods = means,
                    legend_name = "mean GT", filename = "influence_meanGT_c.png",
                    method = "EpiEstim", variation = "with different mean generation times",
                    comparing_parameters = TRUE)

# vary sd generation time
sds <- c(0.001, 1.8, 2.3, 2.9, 3.1, 4.0, 4.2, 7.0)
if (exists("estimates_sd")) rm(estimates_sd)
for (sd in sds){
  R_est <- estimate_RKI_R(incid, method = "EpiEstim",
                          window = 7,
                          gt_type = "gamma",
                          gt_mean = 5,
                          gt_sd = sd,
                          delay = 0)
  names(R_est) <- c("date", sd, paste0(sd, ".lower"), paste0(sd, ".upper"))
  if (!exists("estimates_sd")){
    estimates_sd <- R_est
  } else {
    estimates_sd <- estimates_sd %>% full_join(R_est, by = "date")
  }
}
plot_for_comparison(estimates_sd, comp_methods = sds,
                    legend_name = "sd GT", filename = "influence_sdGT_m.png",
                    method = "EpiEstim", variation = "with different sd of generation time",
                    comparing_parameters = TRUE)

sds <- 1:10
if (exists("estimates_sd")) rm(estimates_sd)
for (sd in sds){
  R_est <- estimate_RKI_R(incid, method = "EpiEstim",
                          window = 7,
                          gt_type = "gamma",
                          gt_mean = 5,
                          gt_sd = sd,
                          delay = 0)
  names(R_est) <- c("date", sd, paste0(sd, ".lower"), paste0(sd, ".upper"))
  if (!exists("estimates_sd")){
    estimates_sd <- R_est
  } else {
    estimates_sd <- estimates_sd %>% full_join(R_est, by = "date")
  }
}
plot_for_comparison(estimates_sd, comp_methods = sds,
                    legend_name = "sd GT", filename = "influence_sdGT_c.png",
                    method = "EpiEstim", variation = "with different sd of generation time",
                    comparing_parameters = TRUE)



