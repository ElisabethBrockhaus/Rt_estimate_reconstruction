library(readr)

setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

# read RKI Nowcast data for RKI estimation
RKI_incid <- load_incidence_data(method = "RKI")

# read smoothed RKI incidence data for SDSC estimation
SDSC_incid <-  read_csv("Rt_estimate_reconstruction/incidence_data/SDSC_incid_21_07_10.csv")

# read incidence data used by rtlive (sourced from RKI line-list data) for other estimations
incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_07_10.csv")


# vary window size
windows <- c(1, 3, 4, 7, 13)
if (exists("estimates_window")) rm(estimates_window)
for (window in windows){
  R_est <- estimate_RKI_R(incid, method = "EpiEstim",
                          window = window,
                          gt_type = "gamma",
                          gt_mean = 4,
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
                    legend_name = "window size", filenames = "_influence_window.pdf",
                    method = "EpiEstim", variation = "with different window sizes",
                    comparing_parameters = TRUE)


# vary generation time distribution
distrs <- c("gamma", "gamma", "constant", "lognorm", "gamma", "empirical", "exponential")
means <-  c(3.4,     3.6,     4.0,        4.7,       4.8,     5.6,         7.0)
sds <-    c(1.8,     3.1,     0.0,        2.9,       2.3,     4.2,         7.0)
gtds <- cbind("type"=distrs, "mean"=means, "sd"=sds)
rownames(gtds) <- c("AGES", "epiforecasts", "RKI", "rtlive", "ETH/SDSC", "Ilmenau", "globalrt")
gtd_strs <- c()

if (exists("estimates_gtd")) rm(estimates_gtd)
for (src in rownames(gtds)){
  gtd_str <- paste0(gtds[src, "type"], ",", gtds[src, "mean"], "(", gtds[src, "sd"], ")")
  print(gtd_str)
  gtd_strs <- c(gtd_strs, gtd_str)
  R_est <- estimate_RKI_R(incid, method = "EpiEstim",
                          window = 7,
                          gt_type = gtds[src, "type"],
                          gt_mean = as.numeric(gtds[src, "mean"]),
                          gt_sd = as.numeric(gtds[src, "sd"]),
                          delay = 0)
  names(R_est) <- c("date", gtd_str, paste0(gtd_str, ".lower"), paste0(gtd_str, ".upper"))
  if (!exists("estimates_gtd")){
    estimates_gtd <- R_est
  } else {
    estimates_gtd <- estimates_gtd %>% full_join(R_est, by = "date")
  }
}
plot_for_comparison(estimates_gtd, comp_methods = gtd_strs,
                    legend_name = "GTD", filename = "influence_GTD.pdf",
                    method = "EpiEstim", variation = "with different generation times",
                    comparing_parameters = FALSE)






