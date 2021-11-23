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



# vary input data
rki_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_11_23.csv")
who_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/epiforecasts_incid_21_11_23.csv")
jhu_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/jhu_incid_21_11_23.csv")
incids <- rki_incid %>%
  inner_join(who_incid, by = "date") %>%
  inner_join(jhu_incid, by = "date") %>%
  dplyr::filter(date<=as_date("2021-07-10"))
data_sources <- c("RKI", "WHO", "JHU")
colnames(incids) <- c("date", data_sources)

if (exists("estimates_input")) rm(estimates_input)
for (src in data_sources){
  incid <- incids[c("date", src)] %>% rename(c("I"=src))
  R_est <- estimate_RKI_R(incid, method = "EpiEstim",
                          window = 7,
                          gt_type = "gamma",
                          gt_mean = 4,
                          gt_sd = 4,
                          delay = 0)
  names(R_est) <- c("date", src, paste0(src, ".lower"), paste0(src, ".upper"))
  if (!exists("estimates_input")){
    estimates_input <- R_est
  } else {
    estimates_input <- estimates_input %>% full_join(R_est, by = "date")
  }
}
plot_for_comparison(estimates_input, comp_methods = data_sources,
                    legend_name = "data source", filenames = "_influence_input_data.pdf",
                    method = "EpiEstim", variation = "with different sources of input data",
                    comparing_parameters = FALSE)


# vary preprocessing
rki_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_07_10.csv")
rki_nowcast <- load_RKI_data()
sdsc_smoothed <- load_SDSC_data(country="Germany", data_status="2021-07-10")
sdsc_smoothed$I <- round(sdsc_smoothed$I)
incids <- rki_incid %>%
  inner_join(rki_nowcast, by = "date") %>%
  inner_join(sdsc_smoothed, by = "date")
preprocessing <- c("none", "nowcast RKI", "smoothing SDSC")
colnames(incids) <- c("date", preprocessing)

if (exists("estimates_preprocess")) rm(estimates_preprocess)
for (type in preprocessing){
  incid <- incids[c("date", type)] %>% rename(c("I"=type))
  R_est <- estimate_RKI_R(incid, method = "EpiEstim",
                          window = 7,
                          gt_type = "gamma",
                          gt_mean = 4,
                          gt_sd = 4,
                          delay = 0)
  names(R_est) <- c("date", type, paste0(type, ".lower"), paste0(type, ".upper"))
  if (!exists("estimates_preprocess")){
    estimates_preprocess <- R_est
  } else {
    estimates_preprocess <- estimates_preprocess %>% full_join(R_est, by = "date")
  }
}

incid_for_ETH <- load_incidence_data(method = "ETHZ_sliding_window", source = "_simpleRKI",
                                     new_deconvolution = FALSE)
R_ETH <- estimate_ETH_R(incid_for_ETH,
                        window = 7,
                        gt_type = "gamma",
                        gt_mean = 4,
                        gt_sd = 4)
R_ETH$date <- R_ETH$date + 11 # mean delay of ETH estimates
estimates_preprocess <- estimates_preprocess %>% full_join(R_ETH, by = "date")

plot_for_comparison(estimates_preprocess, comp_methods = c(preprocessing, "deconvolution ETH"),
                    legend_name = "preprocessing", filenames = "_influence_preprocessing.pdf",
                    method = "EpiEstim", variation = "with different preprocessing methods")

