setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

# read RKI Nowcast data for RKI estimation
RKI_nowcast <- read_csv("Rt_estimate_reconstruction/incidence_data/RKI_nowcast_21_07_10.csv")

# read smoothed RKI incidence data for SDSC estimation
SDSC_smoothed <-  read_csv("Rt_estimate_reconstruction/incidence_data/SDSC_incid_21_07_10.csv")

# read deconvolved RKI incidence data for ETH estimation
ETH_incid_deconvoluted <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_for_ETH_21_07_10.csv")

# read incidence data used by rtlive (sourced from RKI line-list data) for other estimations
RKI_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_07_10.csv")


# function to calculate min/max over all columns of a data frame
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMin <- function(data) sapply(data, min, na.rm = TRUE)


####################
# vary window size #
####################
windows <- c(1, 3, 4, 7, 13)
if (exists("estimates_window")) rm(estimates_window)
for (window in windows){
  R_est <- estimate_RKI_R(RKI_incid, method = "EpiEstim",
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
# find ylim
max <- max(colMax(estimates_window %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(ends_with(as.character(windows)))))
min <- min(colMin(estimates_window %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(ends_with(as.character(windows)))))
# plot
plot_for_comparison(estimates_window, comp_methods = windows,
                    legend_name = "window size", filenames = "_influence_window.pdf",
                    sort_numerically = TRUE, ylims = c(min, max))


#####################################
# vary generation time distribution #
#####################################
distrs <- c("gamma", "gamma", "constant", "gamma", "lognorm", "gamma", "empirical", "exponential")
means <-  c(3.4,     3.6,     4.0,        4.0,     4.7,       4.8,     5.6,         7.0)
sds <-    c(1.8,     3.1,     0.0,        4.0,     2.9,       2.3,     4.2,         7.0)
gtds <- cbind("type"=distrs, "mean"=means, "sd"=sds)
rownames(gtds) <- c("AGES", "epiforecasts", "RKI", "consensus", "rtlive", "ETH/SDSC", "Ilmenau", "globalrt")
gtd_strs <- c()

if (exists("estimates_gtd")) rm(estimates_gtd)
for (src in rownames(gtds)){
  gtd_str <- paste0(gtds[src, "mean"], "(", gtds[src, "sd"], "), ", gtds[src, "type"])
  print(gtd_str)
  gtd_strs <- c(gtd_strs, gtd_str)
  R_est <- estimate_RKI_R(RKI_incid, method = "EpiEstim",
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
# find ylim
max <- max(colMax(estimates_gtd %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(ends_with(gtd_strs))))
min <- min(colMin(estimates_gtd %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(ends_with(gtd_strs))))
# plot
plot_for_comparison(estimates_gtd, comp_methods = gtd_strs,
                    legend_name = "GTD", filenames = "_influence_GTD.pdf",
                    sort_numerically = FALSE, ylims = c(min, max))


###################
# vary input data #
###################
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
for (data_src in data_sources){
  incid <- incids[c("date", data_src)] %>% rename(c("I"=data_src))
  R_est <- estimate_RKI_R(incid, method = "EpiEstim",
                          window = 7,
                          gt_type = "gamma",
                          gt_mean = 4,
                          gt_sd = 4,
                          delay = 0)
  names(R_est) <- c("date", data_src, paste0(data_src, ".lower"), paste0(data_src, ".upper"))
  if (!exists("estimates_input")){
    estimates_input <- R_est
  } else {
    estimates_input <- estimates_input %>% full_join(R_est, by = "date")
  }
}
# find ylim
max <- max(colMax(estimates_input %>% dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>% dplyr::select(ends_with(data_sources))))
min <- min(colMin(estimates_input %>% dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>% dplyr::select(ends_with(data_sources))))
# plot
plot_for_comparison(estimates_input, comp_methods = data_sources,
                    legend_name = "data source", filenames = "_influence_input_data.pdf",
                    sort_numerically = FALSE, ylims = c(min, max))


######################
# vary preprocessing #
######################
incids <- RKI_incid %>%
  inner_join(RKI_nowcast, by = "date") %>%
  inner_join(SDSC_smoothed, by = "date")
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
  
  if (type == "nowcast RKI"){
    R_est$date <- R_est$date + 3 # mean onset-to-reporting delay in RKI line list data
  }
  
  if (!exists("estimates_preprocess")){
    estimates_preprocess <- R_est
  } else {
    estimates_preprocess <- estimates_preprocess %>% full_join(R_est, by = "date")
  }
}

R_ETH <- estimate_ETH_R(ETH_incid_deconvoluted,
                        window = 7,
                        gt_type = "gamma",
                        gt_mean = 4,
                        gt_sd = 4)
R_ETH$date <- R_ETH$date + 11 # mean delay of ETH estimates
estimates_preprocess <- estimates_preprocess %>% full_join(R_ETH, by = "date")

# find ylim
max <- max(colMax(estimates_preprocess %>% dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>% dplyr::select(ends_with(c(preprocessing, "deconvolution ETH")))))
min <- min(colMin(estimates_preprocess %>% dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>% dplyr::select(ends_with(c(preprocessing, "deconvolution ETH")))))
# plot
plot_for_comparison(estimates_preprocess, comp_methods = c(preprocessing, "deconvolution ETH"),
                    legend_name = "preprocessing", filenames = "_influence_preprocessing.pdf",
                    sort_numerically = FALSE, ylims = c(min, max))


###############################################################################
# APPENDIX

######################################
# vary standard deviation of the GTD #
######################################
sds <-    c(1.8, 3.1, 0.001, 4.0, 2.9, 2.3, 4.2, 7.0)

if (exists("estimates_SD_gtd")) rm(estimates_SD_gtd)
for (sd in sds){
  R_est <- estimate_RKI_R(RKI_incid, method = "EpiEstim",
                          window = 7,
                          gt_type = "gamma",
                          gt_mean = 4,
                          gt_sd = sd,
                          delay = 0)
  names(R_est) <- c("date", sd, paste0(sd, ".lower"), paste0(sd, ".upper"))
  if (!exists("estimates_SD_gtd")){
    estimates_SD_gtd <- R_est
  } else {
    estimates_SD_gtd <- estimates_SD_gtd %>% full_join(R_est, by = "date")
  }
}
# find ylim
max <- max(colMax(estimates_SD_gtd %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(ends_with(as.character(sds)))))
min <- min(colMin(estimates_SD_gtd %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(ends_with(as.character(sds)))))
# plot
plot_for_comparison(estimates_SD_gtd, comp_methods = sds,
                    legend_name = "SD of the GTD", filenames = "_influence_SD_GTD.pdf",
                    sort_numerically = TRUE, ylims = c(min, max))
