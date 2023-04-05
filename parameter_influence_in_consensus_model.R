library(lubridate)
library(dplyr)

setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction"
getwd()

# path for saving estimates
path_estimates <- "Rt_estimate_reconstruction/estimates/"


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
methods_window <- c("Ilmenau", "ETH", "SDSC", "RKI")
windows <-         c(1,         3,     4,      7)
names(windows) <- methods_window
window_strs <- c()

if (exists("estimates_window")) rm(estimates_window)
for (src in names(windows)){
  window_str <- paste0(windows[src], ", ", src)
  print(window_str)
  window_strs <- c(window_strs, window_str)
  R_est <- estimate_RKI_R(RKI_incid, method = "EpiEstim",
                          window = windows[src],
                          gt_type = "gamma",
                          gt_mean = 4,
                          gt_sd = 4,
                          delay = 0)
  names(R_est) <- c("date", window_str, paste0(window_str, ".lower"), paste0(window_str, ".upper"))
  if (!exists("estimates_window")){
    estimates_window <- R_est
  } else {
    estimates_window <- estimates_window %>% full_join(R_est, by = "date")
  }
}
# find ylim
ylim_window_l <- min(colMin(estimates_window %>%
                              dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>%
                              dplyr::select(ends_with(as.character(window_strs)))))
ylim_window_u <- max(colMax(estimates_window %>%
                              dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>%
                              dplyr::select(ends_with(as.character(window_strs)))))
# plot
plot_for_comparison(estimates_window, comp_methods = window_strs,
                    col_palette = "YlOrRd", name_consensus = "7, RKI",
                    legend_name = "window size", filenames = "_influence_window.pdf",
                    sort_numerically = TRUE,  plot_diff_matrices=T,
                    ylim_l = ylim_window_l, ylim_u = ylim_window_u)


#####################################
# vary generation time distribution #
#####################################
methods_gtd <- c("epiforecasts", "RKI",      "consensus", "rtlive",  "ETH/SDSC", "Ilmenau", "globalrt", "HZI")
distrs <-      c("gamma",        "constant", "gamma",     "lognorm", "gamma",    "ad hoc",  "gamma",    "?")
means <-       c( 3.6,            4.0,        4.0,         4.7,       4.8,        5.6,       7.0,        10.3)
sds <-         c( 3.1,            0.0,        4.0,         2.9,       2.3,        4.2,       7.0,        7.6)
gtds <- cbind("type"=distrs, "mean"=means, "sd"=sds)
rownames(gtds) <- methods_gtd
gtd_strs <- c()

if (exists("estimates_gtd")) rm(estimates_gtd)
for (src in rownames(gtds)){
  gtd_str <- paste0(gtds[src, "mean"], "(", gtds[src, "sd"], "), ", src)
  print(gtd_str)
  gtd_strs <- c(gtd_strs, gtd_str)
  R_est <- estimate_RKI_R(RKI_incid, method = "EpiEstim",
                          window = 7,
                          gt_type = ifelse(src == "RKI", "constant", "gamma"),
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
ylim_gtd_l <- min(colMin(estimates_gtd %>%
                           dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>%
                           dplyr::select(ends_with(gtd_strs))))
ylim_gtd_u <- max(colMax(estimates_gtd %>%
                           dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>%
                           dplyr::select(ends_with(gtd_strs))))
# plot
plot_for_comparison(estimates_gtd, comp_methods = gtd_strs,
                    col_palette = "YlGn", name_consensus = "4(4), consensus",
                    legend_name = "GTD", filenames = "_influence_GTD.pdf",
                    sort_numerically = TRUE, plot_diff_matrices=T,
                    ylim_l = ylim_gtd_l, ylim_u = ylim_gtd_u)


###################
# vary input data #
###################
rki_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_11_23.csv")
who_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/epiforecasts_incid_21_11_23.csv")
jhu_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/jhu_incid_21_11_23.csv")
rki_incid_nowcast <- read_csv("https://raw.githubusercontent.com/robert-koch-institut/SARS-CoV-2-Nowcasting_und_-R-Schaetzung/main/Archiv/Nowcast_R_2021-11-23.csv") %>%
  dplyr::select(Datum, PS_COVID_Faelle) %>%
  rename(date = Datum)
incids <- rki_incid %>%
  inner_join(who_incid, by = "date") %>%
  inner_join(jhu_incid, by = "date") %>%
  inner_join(rki_incid_nowcast, by = "date") %>%
  dplyr::filter(date<=as_date("2021-07-10"))
data_sources <- c("RKI, positive test", "WHO", "JHU", "RKI, symptom onset")
colnames(incids) <- c("date", data_sources)

if (exists("estimates_input")) rm(estimates_input)
for (data_src in data_sources){
  incid <- incids[c("date", all_of(data_src))] %>% rename(c("I"=data_src))
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
ylim_input_l <- min(colMin(estimates_input %>% 
                             dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>%
                             dplyr::select(ends_with(data_sources))))
ylim_input_u <- max(colMax(estimates_input %>%
                             dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>%
                             dplyr::select(ends_with(data_sources))))
# plot
plot_for_comparison(estimates_input, comp_methods = data_sources,
                    col_palette = "incidence data", name_consensus = "RKI, positive test",
                    legend_name = "data source", filenames = "_influence_input_data.pdf",
                    sort_numerically = FALSE, plot_diff_matrices=T,
                    ylim_l = ylim_input_l, ylim_u = ylim_input_u)


######################
# vary preprocessing #
######################
incids <- RKI_incid %>%
  inner_join(RKI_nowcast, by = "date") %>%
  inner_join(SDSC_smoothed, by = "date")
preprocessing <- c("none", "RKI", "SDSC")
colnames(incids) <- c("date", preprocessing)

if (exists("estimates_preprocess")) rm(estimates_preprocess)
for (type in preprocessing){
  incid <- incids[c("date", all_of(type))] %>% rename(c("I"=type))
  R_est <- estimate_RKI_R(incid, method = "EpiEstim",
                          window = 7,
                          gt_type = "gamma",
                          gt_mean = 4,
                          gt_sd = 4,
                          delay = 0)
  names(R_est) <- c("date", type, paste0(type, ".lower"), paste0(type, ".upper"))
  
  if (type == "RKI"){
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
ylim_preprocess_l <- min(colMin(estimates_preprocess %>%
                                  dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>%
                                  dplyr::select(ends_with(c(preprocessing, "R_calc")))))
ylim_preprocess_u <- max(colMax(estimates_preprocess %>%
                                  dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>%
                                  dplyr::select(ends_with(c(preprocessing, "R_calc")))))
# plot
plot_for_comparison(estimates_preprocess, comp_methods = c(preprocessing, "ETH"),
                    col_palette = "Dark2", name_consensus = "none",
                    legend_name = "preprocessing", filenames = "_influence_preprocessing.pdf",
                    sort_numerically = FALSE, plot_diff_matrices=T,
                    ylim_l = ylim_preprocess_l, ylim_u = ylim_preprocess_u)


###############################################################################
# APPENDIX

######################################
# vary standard deviation of the GTD #
######################################
sds <- c(3.1, 0.001, 4.0, 2.9, 2.3, 4.2, 7.0, 7.6)

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
#ylim_SD_gtd_l <- min(colMin(estimates_SD_gtd %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(ends_with(as.character(sds)))))
#ylim_SD_gtd_u <- max(colMax(estimates_SD_gtd %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(ends_with(as.character(sds)))))
# plot
plot_for_comparison(estimates_SD_gtd, comp_methods = as.character(sds),
                    col_palette = "YlGn", name_consensus = 4.0,
                    legend_name = "SD of the GTD", filenames = "_influence_SD_GTD.pdf",
                    sort_numerically = TRUE, plot_diff_matrices=T,
                    ylim_l = ylim_gtd_l, ylim_u = ylim_gtd_u)


##################
# save estimates #
##################
write_csv(estimates_window, paste0(path_estimates, "R_cmp_Window_2021-07-10.csv"))
write_csv(estimates_gtd, paste0(path_estimates, "R_cmp_GTD_2021-07-10.csv"))
write_csv(estimates_input, paste0(path_estimates, "R_cmp_Input_2021-11-23.csv"))
write_csv(estimates_preprocess, paste0(path_estimates, "R_cmp_Preprocess_2021-07-10.csv"))
write_csv(estimates_SD_gtd, paste0(path_estimates, "R_cmp_GTD_SD_2021-07-10.csv"))
