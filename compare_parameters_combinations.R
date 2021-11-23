setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()


################################
# get functions and parameters #
################################

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

# parameter combinations used in papers
gt_dist <- c("gamma", "constant", "empirical", "gamma", "gamma", "gamma", "lognorm", "exponential")
mean_gt <- c(4.8,      4,          5.6,      4.8,    3.4,     3.6,     4.7,       7)
sd_gt <-   c(2.3,      0,          4.2,      2.3,    1.8,     3.1,     2.9,       7)
delay <-   c(11,       1,          7,        10,     0,       12,      12,        0)

params <- data.frame(gtd=gt_dist, gt_mean=mean_gt, gt_sd=sd_gt, delay=delay)
methods <- c("ETH", "RKI", "Ilmenau", "SDSC", "AGES", "epiforecasts", "rtlive", "globalrt")
rownames(params) <- methods

# bring delays in a format usable for ETH estimation
source("Rt_estimate_reconstruction/ETH/delays_for_ETH_estimation.R")


######################
# get incidence data #
######################

# read RKI Nowcast data for RKI estimation
RKI_incid <- load_incidence_data(method = "RKI")

# read smoothed RKI incidence data for SDSC estimation
SDSC_incid <-  read_csv("Rt_estimate_reconstruction/incidence_data/SDSC_incid_21_07_10.csv")

# deconvolve RKI incidence data for ETH estimation
incid_for_ETH <- load_incidence_data(method = "ETHZ_sliding_window", source = "_simpleRKI",
                                     new_deconvolution = FALSE)

# read incidence data used by rtlive (sourced from RKI line-list data) for other estimations
incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_07_10.csv")

# compare incidence time series used by RKI vs. SDSC vs. rtlive
plot(RKI_incid, type="l")
lines(incid, col="blue")
lines(SDSC_incid, col="red")


###############################
# compare real-time estimates #
###############################
RKI_R_pub <- read_csv("https://raw.githubusercontent.com/robert-koch-institut/SARS-CoV-2-Nowcasting_und_-R-Schaetzung/main/Archiv/Nowcast_R_2021-07-10.csv",
                      col_types = list(Datum = col_date())) %>%
  dplyr::select("Datum", "PS_7_Tage_R_Wert", "UG_PI_7_Tage_R_Wert", "OG_PI_7_Tage_R_Wert") %>%
  rename(c("date" = "Datum", "R_pub" = "PS_7_Tage_R_Wert", "lower" = "UG_PI_7_Tage_R_Wert", "upper" = "OG_PI_7_Tage_R_Wert"))
ETH_R_pub <- load_published_R_estimates("ETHZ_sliding_window")
Ilmenau_R_pub <- load_published_R_estimates("ilmenau")
SDSC_R_pub <- load_published_R_estimates("sdsc")
#Zi_R_pub <- load_published_R_estimates("zidatalab")
globalrt_R_pub <- load_published_R_estimates("globalrt_7d")
epiforecasts_R_pub <- load_published_R_estimates("epiforecasts")
rtlive_R_pub <- load_published_R_estimates("rtlive")

# merge estimates and plot for comparison
estimates_pub <- RKI_R_pub[,c("date", "R_pub")] %>%
  full_join(ETH_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(Ilmenau_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(SDSC_R_pub[,c("date", "R_pub")], by = "date") %>% 
  #full_join(Zi_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(globalrt_R_pub[,c("date", "R_pub")], by = "date") %>%
  full_join(epiforecasts_R_pub[,c("date", "R_pub")], by = "date") %>%
  full_join(rtlive_R_pub[,c("date", "R_pub")], by = "date")

org_methods <- c("RKI", "ETH", "Ilmenau", "SDSC", #"Zi",
                 "globalrt", "epiforecasts", "rtlive")
plot_for_comparison(estimates_pub, org_methods,
                    legend_name = "research group", filenames = "_real-time.pdf")

estimates_pub_ci <- RKI_R_pub %>%
  full_join(ETH_R_pub, by = "date") %>%
  full_join(Ilmenau_R_pub, by = "date") %>%
  full_join(SDSC_R_pub, by = "date") %>%
  full_join(globalrt_R_pub, by = "date") %>%
  full_join(epiforecasts_R_pub, by = "date") %>%
  full_join(rtlive_R_pub, by = "date")
org_methods <- c("RKI", "ETH", "Ilmenau", "SDSC", "globalrt", "epiforecasts", "rtlive")

plot_for_comparison(estimates_pub_ci, org_methods, include_CI=T,
                    legend_name = "research group", filenames = "_CI_real-time.pdf")


#####################
# adjust input data #
#####################
R_raw_EpiEstim_adjInput <- estimate_RKI_R(incid, method = "EpiEstim",
                                          window = 7, delay = 0,
                                          gt_type = "gamma", gt_mean = 4, gt_sd = 4)

R_RKI_adjInput <- estimate_RKI_R(RKI_incid)

R_SDSC_EpiEstim_adjInput <- estimate_SDSC_R(SDSC_incid)

R_ETH_EpiEstim_adjInput <- estimate_ETH_R(incid_for_ETH)

R_AGES_EpiEstim_adjInput <- estimate_AGES_R(incid)

R_Ilmenau_adjInput <- estimate_Ilmenau_R(incid)

# TODO: download from server when estimation done
path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
R_epiforecasts_adjInput <- qread(paste0(path, "R_calc_2021-07-10_final_adjInput.qs"))
R_epiforecasts_adjInput <- R_epiforecasts_adjInput[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_adjInput) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
R_globalrt_smoother_adjInput <- read_csv(paste0(path, "bayesian_smoother_7.csv"))
names(R_globalrt_smoother_adjInput) <- c("date", "R_calc", "lower", "upper")
R_globalrt_filter_adjInput <- read_csv(paste0(path, "bayesian_filter_7.csv"))
names(R_globalrt_filter_adjInput) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
R_rtlive_adjInput <- read.csv(paste0(path, "DE_2021-07-10_all_trace_summary.csv")) %>%
  dplyr::select(X, mean, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "mean", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))

# merge estimates and plot for comparison
estimates_adjInput <- R_raw_EpiEstim_adjInput %>%
  full_join(R_RKI_adjInput, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_adjInput, by = "date") %>% 
  full_join(R_ETH_EpiEstim_adjInput, by = "date") %>% 
  full_join(R_AGES_EpiEstim_adjInput, by = "date") %>% 
  full_join(R_Ilmenau_adjInput, by = "date") %>% 
  #full_join(R_epiforecasts_adjInput, by = "date") %>%
  full_join(R_globalrt_smoother_adjInput, by = "date") %>%
  #full_join(R_globalrt_filter_adjInput, by = "date") %>%
  full_join(R_rtlive_adjInput, by = "date")

comp_methods <- c("raw EpiEstim", "RKI", "SDSC", "ETH",
                  "AGES", "Ilmenau", #"epiforecasts",
                  "globalrt smoother", #"globalrt filter",
                  "rtlive")
plot_for_comparison(estimates_adjInput, comp_methods, filenames = "_adjInput.pdf")

estimates_adjInput_CI <- as.data.frame(R_raw_EpiEstim_adjInput) %>% 
  full_join(R_SDSC_EpiEstim_adjInput, by = "date") %>%
  full_join(R_ETH_EpiEstim_adjInput, by = "date") %>%
  full_join(R_AGES_EpiEstim_adjInput, by = "date") %>%
  full_join(R_Ilmenau_adjInput, by = "date") %>% 
  #full_join(R_epiforecasts_adjInput, by = "date") %>%
  full_join(R_globalrt_smoother_adjInput, by = "date") %>%
  #full_join(R_globalrt_filter_adjInput, by = "date") %>%
  full_join(R_rtlive_adjInput, by = "date")

comp_CI <- c("raw EpiEstim", "SDSC", "ETH",
             "AGES", "Ilmenau", #"epiforecasts",
             "globalrt smoother", #"globalrt filter",
             "rtlive")
plot_for_comparison(estimates_adjInput_CI, comp_CI, include_CI = T, filenames = "_CI_adjInput.pdf")


#####################################
# adjust input data and window size #
#####################################
window <- 7

R_raw_EpiEstim_adjInputWindow <- estimate_RKI_R(incid, method = "EpiEstim",
                                          window = window, delay = 0,
                                          gt_type = "gamma", gt_mean = 4, gt_sd = 4)

R_RKI_adjInputWindow <- estimate_RKI_R(RKI_incid, window = window)

R_SDSC_EpiEstim_adjInputWindow <- estimate_SDSC_R(SDSC_incid, window = window)

R_ETH_EpiEstim_adjInputWindow <- estimate_ETH_R(incid_for_ETH, window = window)

R_AGES_EpiEstim_adjInputWindow <- estimate_AGES_R(incid, window = window)

R_Ilmenau_adjInputWindow <- estimate_Ilmenau_R(incid, window = window)

# TODO: download from server when estimation done
R_epiforecasts_adjInputWindow <- R_epiforecasts_adjInput

R_globalrt_smoother_adjInputWindow <- R_globalrt_smoother_adjInput
R_globalrt_filter_adjInputWindow <- R_globalrt_filter_adjInput

R_rtlive_adjInputWindow <- R_rtlive_adjInput

# merge estimates and plot for comparison
estimates_adjInputWindow <- R_raw_EpiEstim_adjInputWindow %>%
  full_join(R_RKI_adjInputWindow, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_adjInputWindow, by = "date") %>% 
  full_join(R_ETH_EpiEstim_adjInputWindow, by = "date") %>% 
  full_join(R_AGES_EpiEstim_adjInputWindow, by = "date") %>% 
  full_join(R_Ilmenau_adjInputWindow, by = "date") %>% 
  #full_join(R_epiforecasts_adjInputWindow, by = "date") %>%
  full_join(R_globalrt_smoother_adjInputWindow, by = "date") %>%
  #full_join(R_globalrt_filter_adjInputWindow, by = "date") %>%
  full_join(R_rtlive_adjInputWindow, by = "date")

comp_methods <- c("raw EpiEstim", "RKI", "SDSC", "ETH",
                  "AGES", "Ilmenau", #"epiforecasts",
                  "globalrt smoother", #"globalrt filter",
                  "rtlive")
plot_for_comparison(estimates_adjInputWindow, comp_methods, filenames = "_adjInputWindow.pdf")

estimates_adjInputWindow_CI <- as.data.frame(R_raw_EpiEstim_adjInputWindow) %>% 
  full_join(R_SDSC_EpiEstim_adjInputWindow, by = "date") %>%
  full_join(R_ETH_EpiEstim_adjInputWindow, by = "date") %>%
  full_join(R_AGES_EpiEstim_adjInputWindow, by = "date") %>%
  full_join(R_Ilmenau_adjInputWindow, by = "date") %>% 
  #full_join(R_epiforecasts_adjInputWindow, by = "date") %>%
  full_join(R_globalrt_smoother_adjInputWindow, by = "date") %>%
  #full_join(R_globalrt_filter_adjInputWindow, by = "date") %>%
  full_join(R_rtlive_adjInputWindow, by = "date")

comp_CI <- c("raw EpiEstim", "SDSC", "ETH",
             "AGES", "Ilmenau", #"epiforecasts",
             "globalrt smoother", #"globalrt filter",
             "rtlive")
plot_for_comparison(estimates_adjInputWindow_CI, comp_CI, include_CI = T, filenames = "_CI_adjInputWindow.pdf")


######################################################
# adjust input data, window size and generation time #
######################################################
window <- 7
gt_type <- "gamma"
gt_mean <- 4
gt_sd <- 4

R_raw_EpiEstim_adjInputWindowGTD <- estimate_RKI_R(incid, method = "EpiEstim",
                                                   window = window, delay = 0,
                                                   gt_type = gt_type, gt_mean = gt_mean, gt_sd = gt_sd)

R_RKI_adjInputWindowGTD <- estimate_RKI_R(RKI_incid,
                                          window = window,
                                          gt_mean = gt_mean) # RKI always assumes constant generation time

R_SDSC_EpiEstim_adjInputWindowGTD <- estimate_SDSC_R(SDSC_incid,
                                                     window = window,
                                                     gt_type = gt_type, gt_mean = gt_mean, gt_sd = gt_sd)

R_ETH_EpiEstim_adjInputWindowGTD <- estimate_ETH_R(incid_for_ETH,
                                                   window = window,
                                                   gt_type = gt_type, gt_mean = gt_mean, gt_sd = gt_sd)

R_AGES_EpiEstim_adjInputWindowGTD <- estimate_AGES_R(incid,
                                                     window = window,
                                                     gt_type = gt_type, gt_mean = gt_mean, gt_sd = gt_sd)

R_Ilmenau_adjInputWindowGTD <- estimate_Ilmenau_R(incid,
                                                  window = window,
                                                  gt_type = gt_type, gt_mean = gt_mean, gt_sd = gt_sd)

# TODO: download from server when estimation done
path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
R_epiforecasts_adjInputWindowGTD <- qread(paste0(path, "R_calc_2021-07-10_final_adjInputWindowGTD.qs"))
R_epiforecasts_adjInputWindowGTD <- R_epiforecasts_adjInputWindowGTD[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_adjInputWindowGTD) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
R_globalrt_smoother_adjInputWindowGTD <- read_csv(paste0(path, "bayesian_smoother_4.csv"))
names(R_globalrt_smoother_adjInputWindowGTD) <- c("date", "R_calc", "lower", "upper")
R_globalrt_filter_adjInputWindowGTD <- read_csv(paste0(path, "bayesian_filter_4.csv"))
names(R_globalrt_filter_adjInputWindowGTD) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
R_rtlive_adjInputWindowGTD <- read.csv(paste0(path, "gamma_4_de_2021-07-10_all_summary.csv")) %>%
  dplyr::select(X, mean, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "mean", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))

# merge estimates and plot for comparison
estimates_adjInputWindowGTD <- R_raw_EpiEstim_adjInputWindowGTD %>%
  full_join(R_RKI_adjInputWindowGTD, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_adjInputWindowGTD, by = "date") %>% 
  full_join(R_ETH_EpiEstim_adjInputWindowGTD, by = "date") %>% 
  full_join(R_AGES_EpiEstim_adjInputWindowGTD, by = "date") %>% 
  full_join(R_Ilmenau_adjInputWindowGTD, by = "date") %>% 
  #full_join(R_epiforecasts_adjInputWindowGTD, by = "date") %>%
  full_join(R_globalrt_smoother_adjInputWindowGTD, by = "date") %>%
  #full_join(R_globalrt_filter_adjInputWindowGTD, by = "date") %>%
  full_join(R_rtlive_adjInputWindowGTD, by = "date")

comp_methods <- c("raw EpiEstim", "RKI", "SDSC", "ETH",
                  "AGES", "Ilmenau", #"epiforecasts",
                  "globalrt smoother", #"globalrt filter",
                  "rtlive")
plot_for_comparison(estimates_adjInputWindowGTD, comp_methods, filenames = "_adjInputWindowGTD.pdf")

estimates_adjInputWindowGTD_CI <- as.data.frame(R_raw_EpiEstim_adjInputWindowGTD) %>% 
  full_join(R_SDSC_EpiEstim_adjInputWindowGTD, by = "date") %>%
  full_join(R_ETH_EpiEstim_adjInputWindowGTD, by = "date") %>%
  full_join(R_AGES_EpiEstim_adjInputWindowGTD, by = "date") %>%
  full_join(R_Ilmenau_adjInputWindowGTD, by = "date") %>% 
  #full_join(R_epiforecasts_adjInputWindowGTD, by = "date") %>%
  full_join(R_globalrt_smoother_adjInputWindowGTD, by = "date") %>%
  #full_join(R_globalrt_filter_adjInputWindowGTD, by = "date") %>%
  full_join(R_rtlive_adjInputWindowGTD, by = "date")

comp_CI <- c("raw EpiEstim", "SDSC", "ETH",
             "AGES", "Ilmenau", #"epiforecasts",
             "globalrt smoother", #"globalrt filter",
             "rtlive")
plot_for_comparison(estimates_adjInputWindowGTD_CI, comp_CI, include_CI = T, filenames = "_CI_adjInputWindowGTD.pdf")


############################################
# adjust everything                        #
# input data, window size, gtd, mean delay #
############################################
window <- 7
gt_type <- "gamma"
gt_mean <- 4
gt_sd <- 4

R_raw_EpiEstim_adjAll <- estimate_RKI_R(incid, method = "EpiEstim",
                                        window = window,
                                        gt_type = gt_type, gt_mean = gt_mean, gt_sd = gt_sd,
                                        delay = 0)

R_RKI_adjAll <- estimate_RKI_R(RKI_incid,
                               window = window,
                               gt_mean = gt_mean, # RKI always assumes constant generation time
                               delay = -3) # account for delay in preprocessing

R_SDSC_EpiEstim_adjAll <- estimate_SDSC_R(SDSC_incid,
                                          window = window,
                                          gt_type = gt_type, gt_mean = gt_mean, gt_sd = gt_sd,
                                          estimateOffsetting = 0) # TODO: check (Where did 10 in pararms come from?)
#R_SDSC_EpiEstim_adjAll$date <- R_SDSC_EpiEstim_adjAll$date + params["SDSC", "delay"]

R_ETH_EpiEstim_adjAll <- estimate_ETH_R(incid_for_ETH,
                                        window = window,
                                        gt_type = gt_type, gt_mean = gt_mean, gt_sd = gt_sd,
                                        shift = params["ETH", "delay"])

R_AGES_EpiEstim_adjAll <- estimate_AGES_R(incid,
                                          window = window,
                                          gt_type = gt_type, gt_mean = gt_mean, gt_sd = gt_sd) # does not consider any delay

R_Ilmenau_adjAll <- estimate_Ilmenau_R(incid,
                                       window = window,
                                       gt_type = gt_type, gt_mean = gt_mean, gt_sd = gt_sd,
                                       delay = 0)

# TODO: download from server when estimation done
path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- "R_calc_2021-07-10_final.qs"
R_epiforecasts_adjAll <- qread(paste0(path, file))
R_epiforecasts_adjAll <- R_epiforecasts_adjAll[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_adjAll) <- c("date", "R_calc", "lower", "upper")

R_globalrt_smoother_adjAll <- R_globalrt_smoother_adjInputWindowGTD
R_globalrt_filter_adjAll <- R_globalrt_filter_adjInputWindowGTD

R_rtlive_adjAll <- R_rtlive_adjInputWindowGTD
R_rtlive_adjAll$date <- R_rtlive_adjAll$date + params["rtlive", "delay"]

# merge estimates and plot for comparison
estimates_adjAll <- R_raw_EpiEstim_adjAll %>%
  full_join(R_RKI_adjAll, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_adjAll, by = "date") %>% 
  full_join(R_ETH_EpiEstim_adjAll, by = "date") %>% 
  full_join(R_AGES_EpiEstim_adjAll, by = "date") %>% 
  full_join(R_Ilmenau_adjAll, by = "date") %>% 
  #full_join(R_epiforecasts_adjAll, by = "date") %>%
  full_join(R_globalrt_smoother_adjAll, by = "date") %>%
  #full_join(R_globalrt_filter_adjAll, by = "date") %>%
  full_join(R_rtlive_adjAll, by = "date")

comp_methods <- c("raw EpiEstim", "RKI", "SDSC", "ETH",
                  "AGES", "Ilmenau", #"epiforecasts",
                  "globalrt smoother", #"globalrt filter",
                  "rtlive")
plot_for_comparison(estimates_adjAll, comp_methods, filenames = "_adjAll.pdf")

estimates_adjAll_CI <- as.data.frame(R_raw_EpiEstim_adjAll) %>% 
  full_join(R_SDSC_EpiEstim_adjAll, by = "date") %>%
  full_join(R_ETH_EpiEstim_adjAll, by = "date") %>%
  full_join(R_AGES_EpiEstim_adjAll, by = "date") %>%
  full_join(R_Ilmenau_adjAll, by = "date") %>% 
  #full_join(R_epiforecasts_adjAll, by = "date") %>%
  full_join(R_globalrt_smoother_adjAll, by = "date") %>%
  #full_join(R_globalrt_filter_adjAll, by = "date") %>%
  full_join(R_rtlive_adjAll, by = "date")

comp_CI <- c("raw EpiEstim", "SDSC", "ETH",
             "AGES", "Ilmenau", #"epiforecasts",
             "globalrt smoother", #"globalrt filter",
             "rtlive")
plot_for_comparison(estimates_adjAll_CI, comp_CI, include_CI = T, filenames = "_CI_adjAll.pdf")


###############################################################################
# variations with single parameter not adjusted

#############################
# adjust everything but gtd #
#############################
window <- 7

R_raw_EpiEstim_notAdjGTD <- estimate_RKI_R(incid, method = "EpiEstim",
                                           window = window,
                                           gt_type = "gamma",
                                           gt_mean = 4,
                                           gt_sd = 4)
R_raw_EpiEstim_notAdjGTD$date <- R_raw_EpiEstim_notAdjGTD$date + params["RKI", "delay"]

R_RKI_notAdjGTD <- estimate_RKI_R(RKI_incid, method = "RKI",
                                  window = window,
                                  gt_type = params["RKI", "gtd"],
                                  gt_mean = params["RKI", "gt_mean"],
                                  gt_sd = params["RKI", "gt_sd"])
R_RKI_notAdjGTD$date <- R_RKI_notAdjGTD$date + params["RKI", "delay"] + 3

R_SDSC_EpiEstim_notAdjGTD <- estimate_SDSC_R(SDSC_incid,
                                             window = window,
                                             gt_type = params["SDSC", "gtd"],
                                             gt_mean = params["SDSC", "gt_mean"],
                                             gt_sd = params["SDSC", "gt_sd"])
R_SDSC_EpiEstim_notAdjGTD$date <- R_SDSC_EpiEstim_notAdjGTD$date + params["SDSC", "delay"]

R_ETH_EpiEstim_notAdjGTD <- estimate_ETH_R(incid_for_ETH,
                                           window = window,
                                           gt_type = params["ETH", "gtd"],
                                           gt_mean = params["ETH", "gt_mean"],
                                           gt_sd = params["ETH", "gt_sd"])
R_ETH_EpiEstim_notAdjGTD$date <- R_ETH_EpiEstim_notAdjGTD$date + params["ETH", "delay"]

R_AGES_EpiEstim_notAdjGTD <- estimate_AGES_R(incid,
                                             window = window,
                                             gt_type = params["AGES", "gtd"],
                                             gt_mean = params["AGES", "gt_mean"],
                                             gt_sd = params["AGES", "gt_sd"],
                                             delay = params["AGES", "delay"])
R_AGES_EpiEstim_notAdjGTD$date <- R_AGES_EpiEstim_notAdjGTD$date + params["AGES", "delay"]

R_Ilmenau_notAdjGTD <- estimate_Ilmenau_R(incid,
                                          window = window,
                                          gt_type = params["Ilmenau", "gtd"],
                                          gt_mean = params["Ilmenau", "gt_mean"],
                                          gt_sd = params["Ilmenau", "gt_sd"],
                                          delay = params["Ilmenau", "delay"])
names(R_Ilmenau_notAdjGTD) <- c("date", "R_calc", "lower", "upper")
R_Ilmenau_notAdjGTD$date <- R_Ilmenau_notAdjGTD$date + params["Ilmenau", "delay"]

# TODO: download from server when estimation done
path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- "R_calc_2021-07-10_final_notAdjGTD.qs"
R_epiforecasts_notAdjGTD <- qread(paste0(path, file))
R_epiforecasts_notAdjGTD <- R_epiforecasts_notAdjGTD[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_notAdjGTD) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
R_globalrt_smoother_notAdjGTD <- read_csv(paste0(path, "bayesian_smoother_7.csv"))
names(R_globalrt_smoother_notAdjGTD) <- c("date", "R_calc", "lower", "upper")
R_globalrt_filter_notAdjGTD <- read_csv(paste0(path, "bayesian_filter_7.csv"))
names(R_globalrt_filter_notAdjGTD) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
file <- "DE_2021-07-10_all_trace_summary.csv"
R_rtlive_notAdjGTD <- read.csv(paste0(path, file)) %>%
  dplyr::select(X, mean, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "mean", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))
R_rtlive_notAdjGTD$date <- R_rtlive_notAdjGTD$date + params["rtlive", "delay"]

# merge estimates and plot for comparison
estimates_notAdjGTD <- R_raw_EpiEstim_notAdjGTD %>%
  full_join(R_RKI_notAdjGTD, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_notAdjGTD, by = "date") %>% 
  full_join(R_ETH_EpiEstim_notAdjGTD, by = "date") %>% 
  full_join(R_AGES_EpiEstim_notAdjGTD, by = "date") %>% 
  full_join(R_Ilmenau_notAdjGTD, by = "date") %>% 
  #full_join(R_epiforecasts_notAdjGTD, by = "date") %>%
  full_join(R_globalrt_smoother_notAdjGTD, by = "date") %>%
  #full_join(R_globalrt_filter_notAdjGTD, by = "date") %>%
  full_join(R_rtlive_notAdjGTD, by = "date")

comp_methods <- c("raw EpiEstim", "RKI", "SDSC", "ETH",
                  "AGES", "Ilmenau", #"epiforecasts",
                  "globalrt smoother", #"globalrt filter",
                  "rtlive")
plot_for_comparison(estimates_notAdjGTD, comp_methods, filenames = "_notAdjGTD.pdf")

estimates_notAdjGTD_CI <- as.data.frame(R_raw_EpiEstim_notAdjGTD) %>% 
  full_join(R_SDSC_EpiEstim_notAdjGTD, by = "date") %>%
  full_join(R_ETH_EpiEstim_notAdjGTD, by = "date") %>%
  full_join(R_AGES_EpiEstim_notAdjGTD, by = "date") %>%
  full_join(R_Ilmenau_notAdjGTD, by = "date") %>% 
  #full_join(R_epiforecasts_notAdjGTD, by = "date") %>%
  full_join(R_globalrt_smoother_notAdjGTD, by = "date") %>%
  #full_join(R_globalrt_filter_notAdjGTD, by = "date") %>%
  full_join(R_rtlive_notAdjGTD, by = "date")

comp_CI <- c("raw EpiEstim", "SDSC", "ETH",
             "AGES", "Ilmenau", #"epiforecasts",
             "globalrt smoother", #"globalrt filter",
             "rtlive")
plot_for_comparison(estimates_notAdjGTD_CI, comp_CI, include_CI = T, filenames = "_CI_notAdjGTD.pdf")




################################
# adjust everything but delays #
################################
window <- 7
gt_type <- "gamma"
gt_mean <- 4
gt_sd <- 4

R_raw_EpiEstim_notAdjDelays <- estimate_RKI_R(incid, method = "EpiEstim",
                                              window = window,
                                              gt_type = gt_type,
                                              gt_mean = gt_mean,
                                              gt_sd = gt_sd)

# there is no gt_type and gt_sd in the original RKI method
R_RKI_notAdjDelays <- estimate_RKI_R(RKI_incid, method = "RKI",
                                     window = window,
                                     gt_type = "constant",
                                     gt_mean = gt_mean,
                                     gt_sd = 0)

R_SDSC_EpiEstim_notAdjDelays <- estimate_SDSC_R(SDSC_incid,
                                                window = window,
                                                gt_type = gt_type,
                                                gt_mean = gt_mean,
                                                gt_sd = gt_sd)

R_ETH_EpiEstim_notAdjDelays <- estimate_ETH_R(incid_for_ETH,
                                              window = window,
                                              gt_type = gt_type,
                                              gt_mean = gt_mean,
                                              gt_sd = gt_sd)

R_AGES_EpiEstim_notAdjDelays <- estimate_AGES_R(incid,
                                                window = window,
                                                gt_type = gt_type,
                                                gt_mean = gt_mean,
                                                gt_sd = gt_sd,
                                                delay = params["AGES", "delay"])

R_Ilmenau_notAdjDelays <- estimate_Ilmenau_R(incid,
                                             window = window,
                                             gt_type = gt_type,
                                             gt_mean = gt_mean,
                                             gt_sd = gt_sd,
                                             delay = params["Ilmenau", "delay"])
names(R_Ilmenau_notAdjDelays) <- c("date", "R_calc", "lower", "upper")

# TODO: download from server when estimation done
path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- "R_calc_2021-07-10_final_notAdjDelay.qs"
R_epiforecasts_notAdjDelays <- qread(paste0(path, file))
R_epiforecasts_notAdjDelays <- R_epiforecasts_notAdjDelays[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_notAdjDelays) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
R_globalrt_smoother_notAdjDelays <- read_csv(paste0(path, "bayesian_smoother_4.csv"))
names(R_globalrt_smoother_notAdjDelays) <- c("date", "R_calc", "lower", "upper")
R_globalrt_filter_notAdjDelays <- read_csv(paste0(path, "bayesian_filter_4.csv"))
names(R_globalrt_filter_notAdjDelays) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
file <- "gamma_4_de_2021-07-10_all_summary.csv"
R_rtlive_notAdjDelays <- read.csv(paste0(path, file)) %>%
  dplyr::select(X, mean, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "mean", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))

# merge estimates and plot for comparison
estimates_notAdjDelays <- R_raw_EpiEstim_notAdjDelays %>%
  full_join(R_RKI_notAdjDelays, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_notAdjDelays, by = "date") %>% 
  full_join(R_ETH_EpiEstim_notAdjDelays, by = "date") %>% 
  full_join(R_AGES_EpiEstim_notAdjDelays, by = "date") %>% 
  full_join(R_Ilmenau_notAdjDelays, by = "date") %>% 
  #full_join(R_epiforecasts_notAdjDelays, by = "date") %>%
  full_join(R_globalrt_smoother_notAdjDelays, by = "date") %>%
  #full_join(R_globalrt_filter_notAdjDelays, by = "date") %>%
  full_join(R_rtlive_notAdjDelays, by = "date")

comp_methods <- c("raw EpiEstim", "RKI", "SDSC", "ETH",
                  "AGES", "Ilmenau", #"epiforecasts",
                  "globalrt smoother", #"globalrt filter",
                  "rtlive")
plot_for_comparison(estimates_notAdjDelays, comp_methods, filenames = "_notAdjDelays.pdf")

estimates_notAdjDelays_CI <- as.data.frame(R_raw_EpiEstim_notAdjDelays) %>% 
  full_join(R_SDSC_EpiEstim_notAdjDelays, by = "date") %>%
  full_join(R_ETH_EpiEstim_notAdjDelays, by = "date") %>%
  full_join(R_AGES_EpiEstim_notAdjDelays, by = "date") %>%
  full_join(R_Ilmenau_notAdjDelays, by = "date") %>% 
  #full_join(R_epiforecasts_notAdjDelays, by = "date") %>%
  full_join(R_globalrt_smoother_notAdjDelays, by = "date") %>%
  #full_join(R_globalrt_filter_notAdjDelays, by = "date") %>%
  full_join(R_rtlive_notAdjDelays, by = "date")

comp_CI <- c("raw EpiEstim", "SDSC", "ETH",
             "AGES", "Ilmenau", #"epiforecasts",
             "globalrt smoother", #"globalrt filter",
             "rtlive")
plot_for_comparison(estimates_notAdjDelays_CI, comp_CI, include_CI = T, filenames = "_CI_notAdjDelays.pdf")



#####################################
# adjust everything but window size #
#####################################
gt_type <- "gamma"
gt_mean <- 4
gt_sd <- 4

R_raw_EpiEstim_notAdjWindow <- estimate_RKI_R(incid, method = "EpiEstim",
                                        gt_type = gt_type,
                                        gt_mean = gt_mean,
                                        gt_sd = gt_sd)
R_raw_EpiEstim_notAdjWindow$date <- R_raw_EpiEstim_notAdjWindow$date + params["RKI", "delay"]

# there is no gt_type and gt_sd in the original RKI method
R_RKI_notAdjWindow <- estimate_RKI_R(RKI_incid, method = "RKI",
                               gt_type = "constant",
                               gt_mean = gt_mean,
                               gt_sd = 0)
R_RKI_notAdjWindow$date <- R_RKI_notAdjWindow$date + params["RKI", "delay"] + 3

R_SDSC_EpiEstim_notAdjWindow <- estimate_SDSC_R(SDSC_incid,
                                          gt_type = gt_type,
                                          gt_mean = gt_mean,
                                          gt_sd = gt_sd)
R_SDSC_EpiEstim_notAdjWindow$date <- R_SDSC_EpiEstim_notAdjWindow$date + params["SDSC", "delay"]

R_ETH_EpiEstim_notAdjWindow <- estimate_ETH_R(incid_for_ETH,
                                        gt_type = gt_type,
                                        gt_mean = gt_mean,
                                        gt_sd = gt_sd)
R_ETH_EpiEstim_notAdjWindow$date <- R_ETH_EpiEstim_notAdjWindow$date + params["ETH", "delay"]

R_AGES_EpiEstim_notAdjWindow <- estimate_AGES_R(incid,
                                          gt_type = gt_type,
                                          gt_mean = gt_mean,
                                          gt_sd = gt_sd,
                                          delay = params["AGES", "delay"])
R_AGES_EpiEstim_notAdjWindow$date <- R_AGES_EpiEstim_notAdjWindow$date + params["AGES", "delay"]

R_Ilmenau_notAdjWindow <- estimate_Ilmenau_R(incid,
                                       gt_type = gt_type,
                                       gt_mean = gt_mean,
                                       gt_sd = gt_sd,
                                       delay = params["Ilmenau", "delay"])
names(R_Ilmenau_notAdjWindow) <- c("date", "R_calc", "lower", "upper")
R_Ilmenau_notAdjWindow$date <- R_Ilmenau_notAdjWindow$date + params["Ilmenau", "delay"]

# TODO: download from server when estimation done
path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- "R_calc_2021-07-10_final.qs"
R_epiforecasts_notAdjWindow <- qread(paste0(path, file))
R_epiforecasts_notAdjWindow <- R_epiforecasts_notAdjWindow[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_notAdjWindow) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
R_globalrt_smoother_notAdjWindow <- read_csv(paste0(path, "bayesian_smoother_4.csv"))
names(R_globalrt_smoother_notAdjWindow) <- c("date", "R_calc", "lower", "upper")
R_globalrt_filter_notAdjWindow <- read_csv(paste0(path, "bayesian_filter_4.csv"))
names(R_globalrt_filter_notAdjWindow) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
file <- "gamma_4_de_2021-07-10_all_summary.csv"
R_rtlive_notAdjWindow <- read.csv(paste0(path, file)) %>%
  dplyr::select(X, mean, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "mean", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))
R_rtlive_notAdjWindow$date <- R_rtlive_notAdjWindow$date + params["rtlive", "delay"]

# merge estimates and plot for comparison
estimates_notAdjWindow <- R_raw_EpiEstim_notAdjWindow %>%
  full_join(R_RKI_notAdjWindow, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_notAdjWindow, by = "date") %>% 
  full_join(R_ETH_EpiEstim_notAdjWindow, by = "date") %>% 
  full_join(R_AGES_EpiEstim_notAdjWindow, by = "date") %>% 
  full_join(R_Ilmenau_notAdjWindow, by = "date") %>% 
  #full_join(R_epiforecasts_notAdjWindow, by = "date") %>%
  full_join(R_globalrt_smoother_notAdjWindow, by = "date") %>%
  #full_join(R_globalrt_filter_notAdjWindow, by = "date") %>%
  full_join(R_rtlive_notAdjWindow, by = "date")

comp_methods <- c("raw EpiEstim", "RKI", "SDSC", "ETH",
                  "AGES", "Ilmenau", #"epiforecasts",
                  "globalrt smoother", #"globalrt filter",
                  "rtlive")
plot_for_comparison(estimates_notAdjWindow, comp_methods, filenames = "_notAdjWindow.pdf",
                    method = method, variation = "not window size")

estimates_notAdjWindow_CI <- as.data.frame(R_raw_EpiEstim_notAdjWindow) %>% 
  full_join(R_SDSC_EpiEstim_notAdjWindow, by = "date") %>%
  full_join(R_ETH_EpiEstim_notAdjWindow, by = "date") %>%
  full_join(R_AGES_EpiEstim_notAdjWindow, by = "date") %>%
  full_join(R_Ilmenau_notAdjWindow, by = "date") %>% 
  #full_join(R_epiforecasts_notAdjWindow, by = "date") %>%
  full_join(R_globalrt_smoother_notAdjWindow, by = "date") %>%
  #full_join(R_globalrt_filter_notAdjWindow, by = "date") %>%
  full_join(R_rtlive_notAdjWindow, by = "date")

comp_CI <- c("raw EpiEstim", "SDSC", "ETH",
             "AGES", "Ilmenau", #"epiforecasts",
             "globalrt smoother", #"globalrt filter",
             "rtlive")
plot_for_comparison(estimates_notAdjWindow_CI, comp_CI, include_CI = T, filenames = "_CI_notAdjWindow.pdf")
