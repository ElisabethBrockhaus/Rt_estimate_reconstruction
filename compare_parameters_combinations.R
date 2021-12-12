setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()


# TODO: save estimates after calculations


################################
# get functions and parameters #
################################

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

# parameter combinations used in papers
gt_dist <- c("gamma", "constant", "ad hoc", "gamma", "gamma", "gamma", "lognorm", "exponential")
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
#RKI_incid <- load_incidence_data(method = "RKI")
#write_csv(RKI_incid, "Rt_estimate_reconstruction/incidence_data/RKI_nowcast_21_07_10.csv")
RKI_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/RKI_nowcast_21_07_10.csv")

# read smoothed RKI incidence data for SDSC estimation
SDSC_incid <-  read_csv("Rt_estimate_reconstruction/incidence_data/SDSC_incid_21_07_10.csv")

# deconvolve RKI incidence data for ETH estimation
#incid_for_ETH <- load_incidence_data(method = "ETHZ_sliding_window", source = "_RKI_2021_07_10",
#                                     new_deconvolution = TRUE)
#write_csv(incid_for_ETH, "Rt_estimate_reconstruction/incidence_data/rtlive_incid_for_ETH_21_07_10.csv")
incid_for_ETH <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_for_ETH_21_07_10.csv")

# read incidence data used by rtlive (sourced from RKI line-list data) for other estimations
incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_07_10.csv")


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
#plot_for_comparison(estimates_pub, org_methods,
#                    legend_name = "research group", filenames = "_real-time.pdf")

estimates_pub_ci <- RKI_R_pub %>%
  full_join(ETH_R_pub, by = "date") %>%
  full_join(Ilmenau_R_pub, by = "date") %>%
  full_join(SDSC_R_pub, by = "date") %>%
  full_join(globalrt_R_pub, by = "date") %>%
  full_join(epiforecasts_R_pub, by = "date") %>%
  full_join(rtlive_R_pub, by = "date")

#plot_for_comparison(estimates_pub_ci, org_methods, include_CI=T,
#                    legend_name = "research group", filenames = "_CI_real-time.pdf")


#####################
# adjust input data #
#####################
R_consensus_adjInput <- estimate_RKI_R(incid, method = "EpiEstim",
                                          window = 7, delay = 0,
                                          gt_type = "gamma", gt_mean = 4, gt_sd = 4)

R_RKI_adjInput <- estimate_RKI_R(RKI_incid)

R_SDSC_EpiEstim_adjInput <- estimate_SDSC_R(SDSC_incid)

R_ETH_EpiEstim_adjInput <- estimate_ETH_R(incid_for_ETH)

R_AGES_EpiEstim_adjInput <- estimate_AGES_R(incid)

R_Ilmenau_adjInput <- estimate_Ilmenau_R(incid)

path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
R_epiforecasts_adjInput <- qread(paste0(path, "R_calc_2021-07-10_final_adjInput.qs"))
R_epiforecasts_adjInput <- R_epiforecasts_adjInput[,c("date", "median", "lower_95", "upper_95")]
names(R_epiforecasts_adjInput) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
R_globalrt_smoother_adjInput <- read_csv(paste0(path, "bayesian_smoother_7.csv"))[,c("Date", "R", "lb_95", "ub_95")]
names(R_globalrt_smoother_adjInput) <- c("date", "R_calc", "lower", "upper")
#R_globalrt_filter_adjInput <- read_csv(paste0(path, "bayesian_filter_7.csv"))
#names(R_globalrt_filter_adjInput) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
R_rtlive_adjInput <- read.csv(paste0(path, "DE_2021-07-10_all_trace_summary.csv")) %>%
  dplyr::select(X, median, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "median", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))

# merge estimates and plot for comparison
estimates_adjInput <- R_consensus_adjInput %>%
  full_join(R_RKI_adjInput, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_adjInput, by = "date") %>% 
  full_join(R_ETH_EpiEstim_adjInput, by = "date") %>% 
  full_join(R_AGES_EpiEstim_adjInput, by = "date") %>% 
  full_join(R_Ilmenau_adjInput, by = "date") %>% 
  full_join(R_epiforecasts_adjInput, by = "date") %>%
  full_join(R_globalrt_smoother_adjInput, by = "date") %>%
  ##full_join(R_globalrt_filter_adjInput, by = "date") %>%
  full_join(R_rtlive_adjInput, by = "date") %>%
  arrange(date)

comp_methods <- c("consensus", "RKI", "SDSC", "ETH",
                  "AGES", "Ilmenau", "epiforecasts",
                  "globalrt", #"globalrt filter",
                  "rtlive")
#plot_for_comparison(estimates_adjInput, comp_methods, filenames = "_adjInput.pdf")

estimates_adjInput_CI <- as.data.frame(R_consensus_adjInput) %>% 
  full_join(R_SDSC_EpiEstim_adjInput, by = "date") %>%
  full_join(R_ETH_EpiEstim_adjInput, by = "date") %>%
  full_join(R_AGES_EpiEstim_adjInput, by = "date") %>%
  full_join(R_Ilmenau_adjInput, by = "date") %>% 
  full_join(R_epiforecasts_adjInput, by = "date") %>%
  full_join(R_globalrt_smoother_adjInput, by = "date") %>%
  #full_join(R_globalrt_filter_adjInput, by = "date") %>%
  full_join(R_rtlive_adjInput, by = "date") %>%
  arrange(date)

comp_CI <- c("consensus", "SDSC", "ETH",
             "AGES", "Ilmenau", "epiforecasts",
             "globalrt", #"globalrt filter",
             "rtlive")
#plot_for_comparison(estimates_adjInput_CI, comp_CI, include_CI = T, filenames = "_CI_adjInput.pdf")


#####################################
# adjust input data and window size #
#####################################
window <- 7

R_consensus_adjInputWindow <- estimate_RKI_R(incid, method = "EpiEstim",
                                          window = window, delay = 0,
                                          gt_type = "gamma", gt_mean = 4, gt_sd = 4)

R_RKI_adjInputWindow <- estimate_RKI_R(RKI_incid, window = window)

R_SDSC_EpiEstim_adjInputWindow <- estimate_SDSC_R(SDSC_incid, window = window)

R_ETH_EpiEstim_adjInputWindow <- estimate_ETH_R(incid_for_ETH, window = window)

R_AGES_EpiEstim_adjInputWindow <- estimate_AGES_R(incid, window = window)

R_Ilmenau_adjInputWindow <- estimate_Ilmenau_R(incid, window = window)

R_epiforecasts_adjInputWindow <- R_epiforecasts_adjInput

R_globalrt_smoother_adjInputWindow <- R_globalrt_smoother_adjInput
#R_globalrt_filter_adjInputWindow <- R_globalrt_filter_adjInput

R_rtlive_adjInputWindow <- R_rtlive_adjInput

# merge estimates and plot for comparison
estimates_adjInputWindow <- R_consensus_adjInputWindow %>%
  full_join(R_RKI_adjInputWindow, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_adjInputWindow, by = "date") %>% 
  full_join(R_ETH_EpiEstim_adjInputWindow, by = "date") %>% 
  full_join(R_AGES_EpiEstim_adjInputWindow, by = "date") %>% 
  full_join(R_Ilmenau_adjInputWindow, by = "date") %>% 
  full_join(R_epiforecasts_adjInputWindow, by = "date") %>%
  full_join(R_globalrt_smoother_adjInputWindow, by = "date") %>%
  #full_join(R_globalrt_filter_adjInputWindow, by = "date") %>%
  full_join(R_rtlive_adjInputWindow, by = "date") %>%
  arrange(date)

comp_methods <- c("consensus", "RKI", "SDSC", "ETH",
                  "AGES", "Ilmenau", "epiforecasts",
                  "globalrt", #"globalrt filter",
                  "rtlive")
#plot_for_comparison(estimates_adjInputWindow, comp_methods, filenames = "_adjInputWindow.pdf")

estimates_adjInputWindow_CI <- as.data.frame(R_consensus_adjInputWindow) %>% 
  full_join(R_SDSC_EpiEstim_adjInputWindow, by = "date") %>%
  full_join(R_ETH_EpiEstim_adjInputWindow, by = "date") %>%
  full_join(R_AGES_EpiEstim_adjInputWindow, by = "date") %>%
  full_join(R_Ilmenau_adjInputWindow, by = "date") %>% 
  full_join(R_epiforecasts_adjInputWindow, by = "date") %>%
  full_join(R_globalrt_smoother_adjInputWindow, by = "date") %>%
  #full_join(R_globalrt_filter_adjInputWindow, by = "date") %>%
  full_join(R_rtlive_adjInputWindow, by = "date") %>%
  arrange(date)

comp_CI <- c("consensus", "SDSC", "ETH",
             "AGES", "Ilmenau", "epiforecasts",
             "globalrt", #"globalrt filter",
             "rtlive")
#plot_for_comparison(estimates_adjInputWindow_CI, comp_CI, include_CI = T, filenames = "_CI_adjInputWindow.pdf")


######################################################
# adjust input data, window size and generation time #
######################################################
window <- 7
gt_type <- "gamma"
gt_mean <- 4
gt_sd <- 4

R_consensus_adjInputWindowGTD <- estimate_RKI_R(incid, method = "EpiEstim",
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

path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
R_epiforecasts_adjInputWindowGTD <- qread(paste0(path, "R_calc_2021-07-10_final_adjInputWindowGTD.qs"))
R_epiforecasts_adjInputWindowGTD <- R_epiforecasts_adjInputWindowGTD[,c("date", "median", "lower_95", "upper_95")]
names(R_epiforecasts_adjInputWindowGTD) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
R_globalrt_smoother_adjInputWindowGTD <- read_csv(paste0(path, "bayesian_smoother_4.csv"))[,c("Date", "R", "lb_95", "ub_95")]
names(R_globalrt_smoother_adjInputWindowGTD) <- c("date", "R_calc", "lower", "upper")
#R_globalrt_filter_adjInputWindowGTD <- read_csv(paste0(path, "bayesian_filter_4.csv"))
#names(R_globalrt_filter_adjInputWindowGTD) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
R_rtlive_adjInputWindowGTD <- read.csv(paste0(path, "gamma_4_de_2021-07-10_all_summary.csv")) %>%
  dplyr::select(X, median, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "median", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))

# merge estimates and plot for comparison
estimates_adjInputWindowGTD <- R_consensus_adjInputWindowGTD %>%
  full_join(R_RKI_adjInputWindowGTD, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_adjInputWindowGTD, by = "date") %>% 
  full_join(R_ETH_EpiEstim_adjInputWindowGTD, by = "date") %>% 
  full_join(R_AGES_EpiEstim_adjInputWindowGTD, by = "date") %>% 
  full_join(R_Ilmenau_adjInputWindowGTD, by = "date") %>% 
  full_join(R_epiforecasts_adjInputWindowGTD, by = "date") %>%
  full_join(R_globalrt_smoother_adjInputWindowGTD, by = "date") %>%
  #full_join(R_globalrt_filter_adjInputWindowGTD, by = "date") %>%
  full_join(R_rtlive_adjInputWindowGTD, by = "date") %>%
  arrange(date)

comp_methods <- c("consensus", "RKI", "SDSC", "ETH",
                  "AGES", "Ilmenau", "epiforecasts",
                  "globalrt", #"globalrt filter",
                  "rtlive")
#plot_for_comparison(estimates_adjInputWindowGTD, comp_methods, filenames = "_adjInputWindowGTD.pdf")

estimates_adjInputWindowGTD_CI <- as.data.frame(R_consensus_adjInputWindowGTD) %>% 
  full_join(R_SDSC_EpiEstim_adjInputWindowGTD, by = "date") %>%
  full_join(R_ETH_EpiEstim_adjInputWindowGTD, by = "date") %>%
  full_join(R_AGES_EpiEstim_adjInputWindowGTD, by = "date") %>%
  full_join(R_Ilmenau_adjInputWindowGTD, by = "date") %>% 
  full_join(R_epiforecasts_adjInputWindowGTD, by = "date") %>%
  full_join(R_globalrt_smoother_adjInputWindowGTD, by = "date") %>%
  #full_join(R_globalrt_filter_adjInputWindowGTD, by = "date") %>%
  full_join(R_rtlive_adjInputWindowGTD, by = "date") %>%
  arrange(date)

comp_CI <- c("consensus", "SDSC", "ETH",
             "AGES", "Ilmenau", "epiforecasts",
             "globalrt", #"globalrt filter",
             "rtlive")
#plot_for_comparison(estimates_adjInputWindowGTD_CI, comp_CI, include_CI = T, filenames = "_CI_adjInputWindowGTD.pdf")


############################################
# adjust everything                        #
# input data, window size, gtd, mean delay #
############################################
window <- 7
gt_type <- "gamma"
gt_mean <- 4
gt_sd <- 4

R_consensus_adjAll <- estimate_RKI_R(incid, method = "EpiEstim",
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

path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- "R_calc_2021-07-10_final_AdjAll.qs"
R_epiforecasts_adjAll <- qread(paste0(path, file))
R_epiforecasts_adjAll <- R_epiforecasts_adjAll[,c("date", "median", "lower_95", "upper_95")]
names(R_epiforecasts_adjAll) <- c("date", "R_calc", "lower", "upper")

R_globalrt_smoother_adjAll <- R_globalrt_smoother_adjInputWindowGTD
#R_globalrt_filter_adjAll <- R_globalrt_filter_adjInputWindowGTD

R_rtlive_adjAll <- R_rtlive_adjInputWindowGTD
R_rtlive_adjAll$date <- R_rtlive_adjAll$date + params["rtlive", "delay"]

# merge estimates and plot for comparison
estimates_adjAll <- R_consensus_adjAll %>%
  full_join(R_RKI_adjAll, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_adjAll, by = "date") %>% 
  full_join(R_ETH_EpiEstim_adjAll, by = "date") %>% 
  full_join(R_AGES_EpiEstim_adjAll, by = "date") %>% 
  full_join(R_Ilmenau_adjAll, by = "date") %>% 
  full_join(R_epiforecasts_adjAll, by = "date") %>%
  full_join(R_globalrt_smoother_adjAll, by = "date") %>%
  #full_join(R_globalrt_filter_adjAll, by = "date") %>%
  full_join(R_rtlive_adjAll, by = "date") %>%
  arrange(date)

comp_methods <- c("consensus", "RKI", "SDSC", "ETH",
                  "AGES", "Ilmenau", "epiforecasts",
                  "globalrt", #"globalrt filter",
                  "rtlive")
#plot_for_comparison(estimates_adjAll, comp_methods, filenames = "_adjAll.pdf")

estimates_adjAll_CI <- as.data.frame(R_consensus_adjAll) %>% 
  full_join(R_SDSC_EpiEstim_adjAll, by = "date") %>%
  full_join(R_ETH_EpiEstim_adjAll, by = "date") %>%
  full_join(R_AGES_EpiEstim_adjAll, by = "date") %>%
  full_join(R_Ilmenau_adjAll, by = "date") %>% 
  full_join(R_epiforecasts_adjAll, by = "date") %>%
  full_join(R_globalrt_smoother_adjAll, by = "date") %>%
  #full_join(R_globalrt_filter_adjAll, by = "date") %>%
  full_join(R_rtlive_adjAll, by = "date") %>%
  arrange(date)

comp_CI <- c("consensus", "SDSC", "ETH",
             "AGES", "Ilmenau", "epiforecasts",
             "globalrt", #"globalrt filter",
             "rtlive")
#plot_for_comparison(estimates_adjAll_CI, comp_CI, include_CI = T, filenames = "_CI_adjAll.pdf")

