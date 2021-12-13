setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()

# path for saving the estimates
path_estimates <- "Rt_estimate_reconstruction/estimates/"


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
RKI_R_adjInput <- read_csv("https://raw.githubusercontent.com/robert-koch-institut/SARS-CoV-2-Nowcasting_und_-R-Schaetzung/main/Archiv/Nowcast_R_2021-07-10.csv",
                      col_types = list(Datum = col_date())) %>%
  dplyr::select("Datum", "PS_7_Tage_R_Wert", "UG_PI_7_Tage_R_Wert", "OG_PI_7_Tage_R_Wert") %>%
  rename(c("date" = "Datum", "R_adjInput" = "PS_7_Tage_R_Wert", "lower" = "UG_PI_7_Tage_R_Wert", "upper" = "OG_PI_7_Tage_R_Wert"))
ETH_R_adjInput <- load_adjInputlished_R_estimates("ETHZ_sliding_window")
Ilmenau_R_adjInput <- load_adjInputlished_R_estimates("ilmenau")
SDSC_R_adjInput <- load_adjInputlished_R_estimates("sdsc")
globalrt_R_adjInput <- load_adjInputlished_R_estimates("globalrt_7d")
epiforecasts_R_adjInput <- load_adjInputlished_R_estimates("epiforecasts")
rtlive_R_adjInput <- load_adjInputlished_R_estimates("rtlive")

# merge estimates and plot for comparison
estimates_adjInput <- RKI_R_adjInput[,c("date", "R_adjInput")] %>% rename(RKI = R_adjInput) %>%
  full_join(ETH_R_adjInput[,c("date", "R_adjInput")] %>% rename(ETH = R_adjInput), by = "date") %>% 
  full_join(Ilmenau_R_adjInput[,c("date", "R_adjInput")] %>% rename(Ilmenau = R_adjInput), by = "date") %>% 
  full_join(SDSC_R_adjInput[,c("date", "R_adjInput")] %>% rename(SDSC = R_adjInput), by = "date") %>% 
  full_join(globalrt_R_adjInput[,c("date", "R_adjInput")] %>% rename(globalrt = R_adjInput), by = "date") %>%
  full_join(epiforecasts_R_adjInput[,c("date", "R_adjInput")] %>% rename(epiforecasts = R_adjInput), by = "date") %>%
  full_join(rtlive_R_adjInput[,c("date", "R_adjInput")] %>% rename(rtlive = R_adjInput), by = "date")

write_csv(estimates_adjInput, paste0(path_estimates, "R_adjInput_2021-07-10.csv"))

estimates_adjInput_ci <- RKI_R_adjInput %>% rename(R.RKI = R_adjInput, lower.RKI = lower, upper.RKI = upper) %>%
  full_join(ETH_R_adjInput %>% rename(R.ETH = R_adjInput, lower.ETH = lower, upper.ETH = upper), by = "date") %>%
  full_join(Ilmenau_R_adjInput %>% rename(R.Ilmenau = R_adjInput, lower.Ilmenau = lower, upper.Ilmenau = upper), by = "date") %>%
  full_join(SDSC_R_adjInput %>% rename(R.SDSC = R_adjInput, lower.SDSC = lower, upper.SDSC = upper), by = "date") %>%
  full_join(globalrt_R_adjInput %>% rename(R.globalrt = R_adjInput, lower.globalrt = lower, upper.globalrt = upper), by = "date") %>%
  full_join(epiforecasts_R_adjInput %>% rename(R.epiforecasts = R_adjInput, lower.epiforecasts = lower, upper.epiforecasts = upper), by = "date") %>%
  full_join(rtlive_R_adjInput %>% rename(R.rtlive = R_adjInput, lower.rtlive = lower, upper.rtlive = upper), by = "date")

write_csv(estimates_adjInput_ci, paste0(path_estimates, "R_adjInput_CI_2021-07-10.csv"))


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

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
R_rtlive_adjInput <- read.csv(paste0(path, "DE_2021-07-10_all_trace_summary.csv")) %>%
  dplyr::select(X, median, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "median", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))

estimates_adjInput <- R_consensus_adjInput[,c("date", "R_calc")] %>% rename(consensus = R_calc) %>%
  full_join(R_RKI_adjInput[,c("date", "R_calc")] %>% rename(RKI = R_calc), by = "date") %>% 
  full_join(R_SDSC_EpiEstim_adjInput[,c("date", "R_calc")] %>% rename(SDSC = R_calc), by = "date") %>% 
  full_join(R_ETH_EpiEstim_adjInput[,c("date", "R_calc")] %>% rename(ETH = R_calc), by = "date") %>% 
  full_join(R_AGES_EpiEstim_adjInput[,c("date", "R_calc")] %>% rename(AGES = R_calc), by = "date") %>% 
  full_join(R_Ilmenau_adjInput[,c("date", "R_calc")] %>% rename(Ilmenau = R_calc), by = "date") %>% 
  full_join(R_epiforecasts_adjInput[,c("date", "R_calc")] %>% rename(epiforecasts = R_calc), by = "date") %>%
  full_join(R_globalrt_smoother_adjInput[,c("date", "R_calc")] %>% rename(globalrt = R_calc), by = "date") %>%
  full_join(R_rtlive_adjInput[,c("date", "R_calc")] %>% rename(rtlive = R_calc), by = "date") %>%
  arrange(date)

write_csv(estimates_adjInput, paste0(path_estimates, "R_adjInput_2021-07-10.csv"))

estimates_adjInput_ci <- R_consensus_adjInput %>% rename(R.consensus = R_calc, lower.consensus = lower, upper.consensus = upper) %>%
  full_join(R_SDSC_EpiEstim_adjInput %>% rename(R.SDSC = R_calc, lower.SDSC = lower, upper.SDSC = upper), by = "date") %>%
  full_join(R_ETH_EpiEstim_adjInput %>% rename(R.ETH = R_calc, lower.ETH = lower, upper.ETH = upper), by = "date") %>%
  full_join(R_AGES_EpiEstim_adjInput %>% rename(R.AGES = R_calc, lower.AGES = lower, upper.AGES = upper), by = "date") %>%
  full_join(R_Ilmenau_adjInput %>% rename(R.Ilmenau = R_calc, lower.Ilmenau = lower, upper.Ilmenau = upper), by = "date") %>%
  full_join(R_epiforecasts_adjInput %>% rename(R.epiforecasts = R_calc, lower.epiforecasts = lower, upper.epiforecasts = upper), by = "date") %>%
  full_join(R_globalrt_smoother_adjInput %>% rename(R.globalrt = R_calc, lower.globalrt = lower, upper.globalrt = upper), by = "date") %>%
  full_join(R_rtlive_adjInput %>% rename(R.rtlive = R_calc, lower.rtlive = lower, upper.rtlive = upper), by = "date")

write_csv(estimates_adjInput_ci, paste0(path_estimates, "R_adjInput_CI_2021-07-10.csv"))


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

R_rtlive_adjInputWindow <- R_rtlive_adjInput

estimates_adjInputWindow <- R_consensus_adjInputWindow[,c("date", "R_calc")] %>% rename(consensus = R_calc) %>%
  full_join(R_RKI_adjInputWindow[,c("date", "R_calc")] %>% rename(RKI = R_calc), by = "date") %>% 
  full_join(R_SDSC_EpiEstim_adjInputWindow[,c("date", "R_calc")] %>% rename(SDSC = R_calc), by = "date") %>% 
  full_join(R_ETH_EpiEstim_adjInputWindow[,c("date", "R_calc")] %>% rename(ETH = R_calc), by = "date") %>% 
  full_join(R_AGES_EpiEstim_adjInputWindow[,c("date", "R_calc")] %>% rename(AGES = R_calc), by = "date") %>% 
  full_join(R_Ilmenau_adjInputWindow[,c("date", "R_calc")] %>% rename(Ilmenau = R_calc), by = "date") %>% 
  full_join(R_epiforecasts_adjInputWindow[,c("date", "R_calc")] %>% rename(epiforecasts = R_calc), by = "date") %>%
  full_join(R_globalrt_smoother_adjInputWindow[,c("date", "R_calc")] %>% rename(globalrt = R_calc), by = "date") %>%
  full_join(R_rtlive_adjInputWindow[,c("date", "R_calc")] %>% rename(rtlive = R_calc), by = "date") %>%
  arrange(date)

write_csv(estimates_adjInputWindow, paste0(path_estimates, "R_adjInputWindow_2021-07-10.csv"))

estimates_adjInputWindow_ci <- R_consensus_adjInputWindow %>% rename(R.consensus = R_calc, lower.consensus = lower, upper.consensus = upper) %>%
  full_join(R_SDSC_EpiEstim_adjInputWindow %>% rename(R.SDSC = R_calc, lower.SDSC = lower, upper.SDSC = upper), by = "date") %>%
  full_join(R_ETH_EpiEstim_adjInputWindow %>% rename(R.ETH = R_calc, lower.ETH = lower, upper.ETH = upper), by = "date") %>%
  full_join(R_AGES_EpiEstim_adjInputWindow %>% rename(R.AGES = R_calc, lower.AGES = lower, upper.AGES = upper), by = "date") %>%
  full_join(R_Ilmenau_adjInputWindow %>% rename(R.Ilmenau = R_calc, lower.Ilmenau = lower, upper.Ilmenau = upper), by = "date") %>%
  full_join(R_epiforecasts_adjInputWindow %>% rename(R.epiforecasts = R_calc, lower.epiforecasts = lower, upper.epiforecasts = upper), by = "date") %>%
  full_join(R_globalrt_smoother_adjInputWindow %>% rename(R.globalrt = R_calc, lower.globalrt = lower, upper.globalrt = upper), by = "date") %>%
  full_join(R_rtlive_adjInputWindow %>% rename(R.rtlive = R_calc, lower.rtlive = lower, upper.rtlive = upper), by = "date")

write_csv(estimates_adjInputWindow_ci, paste0(path_estimates, "R_adjInputWindow_CI_2021-07-10.csv"))


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

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
R_rtlive_adjInputWindowGTD <- read.csv(paste0(path, "gamma_4_de_2021-07-10_all_summary.csv")) %>%
  dplyr::select(X, median, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "median", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))

estimates_adjInputWindowGTD <- R_consensus_adjInputWindowGTD[,c("date", "R_calc")] %>% rename(consensus = R_calc) %>%
  full_join(R_RKI_adjInputWindowGTD[,c("date", "R_calc")] %>% rename(RKI = R_calc), by = "date") %>% 
  full_join(R_SDSC_EpiEstim_adjInputWindowGTD[,c("date", "R_calc")] %>% rename(SDSC = R_calc), by = "date") %>% 
  full_join(R_ETH_EpiEstim_adjInputWindowGTD[,c("date", "R_calc")] %>% rename(ETH = R_calc), by = "date") %>% 
  full_join(R_AGES_EpiEstim_adjInputWindowGTD[,c("date", "R_calc")] %>% rename(AGES = R_calc), by = "date") %>% 
  full_join(R_Ilmenau_adjInputWindowGTD[,c("date", "R_calc")] %>% rename(Ilmenau = R_calc), by = "date") %>% 
  full_join(R_epiforecasts_adjInputWindowGTD[,c("date", "R_calc")] %>% rename(epiforecasts = R_calc), by = "date") %>%
  full_join(R_globalrt_smoother_adjInputWindowGTD[,c("date", "R_calc")] %>% rename(globalrt = R_calc), by = "date") %>%
  full_join(R_rtlive_adjInputWindowGTD[,c("date", "R_calc")] %>% rename(rtlive = R_calc), by = "date") %>%
  arrange(date)

write_csv(estimates_adjInputWindowGTD, paste0(path_estimates, "R_adjInputWindowGTD_2021-07-10.csv"))

estimates_adjInputWindowGTD_ci <- R_consensus_adjInputWindowGTD %>% rename(R.consensus = R_calc, lower.consensus = lower, upper.consensus = upper) %>%
  full_join(R_SDSC_EpiEstim_adjInputWindowGTD %>% rename(R.SDSC = R_calc, lower.SDSC = lower, upper.SDSC = upper), by = "date") %>%
  full_join(R_ETH_EpiEstim_adjInputWindowGTD %>% rename(R.ETH = R_calc, lower.ETH = lower, upper.ETH = upper), by = "date") %>%
  full_join(R_AGES_EpiEstim_adjInputWindowGTD %>% rename(R.AGES = R_calc, lower.AGES = lower, upper.AGES = upper), by = "date") %>%
  full_join(R_Ilmenau_adjInputWindowGTD %>% rename(R.Ilmenau = R_calc, lower.Ilmenau = lower, upper.Ilmenau = upper), by = "date") %>%
  full_join(R_epiforecasts_adjInputWindowGTD %>% rename(R.epiforecasts = R_calc, lower.epiforecasts = lower, upper.epiforecasts = upper), by = "date") %>%
  full_join(R_globalrt_smoother_adjInputWindowGTD %>% rename(R.globalrt = R_calc, lower.globalrt = lower, upper.globalrt = upper), by = "date") %>%
  full_join(R_rtlive_adjInputWindowGTD %>% rename(R.rtlive = R_calc, lower.rtlive = lower, upper.rtlive = upper), by = "date")

write_csv(estimates_adjInputWindowGTD_ci, paste0(path_estimates, "R_adjInputWindowGTD_CI_2021-07-10.csv"))


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
                                          estimateOffsetting = 0)

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

R_rtlive_adjAll <- R_rtlive_adjInputWindowGTD
R_rtlive_adjAll$date <- R_rtlive_adjAll$date + params["rtlive", "delay"]

estimates_adjAll <- R_consensus_adjAll[,c("date", "R_calc")] %>% rename(consensus = R_calc) %>%
  full_join(R_RKI_adjAll[,c("date", "R_calc")] %>% rename(RKI = R_calc), by = "date") %>% 
  full_join(R_SDSC_EpiEstim_adjAll[,c("date", "R_calc")] %>% rename(SDSC = R_calc), by = "date") %>% 
  full_join(R_ETH_EpiEstim_adjAll[,c("date", "R_calc")] %>% rename(ETH = R_calc), by = "date") %>% 
  full_join(R_AGES_EpiEstim_adjAll[,c("date", "R_calc")] %>% rename(AGES = R_calc), by = "date") %>% 
  full_join(R_Ilmenau_adjAll[,c("date", "R_calc")] %>% rename(Ilmenau = R_calc), by = "date") %>% 
  full_join(R_epiforecasts_adjAll[,c("date", "R_calc")] %>% rename(epiforecasts = R_calc), by = "date") %>%
  full_join(R_globalrt_smoother_adjAll[,c("date", "R_calc")] %>% rename(globalrt = R_calc), by = "date") %>%
  full_join(R_rtlive_adjAll[,c("date", "R_calc")] %>% rename(rtlive = R_calc), by = "date") %>%
  arrange(date)

write_csv(estimates_adjAll, paste0(path_estimates, "R_adjInputWindowGTDDelays_2021-07-10.csv"))

estimates_adjAll_ci <- R_consensus_adjAll %>% rename(R.consensus = R_calc, lower.consensus = lower, upper.consensus = upper) %>%
  full_join(R_SDSC_EpiEstim_adjAll %>% rename(R.SDSC = R_calc, lower.SDSC = lower, upper.SDSC = upper), by = "date") %>%
  full_join(R_ETH_EpiEstim_adjAll %>% rename(R.ETH = R_calc, lower.ETH = lower, upper.ETH = upper), by = "date") %>%
  full_join(R_AGES_EpiEstim_adjAll %>% rename(R.AGES = R_calc, lower.AGES = lower, upper.AGES = upper), by = "date") %>%
  full_join(R_Ilmenau_adjAll %>% rename(R.Ilmenau = R_calc, lower.Ilmenau = lower, upper.Ilmenau = upper), by = "date") %>%
  full_join(R_epiforecasts_adjAll %>% rename(R.epiforecasts = R_calc, lower.epiforecasts = lower, upper.epiforecasts = upper), by = "date") %>%
  full_join(R_globalrt_smoother_adjAll %>% rename(R.globalrt = R_calc, lower.globalrt = lower, upper.globalrt = upper), by = "date") %>%
  full_join(R_rtlive_adjAll %>% rename(R.rtlive = R_calc, lower.rtlive = lower, upper.rtlive = upper), by = "date")

write_csv(estimates_adjAll_ci, paste0(path_estimates, "R_adjInputWindowGTDDelays_CI_2021-07-10.csv"))


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

