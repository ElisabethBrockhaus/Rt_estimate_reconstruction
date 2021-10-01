setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

method <- "ETH"

# parameter combinations used in papers
gt_dist <- c("gamma", "constant", "ad hoc", "gamma", "gamma", "gamma", "gamma", "lognorm", "exponential")
mean_gt <- c(4.8, 4, 5.6, 4.8, 5, 3.4, 3.6, 4.7, 7)
sd_gt <- c(2.3, 0, 4.2, 2.3, 4, 1.8, 3.1, 2.9, 7)
delay <- c(11, 1, 7, 10, 0, 0, 2, 12, 0)
source("Rt_estimate_reconstruction/ETH/delays_for_ETH_estimation.R")

params <- data.frame(gtd=gt_dist, gt_mean=mean_gt, gt_sd=sd_gt, delay=delay)
methods <- c("ETH", "RKI", "Ilmenau", "SDSC", "Zi", "AGES", "epiforecasts", "rtlive", "globalrt")
rownames(params) <- methods

# use RKI data for all methods (not line list!)
incid <- load_incidence_data(method = "RKI")

incid_for_ETH <- load_incidence_data(method = "ETHZ_sliding_window", source = "_simpleRKI",
                                     new_deconvolution = if (method == "ETH") FALSE else TRUE,
                                     delays = if (method == "ETH") list() else delays_ETH[[method]])

# save incidence data for epiforecast estimation
write_csv(incid, "Rt_estimate_reconstruction/incidence_data/RKI_incid.csv")

# for comparison of methods use default window size of each method
R_raw_EpiEstim <- estimate_RKI_R(incid, method = "EpiEstim",
                                 window = 7,
                                 gt_type = params[method, "gtd"],
                                 gt_mean = params[method, "gt_mean"],
                                 gt_sd = params[method, "gt_sd"],
                                 delay = params[method, "delay"])

R_ETH_EpiEstim <- estimate_ETH_R(incid_for_ETH,
                                 gt_type = params[method, "gtd"],
                                 gt_mean = params[method, "gt_mean"],
                                 gt_sd = params[method, "gt_sd"])

R_AGES_EpiEstim <- estimate_AGES_R(incid,
                                   gt_type = params[method, "gtd"],
                                   gt_mean=params[method, "gt_mean"],
                                   gt_sd=params[method, "gt_sd"],
                                   delay = params[method, "delay"])

R_Ilmenau <- estimate_Ilmenau_R(incid,
                                gt_type = params[method, "gtd"],
                                gt_mean=params[method, "gt_mean"],
                                gt_sd=params[method, "gt_sd"],
                                delay = params[method, "delay"])[,c("date", "0.5")]
names(R_Ilmenau)[2] <- "R_calc"

R_epiforecasts <- qread("Rt_estimate_reconstruction/epiforecasts/estimates/R_calc_2021-10-01_ETHParams.qs")
names(R_epiforecasts) <- c("date", "type", "R_calc")

R_globalrt <- read_csv(paste0("Rt_estimate_reconstruction/ArroyoMarioli/estimates/estimated_R_", method, ".csv"))
R_globalrt <- R_globalrt[R_globalrt$`Country/Region` == "Germany", c("Date", "R")]
names(R_globalrt) <- c("date", "R_calc")

# merge estimates and plot for comparison
estimates <- R_raw_EpiEstim[,c("date", "R_calc")] %>%
  full_join(R_ETH_EpiEstim[,c("date", "R_calc")], by = "date") %>% 
  full_join(R_AGES_EpiEstim[,c("date", "R_calc")], by = "date") %>% 
  full_join(R_Ilmenau[,c("date", "R_calc")], by = "date") %>% 
  full_join(R_epiforecasts[,c("date", "R_calc")], by = "date") %>%
  full_join(R_globalrt[,c("date", "R_calc")], by = "date")

plot_multiple_estimates(estimates[estimates$date > "2021-04-01",],
                        methods = c("raw EpiEstim", "ETH EpiEstim", "AGES EpiEstim",
                                    "Ilmenau", "epiforecasts", "globalrt"))


####################################
# in contrast: real-time estimates #
####################################
RKI_R_pub <- load_published_R_estimates(source = "RKI_7day")
ETH_R_pub <- load_published_R_estimates("ETHZ_sliding_window")
Ilmenau_R_pub <- load_published_R_estimates("ilmenau")
SDSC_R_pub <- load_published_R_estimates("sdsc")
globalrt_R_pub <- read_csv("https://raw.githubusercontent.com/crondonm/TrackingR/main/Estimates-Database/database.csv")
globalrt_R_pub <- globalrt_R_pub[globalrt_R_pub$`Country/Region` == "Germany" &
                                   globalrt_R_pub$days_infectious == 7, c("Date", "R")]
names(globalrt_R_pub)[1] <- "date"
epinow_R_pub <- load_published_R_estimates("epiforecasts")

# merge estimates and plot for comparison
estimates <- RKI_R_pub[,c("date", "R_pub")] %>%
  full_join(ETH_R_pub[,c("date", "R_pub")], by = "date") %>% 
  #full_join(Ilmenau_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(SDSC_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(globalrt_R_pub[,c("date", "R")], by = "date") %>%
  full_join(epinow_R_pub[,c("date", "R_pub")], by = "date")

plot_multiple_estimates(estimates[estimates$date > "2020-03-05",],
                        methods = c("RKI", "ETH", "SDSC", "globalrt", "epiforecasts"))


###########################################
# compare estimates of Cori-based methods #
###########################################

# estimations
RKI_R <- estimate_RKI_R(simple_incid,
                        window=params[method, "window"],
                        gt_type = params[method, "gtd"],
                        gt_mean=params[method, "gt_mean"],
                        gt_sd=params[method, "gt_sd"],
                        shift=1) #params["RKI", "shift"])

ETH_R <- estimate_ETH_R(incid,
                        window=params[method, "window"],
                        gt_type = params[method, "gtd"],
                        gt_mean=params[method, "gt_mean"],
                        gt_sd=params[method, "gt_sd"],
                        shift=params["ETH", "shift"])

Ilmenau_R <- estimate_Ilmenau_R(simple_incid,
                                window=params[method, "window"],
                                gt_type = params[method, "gtd"],
                                gt_mean=params[method, "gt_mean"],
                                gt_sd=params[method, "gt_sd"],
                                shift=params["Ilmenau", "shift"])

AGES_R <- estimate_AGES_R(simple_incid,
                          window=params[method, "window"],
                          gt_type = params[method, "gtd"],
                          gt_mean=params[method, "gt_mean"],
                          gt_sd=params[method, "gt_sd"],
                          shift=params["AGES", "shift"])

SDSC_R <- estimate_SDSC_R(simple_incid,
                          window=params[method, "window"],
                          gt_type = params[method, "gtd"],
                          gt_mean=params[method, "gt_mean"],
                          gt_sd=params[method, "gt_sd"],
                          shift=params["SDSC", "shift"])

# merge estimates and plot for comparison
estimates <- RKI_R %>%
  full_join(ETH_R, by = "date") %>% 
  full_join(Ilmenau_R, by = "date") %>%
  full_join(AGES_R, by = "date") %>%
  full_join(SDSC_R, by = "date")
names(estimates) <- c("date", "RKI", "ETH", "Ilmenau", "AGES", "SDSC")

plot_multiple_estimates(estimates)






