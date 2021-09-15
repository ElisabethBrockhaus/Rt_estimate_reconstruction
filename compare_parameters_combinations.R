setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

# parameter combinations used in papers
window_size <- c(3, 7, 1, 4, 13, 1, 1, 1)
gt_dist <- c("gamma", "constant", "ad hoc", "gamma", "gamma", "gamma", "log-normal", "exponential")
mean_gt <- c(4.8, 4.0, 5.6, 4.8, 3.4, 3.6, 4.7, 7)
sd_gt <- c(2.3, 0.0, 4.2, 2.3, 1.8, 3.1, 2.9, 7)
shift <- c(0, 4, 7, 7, 0, 2, 5, 0)

params <- data.frame(window=window_size, gtd=gt_dist, gt_mean=mean_gt, gt_sd=sd_gt, shift=shift)
rownames(params) <- c("ETH", "RKI", "Ilmenau", "SDSC", "AGES", "epiforecasts", "rtlive", "globalrt")

# use same data for all methods
incid <- load_incidence_data(method = "ETHZ_sliding_window")
simple_incid <- incid[, c("date", "value")]
simple_incid <- aggregate(simple_incid$value, by=list(simple_incid$date), FUN=median)
names(simple_incid) <- c("date", "I")

method <- "globalrt"

# for comparison of methods use window = 7
R_raw_EpiEstim <- estimate_RKI_R(simple_incid, method = "EpiEstim",
                                 window = 7,
                                 gt_type = params[method, "gtd"],
                                 gt_mean=params[method, "gt_mean"],
                                 gt_sd=params[method, "gt_sd"])

R_ETH_EpiEstim <- estimate_ETH_R(incid,
                                 window = 7,
                                 gt_type = params[method, "gtd"],
                                 gt_mean=params[method, "gt_mean"],
                                 gt_sd=params[method, "gt_sd"])

R_epiforecasts <- qread("Rt_estimate_reconstruction/epiforecasts/estimates/R_calc_2021-09-15_globalrtParams.qs")
names(R_epiforecasts) <- c("date", "type", "R_calc")

R_globalrt <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/estimates/estimated_R_STAN.csv")
R_globalrt <- R_globalrt[R_globalrt$`Country/Region` == "Germany", c("Date", "R")]
names(R_globalrt) <- c("date", "R_calc")

# merge estimates and plot for comparison
estimates <- R_raw_EpiEstim[,c("date", "R_calc")] %>%
  full_join(R_ETH_EpiEstim[,c("date", "R_calc")], by = "date") %>% 
  full_join(R_epiforecasts[,c("date", "R_calc")], by = "date") %>%
  full_join(R_globalrt[,c("date", "R_calc")], by = "date")

plot_multiple_estimates(estimates, methods = c("raw EpiEstim", "ETH EpiEstim", "epiforecasts", "globalrt"))




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






