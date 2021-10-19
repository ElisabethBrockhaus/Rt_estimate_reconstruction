library(readr)

setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()


###########################
# run globalrt estimation #
###########################

library(reticulate)

# refresh data for globalrt
source("Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/format_data_from_various_sources.R")
envVars <- ls(name = 1)
keepObjects <- c()
rm(list = setdiff(envVars, keepObjects), pos = 1)

py_run_file("Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/construct_dataset_adj.py")
py_run_file("Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/estimate_R_KF_adj.py")



################################
# get functions and parameters #
################################

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

# decide what evaluation steps should be done
plot_all_differences <- FALSE

# parameter combinations used in papers
gt_dist <- c("gamma", "constant", "ad hoc", "gamma", "gamma", "gamma", "gamma", "lognorm", "exponential")
mean_gt <- c(4.8,      4,          5.6,      4.8,     5,       3.4,     3.6,     4.7,       7)
sd_gt <-   c(2.3,      0,          4.2,      2.3,     4,       1.8,     3.1,     2.9,       7)
delay <-   c(11,       1,          7,        10,      0,       0,       12,      12,        0)

params <- data.frame(gtd=gt_dist, gt_mean=mean_gt, gt_sd=sd_gt, delay=delay)
methods <- c("ETH", "RKI", "Ilmenau", "SDSC", "Zi", "AGES", "epiforecasts", "rtlive", "globalrt")
rownames(params) <- methods

source("Rt_estimate_reconstruction/ETH/delays_for_ETH_estimation.R")

# use RKI Nowcast data for all methods
RKI_incid <- load_incidence_data(method = "RKI")
RKI_incid <- RKI_incid[RKI_incid$date < "2021-10-01",]

# use incidence data used by rtlive (sourced from RKI line-list data)
rtlive_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_data_21_10_19.csv")
rtlive_incid <- rtlive_incid[rtlive_incid$region=="all", c("date", "new_cases")]
names(rtlive_incid) <- c("date", "I")

# compare incidence time series used by RKI vs. rtlive
plot(RKI_incid, type="l")
lines(rtlive_incid, col="blue")

# save incidence data for epiforecast estimation
write_csv(RKI_incid, "Rt_estimate_reconstruction/incidence_data/RKI_incid.csv")
write_csv(rtlive_incid, "Rt_estimate_reconstruction/incidence_data/rtlive_incid.csv")

# choose data for comparison
incid <- rtlive_incid

# choose method for comparison
method <- "globalrt"

incid_for_ETH <- load_incidence_data(method = "ETHZ_sliding_window", source = "_simpleRKI",
                                     new_deconvolution = if (method == "ETH") FALSE else TRUE,
                                     delays = delays_ETH[[method]],
                                     delay_source = paste0("_", method))

# for comparison of methods use default window size of each method
R_raw_EpiEstim <- estimate_RKI_R(incid, method = "EpiEstim",
                                 window = 7,
                                 gt_type = params[method, "gtd"],
                                 gt_mean = params[method, "gt_mean"],
                                 gt_sd = params[method, "gt_sd"],
                                 delay = params[method, "delay"])

incid_ma7 <- incid
incid_ma7$I <- as.numeric(stats::filter(incid, rep(1/7, 7), method = "convolution", sides = 1)[,2])

R_raw_EpiEstim_ma7 <- estimate_RKI_R(incid_ma7[7:583,], method = "EpiEstim",
                                     window = 1,
                                     gt_type = params[method, "gtd"],
                                     gt_mean = params[method, "gt_mean"],
                                     gt_sd = params[method, "gt_sd"],
                                     delay = params[method, "delay"])

plot(R_raw_EpiEstim$date, R_raw_EpiEstim$R_calc, type="l")
lines(R_raw_EpiEstim_ma7$date, R_raw_EpiEstim_ma7$R_calc, col="red")

R_ma7 <- estimate_RKI_R(incid, method = "EpiEstim",
                        window = 1,
                        gt_type = params[method, "gtd"],
                        gt_mean = params[method, "gt_mean"],
                        gt_sd = params[method, "gt_sd"],
                        delay = params[method, "delay"])[,c("date", "R_calc")]
R_ma7$R_calc <- as.numeric(stats::filter(R_ma7, rep(1/7, 7), method = "convolution", sides = 1)[,2])
lines(R_ma7$date, R_ma7$R_calc, col="blue")


lines(R_raw_EpiEstim$date, R_raw_EpiEstim$lower, alpha=.5)
lines(R_raw_EpiEstim_ma7$date, R_raw_EpiEstim_ma7$lower, col="red")

x <- rpois(100, 10)
plot(x)
# apply a kernel smoothing:
smooth_kernel <- function(x, kernel){
  y <- rep(NA, length(x))
  for(i in length(kernel):length(x)){
    y[i] <- sum(kernel*x[i - 0:(length(kernel) - 1)])
  }
  y
}
# kernel smoothing is commutative
y1 <- smooth_kernel(x, kernel = c(0.4, 0.3, 0.2, 0.1))
y2 <- smooth_kernel(x, kernel = rep(1/7, 7))
z1 <- smooth_kernel(y2, kernel = c(0.4, 0.3, 0.2, 0.1))
z2 <- smooth_kernel(y1, kernel = rep(1/7, 7))
plot(z1, type = "l")
lines(z2, type = "l", col = "red", lty = 2)
plot(z1 - z2)

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
                                delay = params[method, "delay"])
names(R_Ilmenau) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- max(list.files(path, pattern = paste0("R_calc_\\d{4}-\\d{2}-\\d{2}_", method, "Params.qs")))
print(paste("Most recent epiforecasts estimates with matching parameters:", file))
R_epiforecasts <- qread(paste0(path, file))
R_epiforecasts <- R_epiforecasts[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
file <- paste0("estimated_R_", method, ".csv")
R_globalrt <- read_csv(paste0(path, file))
R_globalrt <- R_globalrt[R_globalrt$`Country/Region` == "Germany", c("Date", "R", "ci_95_l", "ci_95_u")]
names(R_globalrt) <- c("date", "R_calc", "lower", "upper")

#R_globalrt_STAN <- read_csv(paste0(path, "estimated_R_STAN.csv"))
#R_globalrt_KF <- read_csv(paste0(path, "estimated_R_KF.csv"))

# merge estimates and plot for comparison
estimates_allPars <- R_raw_EpiEstim %>%
  full_join(R_ETH_EpiEstim, by = "date") %>% 
  full_join(R_AGES_EpiEstim, by = "date") %>% 
  full_join(R_Ilmenau, by = "date") %>% 
  #full_join(R_epiforecasts, by = "date") %>%
  full_join(R_globalrt, by = "date")

comp_methods <- c("raw EpiEstim", "ETH EpiEstim", "AGES EpiEstim", "Ilmenau", #"epiforecasts",
                  "globalrt")
plot_for_comparison(estimates_allPars, comp_methods, method = method, variation = "parameters")

estimates_allPars_CI <- as.data.frame(R_ETH_EpiEstim) %>% 
  full_join(R_Ilmenau, by = "date") %>% 
  #full_join(R_epiforecasts, by = "date") %>%
  full_join(R_globalrt, by = "date")

comp_CI <- c("ETH EpiEstim", "Ilmenau", #"epiforecasts",
             "globalrt")
plot_for_comparison(estimates_allPars_CI, comp_CI, include_CI = T,
                    method = method, variation = "parameters")



####################################
# in contrast: real-time estimates #
####################################
RKI_R_pub <- load_published_R_estimates(source = "RKI_7day")
ETH_R_pub <- load_published_R_estimates("ETHZ_sliding_window")
Ilmenau_R_pub <- load_published_R_estimates("ilmenau")
SDSC_R_pub <- load_published_R_estimates("sdsc")
globalrt_R_pub <- load_published_R_estimates("globalrt_7d")
epinow_R_pub <- load_published_R_estimates("epiforecasts")
rtlive_R_pub <- load_published_R_estimates("rtlive")

# merge estimates and plot for comparison
estimates_pub <- RKI_R_pub[,c("date", "R_pub")] %>%
  full_join(ETH_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(Ilmenau_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(SDSC_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(globalrt_R_pub[,c("date", "R_pub")], by = "date") %>%
  #full_join(epinow_R_pub[,c("date", "R_pub")], by = "date") %>%
  full_join(rtlive_R_pub[,c("date", "R_pub")], by = "date")
org_methods <- c("RKI", "ETH", "Ilmenau", "SDSC", "globalrt",
                 #"epiforecasts",
                 "rtlive")

plot_for_comparison(estimates_pub, org_methods, method = "original", variation = "parameters")

estimates_pub_ci <- ETH_R_pub %>% 
  full_join(Ilmenau_R_pub, by = "date") %>% 
  full_join(globalrt_R_pub, by = "date") %>%
  full_join(epinow_R_pub, by = "date") #%>%
#full_join(rtlive_R_pub, by = "date")
org_methods <- c("ETH", "Ilmenau", "globalrt", "epiforecasts") #, "rtlive")

source("Rt_estimate_reconstruction/prepared_plots.R")
plot_for_comparison(estimates_pub_ci, org_methods, include_CI=T,
                    method = "original", variation = "parameters")

######################
# adjust delays only #
######################

# for comparison of methods use default window size of each method
R_raw_EpiEstim_d <- estimate_RKI_R(incid, method = "EpiEstim",
                                 window = 7,
                                 delay = params[method, "delay"])

R_ETH_EpiEstim_d <- estimate_ETH_R(incid_for_ETH)

R_AGES_EpiEstim_d <- estimate_AGES_R(incid,
                                     delay = params[method, "delay"])

R_Ilmenau_d <- estimate_Ilmenau_R(incid,
                                delay = params[method, "delay"])
names(R_Ilmenau_d) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- max(list.files(path, pattern = paste0("R_calc_\\d{4}-\\d{2}-\\d{2}_", method, "_delays.qs")))
print(paste("Most recent epiforecasts estimates with matching parameters:", file))
R_epiforecasts_d <- qread(paste0(path, file))
R_epiforecasts_d <- R_epiforecasts_d[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_d) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
file <- paste0("estimated_R_", method, "_delays.csv")
R_globalrt_d <- read_csv(paste0(path, file))
R_globalrt_d <- R_globalrt_d[R_globalrt_d$`Country/Region` == "Germany", c("Date", "R", "ci_95_l", "ci_95_u")]
names(R_globalrt_d) <- c("date", "R_calc", "lower", "upper")

# merge estimates and plot for comparison
estimates_delays <- R_raw_EpiEstim_d %>%
  full_join(R_ETH_EpiEstim_d, by = "date") %>% 
  full_join(R_AGES_EpiEstim_d, by = "date") %>% 
  full_join(R_Ilmenau_d, by = "date") %>% 
  full_join(R_epiforecasts_d, by = "date") %>%
  full_join(R_globalrt_d, by = "date")

comp_methods <- c("raw EpiEstim", "ETH EpiEstim", "AGES EpiEstim", "Ilmenau", "epiforecasts", "globalrt")
plot_for_comparison(estimates_delays, comp_methods, method = method, variation = "delays")

estimates_delays_CI <- as.data.frame(R_ETH_EpiEstim_d) %>% 
  full_join(R_Ilmenau_d, by = "date") %>% 
  full_join(R_epiforecasts_d, by = "date") %>%
  full_join(R_globalrt_d, by = "date")

comp_CI <- c("ETH EpiEstim", "Ilmenau", "epiforecasts", "globalrt")
plot_for_comparison(estimates_delays_CI, comp_CI, include_CI = T,
                    method = method, variation = "delays")



###############################
# adjust generation time only #
###############################

# for comparison of methods use default window size of each method
R_raw_EpiEstim_gt <- estimate_RKI_R(incid, method = "EpiEstim",
                                    window = 7,
                                    gt_type = params[method, "gtd"],
                                    gt_mean = params[method, "gt_mean"],
                                    gt_sd = params[method, "gt_sd"],
                                    delay = params["RKI", "delay"])

incid_ETH <- load_incidence_data(method = "ETHZ_sliding_window", source = "_simpleRKI",
                                 #new_deconvolution = TRUE,
                                 delays = delays_ETH[["ETH"]])

R_ETH_EpiEstim_gt <- estimate_ETH_R(incid_ETH,
                                    gt_type = params[method, "gtd"],
                                    gt_mean = params[method, "gt_mean"],
                                    gt_sd = params[method, "gt_sd"])

R_AGES_EpiEstim_gt <- estimate_AGES_R(incid,
                                      gt_type = params[method, "gtd"],
                                      gt_mean=params[method, "gt_mean"],
                                      gt_sd=params[method, "gt_sd"],
                                      delay = params["AGES", "delay"])

R_Ilmenau_gt <- estimate_Ilmenau_R(incid,
                                   gt_type = params[method, "gtd"],
                                   gt_mean=params[method, "gt_mean"],
                                   gt_sd=params[method, "gt_sd"],
                                   delay = params["Ilmenau", "delay"])
names(R_Ilmenau_gt) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- max(list.files(path, pattern = paste0("R_calc_\\d{4}-\\d{2}-\\d{2}_", method, "_GTD.qs")))
print(paste("Most recent epiforecasts estimates with matching parameters:", file))
R_epiforecasts_gt <- qread(paste0(path, file))
R_epiforecasts_gt <- R_epiforecasts_gt[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_gt) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
file <- paste0("estimated_R_", method, "_GTD.csv")
R_globalrt_gt <- read_csv(paste0(path, file))
R_globalrt_gt <- R_globalrt_gt[R_globalrt_gt$`Country/Region` == "Germany", c("Date", "R", "ci_95_l", "ci_95_u")]
names(R_globalrt_gt) <- c("date", "R_calc", "lower", "upper")

# merge estimates and plot for comparison
estimates_GTD <- R_raw_EpiEstim_gt %>%
  full_join(R_ETH_EpiEstim_gt, by = "date") %>% 
  full_join(R_AGES_EpiEstim_gt, by = "date") %>% 
  full_join(R_Ilmenau_gt, by = "date") %>% 
  full_join(R_epiforecasts_gt, by = "date") %>%
  full_join(R_globalrt_gt, by = "date")

comp_methods <- c("raw EpiEstim", "ETH EpiEstim", "AGES EpiEstim", "Ilmenau", "epiforecasts", "globalrt")
plot_for_comparison(estimates_GTD, comp_methods, method = method, variation = "generation time")

estimates_GTD_CI <- as.data.frame(R_ETH_EpiEstim_gt) %>% 
  full_join(R_Ilmenau_gt, by = "date") %>% 
  full_join(R_epiforecasts_gt, by = "date") %>%
  full_join(R_globalrt_gt, by = "date")

comp_CI <- c("ETH EpiEstim", "Ilmenau", "epiforecasts", "globalrt")
plot_for_comparison(estimates_GTD_CI, comp_CI, include_CI = T,
                    method = method, variation = "generation time")








