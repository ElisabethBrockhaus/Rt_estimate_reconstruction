library(readr)

setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()


source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")


####################################
# compare real-time estimates #
####################################
RKI_R_pub <- load_published_R_estimates(source = "RKI_7day")
ETH_R_pub <- load_published_R_estimates("ETHZ_sliding_window")
Ilmenau_R_pub <- load_published_R_estimates("ilmenau")
SDSC_R_pub <- load_published_R_estimates("sdsc")
Zi_R_pub <- load_published_R_estimates("zidatalab")
globalrt_R_pub <- load_published_R_estimates("globalrt_7d")
epinow_R_pub <- load_published_R_estimates("epiforecasts")
rtlive_R_pub <- load_published_R_estimates("rtlive")

# merge estimates and plot for comparison
estimates_pub <- RKI_R_pub[,c("date", "R_pub")] %>%
  full_join(ETH_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(Ilmenau_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(SDSC_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(Zi_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(globalrt_R_pub[,c("date", "R_pub")], by = "date") %>%
  full_join(epinow_R_pub[,c("date", "R_pub")], by = "date") %>%
  full_join(rtlive_R_pub[,c("date", "R_pub")], by = "date")

org_methods <- c("RKI", "ETH", "Ilmenau", "SDSC", "Zi",
                 "globalrt", "epiforecasts", "rtlive")
plot_for_comparison(estimates_pub, org_methods, filenames = "_real-time.pdf",
                    method = "original", variation = "parameters")

estimates_pub_ci <- ETH_R_pub %>% 
  full_join(Ilmenau_R_pub, by = "date") %>% 
  full_join(globalrt_R_pub, by = "date") %>%
  full_join(epinow_R_pub, by = "date") #%>%
#full_join(rtlive_R_pub, by = "date")
org_methods <- c("ETH", "Ilmenau", "globalrt", "epiforecasts") #, "rtlive")

plot_for_comparison(estimates_pub_ci, org_methods, include_CI=T,
                    method = "original", variation = "parameters")



################################
# get functions and parameters #
################################

# parameter combinations used in papers
gt_dist <- c("gamma", "constant", "ad hoc", "gamma", "gamma", "gamma", "gamma", "lognorm", "exponential")
mean_gt <- c(4.8,      4,          5.6,      4.8,     5,       3.4,     3.6,     4.7,       7)
sd_gt <-   c(2.3,      0,          4.2,      2.3,     4,       1.8,     3.1,     2.9,       7)
delay <-   c(11,       1,          7,        10,      0,       0,       12,      12,        0)

params <- data.frame(gtd=gt_dist, gt_mean=mean_gt, gt_sd=sd_gt, delay=delay)
methods <- c("ETH", "RKI", "Ilmenau", "SDSC", "Zi", "AGES", "epiforecasts", "rtlive", "globalrt")
rownames(params) <- methods

source("Rt_estimate_reconstruction/ETH/delays_for_ETH_estimation.R")

# save generation time distributions for rtlive
gt_main_analysis <- get_infectivity_profile(gt_type="gamma", gt_mean = 4, gt_sd = 4, n_days = 30)
gt_lower_example <- get_infectivity_profile(gt_type="gamma",
                                            gt_mean = params["AGES", "gt_mean"],
                                            gt_sd = params["AGES", "gt_sd"], n_days = 30)
gt_upper_example <- get_infectivity_profile(gt_type="gamma",
                                            gt_mean = params["Ilmenau", "gt_mean"],
                                            gt_sd = params["Ilmenau", "gt_sd"], n_days = 30)
gt_distributions <- cbind(0:30, gt_main_analysis, gt_lower_example, gt_upper_example)
colnames(gt_distributions) <- c("days_after_transmission", "gamma 4 (4)", "gamma 3.4 (1.8)", "gamma 5.6 (4.2)")

write.csv(gt_distributions, "Rt_estimate_reconstruction/gt_distributions.csv", row.names = F)

plot(gt_main_analysis, type="l", ylim = c(0, 0.32))
lines(gt_lower_example, col="red")
lines(gt_upper_example, col="blue")
#lines(get_infectivity_profile(gt_type="exponential", gt_mean = 4.7, gt_sd = 2.9, n_days = 30), col = "green")
lines(c(0.00000000e+00, 7.33148180e-03, 1.03868478e-01, 1.95147643e-01,
        1.93797118e-01, 1.52828476e-01, 1.09576993e-01, 7.54229927e-02,
        5.10921821e-02, 3.44849062e-02, 2.33407210e-02, 1.58958932e-02,
        1.09121047e-02, 7.55711343e-03, 5.28165515e-03, 3.72529635e-03,
        2.65131301e-03, 1.90352772e-03, 1.37822529e-03, 1.00599866e-03,
        7.40014799e-04, 5.48403039e-04, 4.09289125e-04, 3.07532780e-04,
        2.32567920e-04, 1.76961382e-04, 1.35443284e-04, 1.04249379e-04,
        8.06714362e-05, 6.27474302e-05), col="orange")
legend(x="topright", legend=c("gamma_4", "gamma_3.4", "gamma_5.6", "lnorm_discr_4.7"),
       lty=1, col=c("black", "red", "blue", "orange"))

# plot assumed mean and sd of generation times
gtd_scatterplot <- ggplot(data = params, aes(x = gt_mean, y = gt_sd)) +
  labs(x = "mean", y = "standard deviation") +
  theme_minimal() +
  theme(
    axis.text=element_text(size=16),
    axis.title=element_text(size=18),
    axis.line = element_line(),
    axis.line.y.right = element_line(),
    axis.line.x.top = element_line(),
    legend.position = "bottom",
    panel.grid.minor.x = element_line(colour = "transparent"),
    panel.background = element_rect(fill = "transparent")
    ) +
  geom_point(size = 2) +
  geom_text(label=c("ETH/SDSC", "RKI", "Ilmenau", "SDSC", "Zi", "AGES", "epiforecasts", "rtlive", "globalrt"),
            nudge_y = 0.3, check_overlap = TRUE, size = 5) + 
  geom_point(data = data.frame(gt_mean=4, gt_sd=4), colour = "red", size = 3)
print(gtd_scatterplot)
ggsave(gtd_scatterplot, filename = "gtd_scatterplot.pdf",  bg = "transparent",
       width = 8, height = 5)

# read RKI Nowcast data for RKI estimation
RKI_incid <- load_incidence_data(method = "RKI")

# read smoothed RKI incidence data for SDSC estimation
SDSC_incid <-  read_csv("Rt_estimate_reconstruction/incidence_data/SDSC_incid_21_07_10.csv")

# read incidence data used by rtlive (sourced from RKI line-list data) for other estimations
incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_07_10.csv")

# compare incidence time series used by RKI vs. rtlive
plot(RKI_incid, type="l")
lines(incid, col="blue")
lines(SDSC_incid, col="red")

# choose method for comparison
method <- "globalrt"
params["globalrt", "gt_mean"] <- 4
params["globalrt", "gt_sd"] <- 4

incid_for_ETH <- load_incidence_data(method = "ETHZ_sliding_window", source = "_simpleRKI",
                                     new_deconvolution = FALSE)


##########################
# adjust input data only #
##########################
R_raw_EpiEstim_input <- estimate_RKI_R(incid, method = "EpiEstim",
                                       window = 7,
                                       gt_type = params["RKI", "gtd"],
                                       gt_mean = params["RKI", "gt_mean"],
                                       gt_sd = params["RKI", "gt_sd"])

R_RKI_input <- estimate_RKI_R(RKI_incid, method = "RKI",
                              gt_type = params["RKI", "gtd"],
                              gt_mean = params["RKI", "gt_mean"],
                              gt_sd = params["RKI", "gt_sd"])

R_SDSC_EpiEstim_input <- estimate_SDSC_R(SDSC_incid,
                                         gt_type = params["SDSC", "gtd"],
                                         gt_mean = params["SDSC", "gt_mean"],
                                         gt_sd = params["SDSC", "gt_sd"])

R_ETH_EpiEstim_input <- estimate_ETH_R(incid_for_ETH,
                                       gt_type = params["ETH", "gtd"],
                                       gt_mean = params["ETH", "gt_mean"],
                                       gt_sd = params["ETH", "gt_sd"])

R_AGES_EpiEstim_input <- estimate_AGES_R(incid,
                                         gt_type = params["AGES", "gtd"],
                                         gt_mean=params["AGES", "gt_mean"],
                                         gt_sd=params["AGES", "gt_sd"])

R_Ilmenau_input <- estimate_Ilmenau_R(incid,
                                      gt_type = params["Ilmenau", "gtd"],
                                      gt_mean=params["Ilmenau", "gt_mean"],
                                      gt_sd=params["Ilmenau", "gt_sd"])
names(R_Ilmenau_input) <- c("date", "R_calc", "lower", "upper")

# TODO: repeat estimation with data from 2021-07-10
path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- max(list.files(path, pattern = paste0("R_calc_\\d{4}-\\d{2}-\\d{2}_epiforecasts.qs")))
print(paste("Most recent epiforecasts estimates with matching parameters:", file))
R_epiforecasts_input <- qread(paste0(path, file))
R_epiforecasts_input <- R_epiforecasts_input[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_input) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
file <- paste0("bayesian_smoother_gt7.csv")
R_globalrt_smoother_input <- read_csv(paste0(path, file))
names(R_globalrt_smoother_input) <- c("date", "R_calc", "lower", "upper")

file <- paste0("bayesian_filter_gt7.csv")
R_globalrt_filter_input <- read_csv(paste0(path, file))
names(R_globalrt_filter_input) <- c("date", "R_calc", "lower", "upper")

#plot(R_globalrt_input2$date, R_globalrt_input2$R_calc, type="l")
#lines(R_globalrt_smoother$Date, R_globalrt_smoother$R, col = "red")
#lines(R_globalrt_filter$Date, R_globalrt_filter$R, col = "blue")
#legend(x="topright", legend=c("estimates from KF", "Bayesian smoothed", "Bayesian filtered"),
#       lty=1, col=c("black", "red", "blue"))
#lines(R_globalrt_smoother$Date, R_globalrt_smoother$lb_95, col = "red", lty=2)
#lines(R_globalrt_smoother$Date, R_globalrt_smoother$ub_95, col = "red", lty=2)
#lines(R_globalrt_filter$Date, R_globalrt_filter$lb_95, col = "blue", lty=2)
#lines(R_globalrt_filter$Date, R_globalrt_filter$ub_95, col = "blue", lty=2)

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
file <- "DE_2021-07-10_all_trace_summary.csv"
R_rtlive_input <- read.csv(paste0(path, file)) %>%
  dplyr::select(X, mean, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "mean", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))

# merge estimates and plot for comparison
estimates_input <- R_raw_EpiEstim_input %>%
  full_join(R_RKI_input, by = "date") %>%
  full_join(R_SDSC_EpiEstim_input, by = "date") %>%
  full_join(R_ETH_EpiEstim_input, by = "date") %>% 
  full_join(R_AGES_EpiEstim_input, by = "date") %>% 
  full_join(R_Ilmenau_input, by = "date") %>% 
  full_join(R_epiforecasts_input, by = "date") %>%
  full_join(R_globalrt_smoother_input, by = "date") %>%
  full_join(R_globalrt_filter_input, by = "date") %>%
  full_join(R_rtlive_input, by = "date")

comp_methods <- c("raw EpiEstim", "RKI", "SDSC EpiEstim", "ETH EpiEstim",
                  "AGES EpiEstim", "Ilmenau", "epiforecasts",
                  "globalrt smoother", "globalrt filter", "rtlive")
plot_for_comparison(estimates_input, comp_methods, filenames = "_input-data.pdf",
                    method = "rtlive", variation = "input data")



#######################################
# additionally adjust generation time #
#######################################
R_raw_EpiEstim_gt <- estimate_RKI_R(incid, method = "EpiEstim",
                                    window = 7,
                                    gt_type = "gamma",
                                    gt_mean = params[method, "gt_mean"],
                                    gt_sd = params[method, "gt_sd"])

R_RKI_gt <- estimate_RKI_R(RKI_incid, method = "RKI",
                           gt_type = params["RKI", "gtd"],
                           gt_mean = params[method, "gt_mean"],
                           gt_sd = params[method, "gt_sd"])

R_SDSC_EpiEstim_gt <- estimate_SDSC_R(SDSC_incid,
                                      gt_type = "gamma",
                                      gt_mean = params[method, "gt_mean"],
                                      gt_sd = params[method, "gt_sd"])

R_ETH_EpiEstim_gt <- estimate_ETH_R(incid_for_ETH,
                                    gt_type = "gamma",
                                    gt_mean = params[method, "gt_mean"],
                                    gt_sd = params[method, "gt_sd"])

R_AGES_EpiEstim_gt <- estimate_AGES_R(incid,
                                      gt_type = "gamma",
                                      gt_mean=params[method, "gt_mean"],
                                      gt_sd=params[method, "gt_sd"])

R_Ilmenau_gt <- estimate_Ilmenau_R(incid,
                                   gt_type = "gamma",
                                   gt_mean=params[method, "gt_mean"],
                                   gt_sd=params[method, "gt_sd"])
names(R_Ilmenau_gt) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- max(list.files(path, pattern = paste0("R_calc_\\d{4}-\\d{2}-\\d{2}_", method, "_GTD.qs")))
print(paste("Most recent epiforecasts estimates with matching parameters:", file))
R_epiforecasts_gt <- qread(paste0(path, file))
R_epiforecasts_gt <- R_epiforecasts_gt[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_gt) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
file <- paste0("bayesian_smoother_gt4.csv")
R_globalrt_smoother_gt <- read_csv(paste0(path, file))
names(R_globalrt_smoother_gt) <- c("date", "R_calc", "lower", "upper")

file <- paste0("bayesian_filter_gt4.csv")
R_globalrt_filter_gt <- read_csv(paste0(path, file))
names(R_globalrt_filter_gt) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
file <- "gamma_4_de_2021-07-10_all_summary.csv"
R_rtlive_gt <- read.csv(paste0(path, file)) %>%
  dplyr::select(X, mean, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "mean", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))

# merge estimates and plot for comparison
estimates_GTD <- R_raw_EpiEstim_gt %>%
  full_join(R_RKI_gt, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_gt, by = "date") %>% 
  full_join(R_ETH_EpiEstim_gt, by = "date") %>% 
  full_join(R_AGES_EpiEstim_gt, by = "date") %>% 
  full_join(R_Ilmenau_gt, by = "date") %>% 
  full_join(R_epiforecasts_gt, by = "date") %>%
  full_join(R_globalrt_smoother_gt, by = "date") %>%
  full_join(R_globalrt_filter_gt, by = "date") %>%
  full_join(R_rtlive_gt, by = "date")
  #full_join(R_rtlive_input, by = "date")

comp_methods <- c("raw EpiEstim", "RKI", "SDSC EpiEstim", "ETH EpiEstim",
                  "AGES EpiEstim", "Ilmenau", "epiforecasts",
                  "globalrt smoother", "globalrt filter", "rtlive")
plot_for_comparison(estimates_GTD, comp_methods, filenames = "_gtd.pdf",
                    method = method, variation = "generation time")

estimates_GTD_CI <- as.data.frame(R_SDSC_EpiEstim_gt) %>% 
  full_join(R_ETH_EpiEstim_gt, by = "date") %>% 
  full_join(R_Ilmenau_gt, by = "date") %>% 
  #full_join(R_epiforecasts_gt, by = "date") %>%
  full_join(R_globalrt_gt, by = "date") %>%
  full_join(R_rtlive_gt, by = "date")

comp_CI <- c("SDSC EpiEstim", "ETH EpiEstim", "Ilmenau", #"epiforecasts",
             "globalrt", "rtlive")
plot_for_comparison(estimates_GTD_CI, comp_CI, include_CI = T,
                    method = method, variation = "generation time")

dataset_globalrt_org <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/dataset.csv")
dataset_globalrt_org <- dataset_globalrt_org %>% dplyr::filter(`Country/Region`=="Germany")
plot(dataset_globalrt$Date, dataset_globalrt$new_cases, type="l")
lines(dataset_globalrt_org$Date, dataset_globalrt_org $new_cases, col="red")
plot(dataset_globalrt$Date, dataset_globalrt$gr_infected_4, type="l")
lines(dataset_globalrt_org$Date, dataset_globalrt_org$gr_infected_4, col="red")






##############################
# additionally adjust delays #
##############################

# for comparison of methods use default window size of each method
R_raw_EpiEstim_d <- R_raw_EpiEstim_gt
R_raw_EpiEstim_d$date <- R_raw_EpiEstim_d$date + params["RKI", "delay"]

R_RKI_d <- R_RKI_gt
R_RKI_d$date <- R_RKI_d$date + params["RKI", "delay"]
# additionally correct for mean difference between reference date and reporting date (Nowcast vs. rtlive aggregation used for other forecasts)
R_RKI_d$date <- R_RKI_d$date + 3 # eigentlich ca 3.37 Tage

R_SDSC_EpiEstim_d <- R_SDSC_EpiEstim_gt
R_SDSC_EpiEstim_d$date <- R_SDSC_EpiEstim_d$date + params["SDSC", "delay"]

R_ETH_EpiEstim_d <- R_ETH_EpiEstim_gt
R_ETH_EpiEstim_d$date <- R_ETH_EpiEstim_d$date + params["ETH", "delay"]

R_AGES_EpiEstim_d <- R_AGES_EpiEstim_gt
R_AGES_EpiEstim_d$date <- R_AGES_EpiEstim_d$date + params["AGES", "delay"]

R_Ilmenau_d <- R_Ilmenau_gt
R_Ilmenau_d$date <- R_Ilmenau_d$date + params["Ilmenau", "delay"]

R_epiforecasts_d <- R_epiforecasts_gt
R_epiforecasts_d$date <- R_epiforecasts_d$date + params["epiforecasts", "delay"]

R_globalrt_smoother_d <- R_globalrt_smoother_gt
R_globalrt_smoother_d$date <- R_globalrt_smoother_d$date + params["globalrt", "delay"]

R_globalrt_filter_d <- R_globalrt_filter_gt
R_globalrt_filter_d$date <- R_globalrt_filter_d$date + params["globalrt", "delay"]

R_rtlive_d <- R_rtlive_gt
#R_rtlive_d <- R_rtlive_input
R_rtlive_d$date <- R_rtlive_d$date + params["rtlive", "delay"]

# merge estimates and plot for comparison
estimates_delays <- R_raw_EpiEstim_d %>%
  full_join(R_RKI_d, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_d, by = "date") %>% 
  full_join(R_ETH_EpiEstim_d, by = "date") %>% 
  full_join(R_AGES_EpiEstim_d, by = "date") %>% 
  full_join(R_Ilmenau_d, by = "date") %>% 
  full_join(R_epiforecasts_d, by = "date") %>%
  full_join(R_globalrt_smoother_d, by = "date") %>%
  full_join(R_globalrt_filter_d, by = "date") %>%
  full_join(R_rtlive_d, by = "date")

comp_methods <- c("raw EpiEstim", "RKI", "SDSC EpiEstim", "ETH EpiEstim",
                  "AGES EpiEstim", "Ilmenau", "epiforecasts",
                  "globalrt smoother", "globalrt filter", "rtlive")
plot_for_comparison(estimates_delays, comp_methods, filenames = "_delays.pdf",
                    method = method, variation = "delays")

estimates_delays_CI <- as.data.frame(R_ETH_EpiEstim_d) %>% 
  full_join(R_ETH_EpiEstim_d, by = "date") %>% 
  full_join(R_Ilmenau_d, by = "date") %>% 
  #full_join(R_epiforecasts_d, by = "date") %>%
  #full_join(R_globalrt_d, by = "date") %>%
  full_join(R_rtlive_d, by = "date")

comp_CI <- c("SDSC EpiEstim", "ETH EpiEstim", "Ilmenau", #"epiforecasts", "globalrt",
             "rtlive")
plot_for_comparison(estimates_delays_CI, comp_CI, include_CI = T,
                    method = method, variation = "delays")

###################################
# additionally adjust window size #
###################################
R_raw_EpiEstim_ws <- R_raw_EpiEstim_d

R_RKI_ws <- R_RKI_d

R_SDSC_EpiEstim_ws <- estimate_SDSC_R(SDSC_incid,
                                      window = 7,
                                      gt_type = "gamma",
                                      gt_mean = params[method, "gt_mean"],
                                      gt_sd = params[method, "gt_sd"])
R_SDSC_EpiEstim_ws$date <- R_SDSC_EpiEstim_ws$date + params["SDSC", "delay"]

R_ETH_EpiEstim_ws <- estimate_ETH_R(incid_for_ETH,
                                    window = 7,
                                    gt_type = "gamma",
                                    gt_mean = params[method, "gt_mean"],
                                    gt_sd = params[method, "gt_sd"])
R_ETH_EpiEstim_ws$date <- R_ETH_EpiEstim_ws$date + params["ETH", "delay"]

R_AGES_EpiEstim_ws <- estimate_AGES_R(incid,
                                      window = 7,
                                      gt_type = "gamma",
                                      gt_mean=params[method, "gt_mean"],
                                      gt_sd=params[method, "gt_sd"],
                                      delay = params["AGES", "delay"])
R_AGES_EpiEstim_ws$date <- R_AGES_EpiEstim_ws$date + params["AGES", "delay"]

R_Ilmenau_ws <- estimate_Ilmenau_R(incid,
                                   window = 7,
                                   gt_type = "gamma",
                                   gt_mean=params[method, "gt_mean"],
                                   gt_sd=params[method, "gt_sd"],
                                   delay = params["Ilmenau", "delay"])
names(R_Ilmenau_ws) <- c("date", "R_calc", "lower", "upper")
R_Ilmenau_ws$date <- R_Ilmenau_ws$date + params["Ilmenau", "delay"]

R_epiforecasts_ws <- R_epiforecasts_d

# without substitute for window size = 7
R_globalrt_smoother_ws <- R_globalrt_smoother_d

R_globalrt_filter_ws <- R_globalrt_filter_d

R_rtlive_ws <- R_rtlive_d

# merge estimates and plot for comparison
estimates_window_Cori <- R_raw_EpiEstim_ws %>%
  full_join(R_RKI_ws, by = "date") %>% 
  full_join(R_SDSC_EpiEstim_ws, by = "date") %>% 
  full_join(R_ETH_EpiEstim_ws, by = "date") %>% 
  full_join(R_AGES_EpiEstim_ws, by = "date") %>% 
  full_join(R_Ilmenau_ws, by = "date") %>% 
  full_join(R_epiforecasts_ws, by = "date") %>%
  full_join(R_globalrt_smoother_ws, by = "date") %>%
  full_join(R_globalrt_filter_ws, by = "date") %>%
  full_join(R_rtlive_ws, by = "date")

comp_methods <- c("raw EpiEstim", "RKI", "SDSC EpiEstim", "ETH EpiEstim",
                  "AGES EpiEstim", "Ilmenau", "epiforecasts",
                  "globalrt smoother", "globalrt filter", "rtlive")
plot_for_comparison(estimates_window_Cori, comp_methods, filenames = "_windowCori.pdf",
                    method = method, variation = "generation time")

estimates_ws_CI <- as.data.frame(R_raw_EpiEstim_ws) %>% 
  full_join(R_ETH_EpiEstim_ws, by = "date") %>%
  full_join(R_Ilmenau_ws, by = "date") %>% 
  full_join(R_epiforecasts_ws, by = "date") %>%
  full_join(R_globalrt_ws, by = "date") %>%
  full_join(R_rtlive_ws, by = "date")

comp_CI <- c("raw EpiEstim", "ETH EpiEstim", "Ilmenau", "epiforecasts", "globalrt", "rtlive")
plot_for_comparison(estimates_ws_CI, comp_CI, include_CI = T,
                    method = method, variation = "delays")



# with preprocessing of incidence data for globalrt method
path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
file <- paste0("estimated_R_", method, "_window_data_from_rtlive.csv")
R_globalrt_ws <- read_csv(paste0(path, file))
R_globalrt_ws <- R_globalrt_ws[R_globalrt_ws$`Country/Region` == "Germany", c("Date", "R", "ci_95_l", "ci_95_u")]
names(R_globalrt_ws) <- c("date", "R_calc", "lower", "upper")

# merge estimates and plot for comparison
estimates_window <- R_raw_EpiEstim_ws %>%
  full_join(R_ETH_EpiEstim_ws, by = "date") %>%
  full_join(R_AGES_EpiEstim_ws, by = "date") %>%
  full_join(R_Ilmenau_ws, by = "date") %>%
  full_join(R_epiforecasts_ws, by = "date") %>%
  full_join(R_globalrt_ws, by = "date") %>%
  full_join(R_rtlive_ws, by = "date")

comp_methods <- c("raw EpiEstim", "ETH EpiEstim", "AGES EpiEstim", "Ilmenau",
                  "epiforecasts", "globalrt", "rtlive")
plot_for_comparison(estimates_window, comp_methods, filenames = "_windowAll.pdf",
                    method = method, variation = "generation time")





# for comparison of methods use default window size of each method
R_raw_EpiEstim <- estimate_RKI_R(incid, method = "EpiEstim",
                                 window = 7,
                                 gt_type = params[method, "gtd"],
                                 gt_mean = params[method, "gt_mean"],
                                 gt_sd = params[method, "gt_sd"],
                                 delay = params[method, "delay"])

incid_ma7 <- incid
incid_ma7$I <- as.numeric(stats::filter(incid, rep(1/7, 7), method = "convolution", sides = 1)[,2])

incid_ms7 <- incid
incid_ms7$I <- as.numeric(stats::filter(incid, rep(7, 7), method = "convolution", sides = 1)[,2])

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
  full_join(R_epiforecasts, by = "date") %>%
  full_join(R_globalrt, by = "date")

comp_methods <- c("raw EpiEstim", "ETH EpiEstim", "AGES EpiEstim", "Ilmenau", "epiforecasts",
                  "globalrt")
plot_for_comparison(estimates_allPars, comp_methods, method = method, variation = "parameters")

estimates_allPars_CI <- as.data.frame(R_ETH_EpiEstim) %>% 
  full_join(R_Ilmenau, by = "date") %>% 
  full_join(R_epiforecasts, by = "date") %>%
  full_join(R_globalrt, by = "date")

comp_CI <- c("ETH EpiEstim", "Ilmenau", "epiforecasts",
             "globalrt")
plot_for_comparison(estimates_allPars_CI, comp_CI, include_CI = T,
                    method = method, variation = "parameters")


















