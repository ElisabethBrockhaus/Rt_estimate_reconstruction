library(readr)

setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()


source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")


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
epinow_R_pub <- load_published_R_estimates("epiforecasts")
rtlive_R_pub <- load_published_R_estimates("rtlive")

# merge estimates and plot for comparison
estimates_pub <- RKI_R_pub[,c("date", "R_pub")] %>%
  full_join(ETH_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(Ilmenau_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(SDSC_R_pub[,c("date", "R_pub")], by = "date") %>% 
  #full_join(Zi_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(globalrt_R_pub[,c("date", "R_pub")], by = "date") %>%
  full_join(epinow_R_pub[,c("date", "R_pub")], by = "date") %>%
  full_join(rtlive_R_pub[,c("date", "R_pub")], by = "date")

org_methods <- c("RKI", "ETH", "Ilmenau", "SDSC", #"Zi",
                 "globalrt", "epiforecasts", "rtlive")
plot_for_comparison(estimates_pub, org_methods,
                    legend_name = "research group", filenames = "_real-time.pdf",
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
gt_dist <- c("gamma", "constant", "ad hoc", "gamma", "gamma", "gamma", "lognorm", "exponential")
mean_gt <- c(4.8,      4,          5.6,      4.8,    3.4,     3.6,     4.7,       7)
sd_gt <-   c(2.3,      0,          4.2,      2.3,    1.8,     3.1,     2.9,       7)
delay <-   c(11,       1,          7,        10,     0,       12,      12,        0)

params <- data.frame(gtd=gt_dist, gt_mean=mean_gt, gt_sd=sd_gt, delay=delay)
methods <- c("ETH", "RKI", "Ilmenau", "SDSC", "AGES", "epiforecasts", "rtlive", "globalrt")
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
  geom_text(label=c("ETH/SDSC", "RKI", "Ilmenau", "SDSC", "AGES", "epiforecasts", "rtlive", "globalrt"),
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



#####################
# adjust everything #
#####################
window <- 7
gt_type <- "gamma"
gt_mean <- 4
gt_sd <- 4

R_raw_EpiEstim_adjAll <- estimate_RKI_R(incid, method = "EpiEstim",
                                        window = window,
                                        gt_type = gt_type,
                                        gt_mean = gt_mean,
                                        gt_sd = gt_sd)
R_raw_EpiEstim_adjAll$date <- R_raw_EpiEstim_adjAll$date + params["RKI", "delay"]

# there is no gt_type and gt_sd in the original RKI method
R_RKI_adjAll <- estimate_RKI_R(RKI_incid, method = "RKI",
                               window = window,
                               gt_type = "constant",
                               gt_mean = gt_mean,
                               gt_sd = 0)
R_RKI_adjAll$date <- R_RKI_adjAll$date + params["RKI", "delay"] + 3

R_SDSC_EpiEstim_adjAll <- estimate_SDSC_R(SDSC_incid,
                                          window = window,
                                          gt_type = gt_type,
                                          gt_mean = gt_mean,
                                          gt_sd = gt_sd)
R_SDSC_EpiEstim_adjAll$date <- R_SDSC_EpiEstim_adjAll$date + params["SDSC", "delay"]

R_ETH_EpiEstim_adjAll <- estimate_ETH_R(incid_for_ETH,
                                        window = window,
                                        gt_type = gt_type,
                                        gt_mean = gt_mean,
                                        gt_sd = gt_sd)
R_ETH_EpiEstim_adjAll$date <- R_ETH_EpiEstim_adjAll$date + params["ETH", "delay"]

R_AGES_EpiEstim_adjAll <- estimate_AGES_R(incid,
                                          window = window,
                                          gt_type = gt_type,
                                          gt_mean = gt_mean,
                                          gt_sd = gt_sd,
                                          delay = params["AGES", "delay"])
R_AGES_EpiEstim_adjAll$date <- R_AGES_EpiEstim_adjAll$date + params["AGES", "delay"]

R_Ilmenau_adjAll <- estimate_Ilmenau_R(incid,
                                       window = window,
                                       gt_type = gt_type,
                                       gt_mean = gt_mean,
                                       gt_sd = gt_sd,
                                       delay = params["Ilmenau", "delay"])
names(R_Ilmenau_adjAll) <- c("date", "R_calc", "lower", "upper")
R_Ilmenau_adjAll$date <- R_Ilmenau_adjAll$date + params["Ilmenau", "delay"]

# TODO: download from server when estimation done
path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- "R_calc_2021-07-10_final.qs"
R_epiforecasts_adjAll <- qread(paste0(path, file))
R_epiforecasts_adjAll <- R_epiforecasts_adjAll[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_adjAll) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
R_globalrt_smoother_adjAll <- read_csv(paste0(path, "bayesian_smoother_4.csv"))
names(R_globalrt_smoother_adjAll) <- c("date", "R_calc", "lower", "upper")
R_globalrt_filter_adjAll <- read_csv(paste0(path, "bayesian_filter_4.csv"))
names(R_globalrt_filter_adjAll) <- c("date", "R_calc", "lower", "upper")

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
file <- "gamma_4_de_2021-07-10_all_summary.csv"
R_rtlive_adjAll <- read.csv(paste0(path, file)) %>%
  dplyr::select(X, mean, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "mean", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))
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
plot_for_comparison(estimates_adjAll, comp_methods, filenames = "_adjAll.pdf",
                    method = method, variation = "everything")

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
plot_for_comparison(estimates_adjAll_CI, comp_CI, include_CI = T, filenames = "_CI_adjAll.pdf",
                    method = method, variation = "everything")



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
plot_for_comparison(estimates_notAdjGTD, comp_methods, filenames = "_notAdjGTD.pdf",
                    method = method, variation = "not generation time")

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
plot_for_comparison(estimates_notAdjGTD_CI, comp_CI, include_CI = T, filenames = "_CI_notAdjGTD.pdf",
                    method = method, variation = "not generation time")




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
plot_for_comparison(estimates_notAdjDelays, comp_methods, filenames = "_notAdjDelays.pdf",
                    method = method, variation = "not delays")

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
plot_for_comparison(estimates_notAdjDelays_CI, comp_CI, include_CI = T, filenames = "_CI_notAdjDelays.pdf",
                    method = method, variation = "not delays")



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
plot_for_comparison(estimates_notAdjWindow_CI, comp_CI, include_CI = T, filenames = "_CI_notAdjWindow.pdf",
                    method = method, variation = "not window size")


###############################################################################
# old variations

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
file <- paste0("bayesian_smoother7.csv")
R_globalrt_smoother_input <- read_csv(paste0(path, file))
names(R_globalrt_smoother_input) <- c("date", "R_calc", "lower", "upper")

file <- paste0("bayesian_filter7.csv")
R_globalrt_filter_input <- read_csv(paste0(path, file))
names(R_globalrt_filter_input) <- c("date", "R_calc", "lower", "upper")

plot(R_globalrt_smoother_input$date, R_globalrt_smoother_input$R_calc, type="l")
lines(R_globalrt_filter_input$date, R_globalrt_filter_input$R_calc, col = "blue")
legend(x="topright", legend=c("Bayesian smoother", "Bayesian filter"),
       lty=1, col=c("black", "blue"))

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

