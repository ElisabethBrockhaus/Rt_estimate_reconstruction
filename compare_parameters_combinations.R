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

# parameter combinations used in papers
gt_dist <- c("gamma", "constant", "ad hoc", "gamma", "gamma", "gamma", "gamma", "lognorm", "exponential")
mean_gt <- c(4.8,      4,          5.6,      4.8,     5,       3.4,     3.6,     4.7,       7)
sd_gt <-   c(2.3,      0,          4.2,      2.3,     4,       1.8,     3.1,     2.9,       7)
delay <-   c(11,       1,          7,        10,      0,       0,       12,      12,        0)

params <- data.frame(gtd=gt_dist, gt_mean=mean_gt, gt_sd=sd_gt, delay=delay)
methods <- c("ETH", "RKI", "Ilmenau", "SDSC", "Zi", "AGES", "epiforecasts", "rtlive", "globalrt")
rownames(params) <- methods

# save generation time distributions for rtlive
gt_main_analysis <- get_infectivity_profile(gt_type="gamma", gt_mean = 4, gt_sd = 4)
gt_lower_example <- get_infectivity_profile(gt_type="gamma",
                                            gt_mean = params["AGES", "gt_mean"],
                                            gt_sd = params["AGES", "gt_sd"])
gt_upper_example <- get_infectivity_profile(gt_type="gamma",
                                            gt_mean = params["Ilmenau", "gt_mean"],
                                            gt_sd = params["Ilmenau", "gt_sd"])
gt_distributions <- cbind(0:1000, gt_main_analysis, gt_lower_example, gt_upper_example)
colnames(gt_distributions) <- c("days_after_transmission", "gamma 4 (4)", "gamma 3.4 (1.8)", "gamma 5.6 (4.2)")

write.csv(gt_distributions, "Rt_estimate_reconstruction/gt_distributions.csv", row.names = F)

plot(gt_main_analysis[1:30], type="l", ylim = c(0, 0.32))
lines(gt_lower_example, col="red")
lines(gt_upper_example, col="blue")

R_Cori_exponential <- estimate_RKI_R(incid, method = "EpiEstim",
                                     window = 7,
                                     gt_type = "exponential",
                                     gt_mean = 4,
                                     gt_sd = 4)

R_Cori_gamma <- estimate_RKI_R(incid, method = "EpiEstim",
                               window = 7,
                               gt_type = "gamma",
                               gt_mean = 4,
                               gt_sd = 4)

plot(R_Cori_exponential, type="l")


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
incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid.csv")

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

path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- max(list.files(path, pattern = paste0("R_calc_\\d{4}-\\d{2}-\\d{2}_epiforecasts.qs")))
print(paste("Most recent epiforecasts estimates with matching parameters:", file))
R_epiforecasts_input <- qread(paste0(path, file))
R_epiforecasts_input <- R_epiforecasts_input[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_input) <- c("date", "R_calc", "lower", "upper")

#path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
#file <- paste0("estimated_R_globalrt.csv")
#R_globalrt_input <- read_csv(paste0(path, file))
#R_globalrt_input <- R_globalrt_input[R_globalrt_input$`Country/Region` == "Germany", c("Date", "R", "ci_95_l", "ci_95_u")]
#names(R_globalrt_input) <- c("date", "R_calc", "lower", "upper")

file <- paste0("estimated_R_globalrt_data_from_rtlive.csv")
R_globalrt_input2 <- read_csv(paste0(path, file))
R_globalrt_input2 <- R_globalrt_input2[R_globalrt_input2$`Country/Region` == "Germany", c("Date", "R", "ci_95_l", "ci_95_u")]
names(R_globalrt_input2) <- c("date", "R_calc", "lower", "upper")

file <- paste0("estimated_R_globalrt_data_from.csv")
R_globalrt_input3 <- read_csv(paste0(path, file))
R_globalrt_input3 <- R_globalrt_input3[R_globalrt_input3$`Country/Region` == "Germany", c("Date", "R", "ci_95_l", "ci_95_u")]
names(R_globalrt_input3) <- c("date", "R_calc", "lower", "upper")

plot(R_globalrt_input3$date, R_globalrt_input3$R_calc, type="l")
lines(R_globalrt_input2$date, R_globalrt_input2$R_calc, col = "red")
legend(x="topright", legend=c("estimates from original data", "estimates from rtlive data"),
       lty=c(1,1), col=c("black", "red"))

# merge estimates and plot for comparison
estimates_input <- R_raw_EpiEstim_input %>%
  full_join(R_ETH_EpiEstim_input, by = "date") %>% 
  full_join(R_AGES_EpiEstim_input, by = "date") %>% 
  full_join(R_Ilmenau_input, by = "date") %>% 
  full_join(R_epiforecasts_input, by = "date") %>%
  full_join(R_globalrt_input2, by = "date")

comp_methods <- c("raw EpiEstim", "ETH EpiEstim", "AGES EpiEstim",
                  "Ilmenau", "epiforecasts", "globalrt")
plot_for_comparison(estimates_input, comp_methods, filenames = "_input-data.pdf",
                    method = "rtlive", variation = "input data")



#######################################
# additionally adjust generation time #
#######################################
R_raw_EpiEstim_gt <- estimate_RKI_R(incid, method = "EpiEstim",
                                    window = 7,
                                    gt_type = params[method, "gtd"],
                                    gt_mean = params[method, "gt_mean"],
                                    gt_sd = params[method, "gt_sd"])

R_ETH_EpiEstim_gt <- estimate_ETH_R(incid_for_ETH,
                                    gt_type = params[method, "gtd"],
                                    gt_mean = params[method, "gt_mean"],
                                    gt_sd = params[method, "gt_sd"])

R_AGES_EpiEstim_gt <- estimate_AGES_R(incid,
                                      gt_type = params[method, "gtd"],
                                      gt_mean=params[method, "gt_mean"],
                                      gt_sd=params[method, "gt_sd"])

R_Ilmenau_gt <- estimate_Ilmenau_R(incid,
                                   gt_type = params[method, "gtd"],
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
file <- paste0("estimated_R_", method, "_GTD_data_from_rtlive.csv")
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
plot_for_comparison(estimates_GTD, comp_methods, filenames = "_gtd.pdf",
                    method = method, variation = "generation time")

estimates_GTD_CI <- as.data.frame(R_ETH_EpiEstim_gt) %>% 
  full_join(R_Ilmenau_gt, by = "date") %>% 
  full_join(R_epiforecasts_gt, by = "date") %>%
  full_join(R_globalrt_gt, by = "date")

comp_CI <- c("ETH EpiEstim", "Ilmenau", "epiforecasts", "globalrt")
plot_for_comparison(estimates_GTD_CI, comp_CI, include_CI = T,
                    method = method, variation = "generation time")

##############################
# additionally adjust delays #
##############################

# for comparison of methods use default window size of each method
R_raw_EpiEstim_d <- R_raw_EpiEstim_gt
R_raw_EpiEstim_d$date <- R_raw_EpiEstim_d$date + params["RKI", "delay"]

R_ETH_EpiEstim_d <- R_ETH_EpiEstim_gt
R_ETH_EpiEstim_d$date <- R_ETH_EpiEstim_d$date + params["ETH", "delay"]

R_AGES_EpiEstim_d <- R_AGES_EpiEstim_gt
R_AGES_EpiEstim_d$date <- R_AGES_EpiEstim_d$date + params["AGES", "delay"]

R_Ilmenau_d <- R_Ilmenau_gt
R_Ilmenau_d$date <- R_Ilmenau_d$date + params["Ilmenau", "delay"]

R_epiforecasts_d <- R_epiforecasts_gt
R_epiforecasts_d$date <- R_epiforecasts_d$date + params["epiforecasts", "delay"]

R_globalrt_d <- R_globalrt_gt
R_globalrt_d$date <- R_globalrt_d$date + params["globalrt", "delay"]

# merge estimates and plot for comparison
estimates_delays <- R_raw_EpiEstim_d %>%
  full_join(R_ETH_EpiEstim_d, by = "date") %>% 
  full_join(R_AGES_EpiEstim_d, by = "date") %>% 
  full_join(R_Ilmenau_d, by = "date") %>% 
  full_join(R_epiforecasts_d, by = "date") %>%
  full_join(R_globalrt_d, by = "date")

comp_methods <- c("raw EpiEstim", "ETH EpiEstim", "AGES EpiEstim", "Ilmenau", "epiforecasts", "globalrt")
plot_for_comparison(estimates_delays, comp_methods, filenames = "_delays.pdf",
                    method = method, variation = "delays")

estimates_delays_CI <- as.data.frame(R_ETH_EpiEstim_d) %>% 
  full_join(R_Ilmenau_d, by = "date") %>% 
  full_join(R_epiforecasts_d, by = "date") %>%
  full_join(R_globalrt_d, by = "date")

comp_CI <- c("ETH EpiEstim", "Ilmenau", "epiforecasts", "globalrt")
plot_for_comparison(estimates_delays_CI, comp_CI, include_CI = T,
                    method = method, variation = "delays")

###################################
# additionally adjust window size #
###################################
R_raw_EpiEstim_ws <- R_raw_EpiEstim_d

R_ETH_EpiEstim_ws <- estimate_ETH_R(incid_for_ETH,
                                    window = 7,
                                    gt_type = params[method, "gtd"],
                                    gt_mean = params[method, "gt_mean"],
                                    gt_sd = params[method, "gt_sd"])
R_ETH_EpiEstim_ws$date <- R_ETH_EpiEstim_ws$date + params["ETH", "delay"]

R_AGES_EpiEstim_ws <- estimate_AGES_R(incid,
                                      window = 7,
                                      gt_type = params[method, "gtd"],
                                      gt_mean=params[method, "gt_mean"],
                                      gt_sd=params[method, "gt_sd"],
                                      delay = params["AGES", "delay"])
R_AGES_EpiEstim_ws$date <- R_AGES_EpiEstim_ws$date + params["AGES", "delay"]

R_Ilmenau_ws <- estimate_Ilmenau_R(incid,
                                   window = 7,
                                   gt_type = params[method, "gtd"],
                                   gt_mean=params[method, "gt_mean"],
                                   gt_sd=params[method, "gt_sd"],
                                   delay = params["Ilmenau", "delay"])
names(R_Ilmenau_ws) <- c("date", "R_calc", "lower", "upper")
R_Ilmenau_ws$date <- R_Ilmenau_ws$date + params["Ilmenau", "delay"]

R_epiforecasts_ws <- R_epiforecasts_d

# without substitute for window size = 7
R_globalrt_ws <- R_globalrt_d

# merge estimates and plot for comparison
estimates_window_Cori <- R_raw_EpiEstim_ws %>%
  full_join(R_ETH_EpiEstim_ws, by = "date") %>% 
  full_join(R_AGES_EpiEstim_ws, by = "date") %>% 
  full_join(R_Ilmenau_ws, by = "date") %>% 
  full_join(R_epiforecasts_ws, by = "date") %>%
  full_join(R_globalrt_ws, by = "date")

comp_methods <- c("raw EpiEstim", "ETH EpiEstim", "AGES EpiEstim", "Ilmenau", "epiforecasts", "globalrt")
plot_for_comparison(estimates_window_Cori, comp_methods, filenames = "_windowCori.pdf",
                    method = method, variation = "generation time")


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
  full_join(R_globalrt_ws, by = "date")

comp_methods <- c("raw EpiEstim", "ETH EpiEstim", "AGES EpiEstim", "Ilmenau", "epiforecasts", "globalrt")
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


















