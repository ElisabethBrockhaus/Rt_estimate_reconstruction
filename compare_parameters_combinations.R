library(pheatmap)
library(viridis)

setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()

# look at uncertainty of incidence time series
data <- read_csv("https://raw.githubusercontent.com/robert-koch-institut/SARS-CoV-2-Nowcasting_und_-R-Schaetzung/main/Nowcast_R_aktuell.csv")
plot(data$Datum, data$PS_COVID_Faelle, type="l", xlab="date", ylab="incidence",
     xlim = c(data$Datum[551], data$Datum[583]), ylim = c(5000, 15000))
axis(1, at = seq(data$Datum[551], data$Datum[583], by="weeks"),
     labels = seq(data$Datum[551], data$Datum[583], by="weeks"))
lines(data$Datum, data$UG_PI_COVID_Faelle, col = "blue")
lines(data$Datum, data$OG_PI_COVID_Faelle, col = "blue")
abline(v = as.Date("2021-10-01"), col = "red")


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

# use RKI data for all methods (not line list!)
incid <- load_incidence_data(method = "RKI")
incid <- incid[incid$date < "2021-10-01",]

rtlive_incid <- read_csv("Rt_estimate_reconstruction/rtlive/rtlive-global/data/rtlive_data.csv")
rtlive_region_incid <- rtlive_incid[rtlive_incid$region!="all",]
rtlive_incid_agg <- aggregate(rtlive_region_incid$new_cases, by = list(rtlive_region_incid$date), FUN = sum)
names(rtlive_incid_agg) <- c("date", "I")
plot(incid, type="l")
lines(rtlive_incid_agg, col="red")
lines(rtlive_incid[rtlive_incid$region=="all", c("date", "new_cases")], col="blue")

# save incidence data for epiforecast estimation
write_csv(incid, "Rt_estimate_reconstruction/incidence_data/RKI_incid.csv")

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

path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- max(list.files(path, pattern = paste0("R_calc_\\d{4}-\\d{2}-\\d{2}", method, "Params.qs")))
print(paste("Most recent epiforecasts estimates with matching parameters:", file))
R_epiforecasts <- qread(paste0(path, file))
R_epiforecasts <- R_epiforecasts[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts) <- c("date", "R_calc", "0.025", "0.975")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
file <- max(list.files(path, pattern = paste0("estimated_R_\\d{4}-\\d{2}-\\d{2}_", method, ".csv")))
print(paste("Most recent globalrt estimates with matching parameters:", file))
R_globalrt <- read_csv(paste0(path, file))
R_globalrt <- R_globalrt[R_globalrt$`Country/Region` == "Germany", c("Date", "R")]
names(R_globalrt) <- c("date", "R_calc")

# merge estimates and plot for comparison
estimates <- R_raw_EpiEstim[,c("date", "R_calc")] %>%
  full_join(R_ETH_EpiEstim[,c("date", "R_calc")], by = "date") %>% 
  full_join(R_AGES_EpiEstim[,c("date", "R_calc")], by = "date") %>% 
  full_join(R_Ilmenau[,c("date", "R_calc")], by = "date") %>% 
  full_join(R_epiforecasts[,c("date", "R_calc")], by = "date") %>%
  full_join(R_globalrt[,c("date", "R_calc")], by = "date")
names(estimates) <- c("date", "rawEpiEstim", "ETH", "AGES", "Ilmenau", "epiforecasts", "globalrt")

plot_multiple_estimates(estimates[estimates$date > "2021-04-12" & estimates$date < "2021-10-01",],
                        methods = c("raw EpiEstim", "ETH EpiEstim", "AGES EpiEstim",
                                    "Ilmenau", "epiforecasts",
                                    "globalrt"))

latest_estimates <- estimates[rowSums(is.na(estimates)) == 0,]
latest_estimates <- latest_estimates[latest_estimates$date > "2021-04-12",]

n <- dim(latest_estimates)[2] - 1
comp_methods <- names(latest_estimates)[2:(n+1)]
matr <- matrix(rep(rep(0,n), n), ncol=n)
corr <- matrix(rep(rep(0,n), n), ncol=n)
colnames(matr) <- rownames(matr) <- colnames(corr) <- rownames(corr) <- comp_methods

par(mfrow=c(6,6))
for (method1 in comp_methods) {
  for (method2 in comp_methods){
    diff <- latest_estimates[,method1] - latest_estimates[,method2]
    if (plot_all_differences){
      plot(diff, type="l", ylab = "diff", ylim=c(-0.5, 0.5),
           main = paste(method1, "vs.", method2))
    }
    matr[method1, method2] <- mean(abs(diff))
    corr[method1, method2] <- cor(latest_estimates[,method1], latest_estimates[,method2])
  }
}
par(mfrow=c(1,1))

pheatmap(matr, color = viridis(100), breaks = seq(0,0.4,0.4/100),
         border_color = NA, display_numbers = TRUE,
         fontsize = 12, fontsize_number=20, number_color = "white",
         angle_col = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
         main = paste("Mean absolute differences between estimates over",
                      max(latest_estimates$date) - min(latest_estimates$date),
                      "days using", method, "parameters"))

pheatmap(corr, color = viridis(100, direction = -1), breaks = seq(0.5,1,0.5/100),
         border_color = NA, display_numbers = TRUE,
         fontsize = 12, fontsize_number=20, number_color = "white",
         angle_col = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
         main = paste("Correlations between estimates over",
                      max(latest_estimates$date) - min(latest_estimates$date),
                      "days using", method, "parameters"))


####################################
# in contrast: real-time estimates #
####################################
RKI_R_pub <- load_published_R_estimates(source = "RKI_7day")
ETH_R_pub <- load_published_R_estimates("ETHZ_sliding_window")
Ilmenau_R_pub <- load_published_R_estimates("ilmenau")
SDSC_R_pub <- load_published_R_estimates("sdsc")
globalrt_R_pub <- load_published_R_estimates("globalrt_7d")
#epinow_R_pub <- load_published_R_estimates("epiforecasts")
rtlive_R_pub <- load_published_R_estimates("rtlive")

# merge estimates and plot for comparison
estimates <- RKI_R_pub[,c("date", "R_pub")] %>%
  full_join(ETH_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(Ilmenau_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(SDSC_R_pub[,c("date", "R_pub")], by = "date") %>% 
  full_join(globalrt_R_pub[,c("date", "R_pub")], by = "date") %>%
  #full_join(epinow_R_pub[,c("date", "R_pub")], by = "date") %>%
  full_join(rtlive_R_pub[,c("date", "R_pub")], by = "date")
names(estimates) <- c("date", "RKI", "ETH", "Ilmenau", "SDSC", "globalrt",
                      #"epiforecasts",
                      "rtlive")

plot_multiple_estimates(estimates[estimates$date > "2021-04-01",],
                        methods = c("RKI", "ETH", "Ilmenau", "SDSC", "globalrt",
                                    #"epiforecasts",
                                    "rtlive"))

latest_estimates <- estimates[rowSums(is.na(estimates)) == 0,]

n <- dim(latest_estimates)[2] - 1
comp_methods <- names(latest_estimates)[2:(n+1)]
matr <- matrix(rep(rep(0,n), n), ncol=n)
corr <- matrix(rep(rep(0,n), n), ncol=n)
colnames(matr) <- rownames(matr) <- colnames(corr) <- rownames(corr) <- comp_methods

par(mfrow=c(6,6))
for (method1 in comp_methods) {
  for (method2 in comp_methods){
    diff <- latest_estimates[,method1] - latest_estimates[,method2]
    if (plot_all_differences){
      plot(diff, type="l", ylab = "diff", ylim=c(-0.5, 0.5),
           main = paste(method1, "vs.", method2))
    }
    matr[method1, method2] <- mean(abs(diff))
    corr[method1, method2] <- cor(latest_estimates[,method1], latest_estimates[,method2])
  }
}
par(mfrow=c(1,1))

pheatmap(matr, color = viridis(100), breaks = seq(0,0.4,0.4/100),
         border_color = NA, display_numbers = TRUE,
         fontsize = 12, fontsize_number=20, number_color = "white",
         angle_col = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
         main = paste("Mean absolute differences between published estimates over",
                      max(latest_estimates$date) - min(latest_estimates$date), "days"))

pheatmap(corr, color = viridis(100, direction = -1), breaks = seq(0.5,1,0.5/100),
         border_color = NA, display_numbers = TRUE,
         fontsize = 12, fontsize_number=20, number_color = "white",
         angle_col = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
         main = paste("Correlations between published estimates over",
                      max(latest_estimates$date) - min(latest_estimates$date),
                      "days"))


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
                                delay = params[method, "delay"])[,c("date", "0.5")]
names(R_Ilmenau_d)[2] <- "R_calc"

path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- max(list.files(path, pattern = paste0("R_calc_\\d{4}-\\d{2}-\\d{2}", method, "_Delays.qs")))
print(paste("Most recent epiforecasts estimates with matching parameters:", file))
R_epiforecasts_d <- qread(paste0(path, file))
R_epiforecasts_d <- R_epiforecasts_d[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_d) <- c("date", "R_calc", "0.025", "0.975")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
file <- max(list.files(path, pattern = paste0("estimated_R_\\d{4}-\\d{2}-\\d{2}_", method, "_Delays.csv")))
print(paste("Most recent globalrt estimates with matching parameters:", file))
R_globalrt_d <- read_csv(paste0(path, file))
R_globalrt_d <- R_globalrt_d[R_globalrt_d$`Country/Region` == "Germany", c("Date", "R")]
names(R_globalrt_d) <- c("date", "R_calc")

# merge estimates and plot for comparison
estimates <- R_raw_EpiEstim_d[,c("date", "R_calc")] %>%
  full_join(R_ETH_EpiEstim_d[,c("date", "R_calc")], by = "date") %>% 
  full_join(R_AGES_EpiEstim_d[,c("date", "R_calc")], by = "date") %>% 
  full_join(R_Ilmenau_d[,c("date", "R_calc")], by = "date") %>% 
  full_join(R_epiforecasts_d[,c("date", "R_calc")], by = "date") %>%
  full_join(R_globalrt_d[,c("date", "R_calc")], by = "date")
names(estimates) <- c("date", "rawEpiEstim", "ETH", "AGES", "Ilmenau",
                      "epiforecasts",
                      "globalrt")

plot_multiple_estimates(estimates[estimates$date > "2021-04-12" & estimates$date < "2021-10-01",],
                        methods = c("raw EpiEstim", "ETH EpiEstim", "AGES EpiEstim",
                                    "Ilmenau",
                                    "epiforecasts",
                                    "globalrt"))

latest_estimates <- estimates[rowSums(is.na(estimates)) == 0,]
latest_estimates <- latest_estimates[latest_estimates$date > "2021-04-12",]

n <- dim(latest_estimates)[2] - 1
comp_methods <- names(latest_estimates)[2:(n+1)]
matr <- matrix(rep(rep(0,n), n), ncol=n)
corr <- matrix(rep(rep(0,n), n), ncol=n)
colnames(matr) <- rownames(matr) <- colnames(corr) <- rownames(corr) <- comp_methods

par(mfrow=c(6,6))
for (method1 in comp_methods) {
  for (method2 in comp_methods){
    diff <- latest_estimates[,method1] - latest_estimates[,method2]
    if (plot_all_differences){
      plot(diff, type="l", ylab = "diff", ylim=c(-0.5, 0.5),
           main = paste(method1, "vs.", method2))
    }
    matr[method1, method2] <- mean(abs(diff))
    corr[method1, method2] <- cor(latest_estimates[,method1], latest_estimates[,method2])
  }
}
par(mfrow=c(1,1))

pheatmap(matr, color = viridis(100), breaks = seq(0,0.4,0.4/100),
         border_color = NA, display_numbers = TRUE,
         fontsize = 12, fontsize_number=20, number_color = "white",
         angle_col = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
         main = paste("Mean absolute differences between estimates over",
                      max(latest_estimates$date) - min(latest_estimates$date),
                      "days using", method, "delays"))

pheatmap(corr, color = viridis(100, direction = -1), breaks = seq(0.5,1,0.5/100),
         border_color = NA, display_numbers = TRUE,
         fontsize = 12, fontsize_number=20, number_color = "white",
         angle_col = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
         main = paste("Correlations between estimates over",
                      max(latest_estimates$date) - min(latest_estimates$date),
                      "days using", method, "delays"))



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
                                   delay = params["Ilmenau", "delay"])[,c("date", "0.5")]
names(R_Ilmenau_gt)[2] <- "R_calc"

path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"
file <- max(list.files(path, pattern = paste0("R_calc_\\d{4}-\\d{2}-\\d{2}", method, "_GTD.qs")))
print(paste("Most recent epiforecasts estimates with matching parameters:", file))
R_epiforecasts_gt <- qread(paste0(path, file))
R_epiforecasts_gt <- R_epiforecasts_gt[,c("date", "mean", "lower_95", "upper_95")]
names(R_epiforecasts_gt) <- c("date", "R_calc", "0.025", "0.975")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"
file <- max(list.files(path, pattern = paste0("estimated_R_\\d{4}-\\d{2}-\\d{2}_", method, "_GTD.csv")))
print(paste("Most recent globalrt estimates with matching parameters:", file))
R_globalrt_gt <- read_csv(paste0(path, file))
R_globalrt_gt <- R_globalrt_gt[R_globalrt_gt$`Country/Region` == "Germany", c("Date", "R")]
names(R_globalrt_gt) <- c("date", "R_calc")

# merge estimates and plot for comparison
estimates <- R_raw_EpiEstim_gt[,c("date", "R_calc")] %>%
  full_join(R_ETH_EpiEstim_gt[,c("date", "R_calc")], by = "date") %>% 
  full_join(R_AGES_EpiEstim_gt[,c("date", "R_calc")], by = "date") %>% 
  full_join(R_Ilmenau_gt[,c("date", "R_calc")], by = "date") %>% 
  #full_join(R_epiforecasts_gt[,c("date", "R_calc")], by = "date") %>%
  full_join(R_globalrt_gt[,c("date", "R_calc")], by = "date")
names(estimates) <- c("date", "rawEpiEstim", "ETH", "AGES", "Ilmenau",
                      #"epiforecasts",
                      "globalrt")

plot_multiple_estimates(estimates[estimates$date > "2021-04-12" & estimates$date < "2021-10-01",],
                        methods = c("raw EpiEstim", "ETH EpiEstim", "AGES EpiEstim",
                                    "Ilmenau",
                                    #"epiforecasts",
                                    "globalrt"))

latest_estimates <- estimates[rowSums(is.na(estimates)) == 0,]
latest_estimates <- latest_estimates[latest_estimates$date > "2021-04-12",]

n <- dim(latest_estimates)[2] - 1
comp_methods <- names(latest_estimates)[2:(n+1)]
matr <- matrix(rep(rep(0,n), n), ncol=n)
corr <- matrix(rep(rep(0,n), n), ncol=n)
colnames(matr) <- rownames(matr) <- colnames(corr) <- rownames(corr) <- comp_methods

par(mfrow=c(6,6))
for (method1 in comp_methods) {
  for (method2 in comp_methods){
    diff <- latest_estimates[,method1] - latest_estimates[,method2]
    if (plot_all_differences){
      plot(diff, type="l", ylab = "diff", ylim=c(-0.5, 0.5),
           main = paste(method1, "vs.", method2))
    }
    matr[method1, method2] <- mean(abs(diff))
    corr[method1, method2] <- cor(latest_estimates[,method1], latest_estimates[,method2])
  }
}
par(mfrow=c(1,1))

pheatmap(matr, color = viridis(100), breaks = seq(0,0.4,0.4/100),
         border_color = NA, display_numbers = TRUE,
         fontsize = 12, fontsize_number=20, number_color = "white",
         angle_col = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
         main = paste("Mean absolute differences between estimates over",
                      max(latest_estimates$date) - min(latest_estimates$date),
                      "days using", method, "generation time"))

pheatmap(corr, color = viridis(100, direction = -1), breaks = seq(0.5,1,0.5/100),
         border_color = NA, display_numbers = TRUE,
         fontsize = 12, fontsize_number=20, number_color = "white",
         angle_col = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
         main = paste("Correlations between estimates over",
                      max(latest_estimates$date) - min(latest_estimates$date),
                      "days using", method, "generation time"))



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






