setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

# parameter combinations used in papers
window_size <- c(3, 7, 1, 4, 13, 1, 1, 1)
gt_dist <- c("gamma", "constant", "ad hoc", "gamma", "gamma", "?", "log-normal", "exponential")
mean_gt <- c(4.8, 4.0, 5.6, 4.8, 3.4, 3.6, 4.7, 7)
sd_gt <- c(2.3, 0.0, 4.2, 2.3, 1.8, 3.1, 2.9, 7)

params <- data.frame(window=window_size, gtd=gt_dist, gt_mean=mean_gt, gt_sd=sd_gt)
rownames(params) <- c("ETH", "RKI", "Ilmenau", "SDSC", "AGES", "epiforecasts", "rtlive", "globalrt")

method <- "ETH"

# use same data for all methods
incid <- load_incidence_data(method = "ETHZ_sliding_window")
simple_incid <- incid[, c("date", "value")]
simple_incid <- aggregate(simple_incid$value, by=list(simple_incid$date), FUN=mean)
names(simple_incid) <- c("date", "I")


# estimations
RKI_R <- estimate_RKI_R(simple_incid,
                        window=params[method, "window"],
                        gt_mean=params[method, "gt_mean"],
                        gt_sd=params[method, "gt_sd"])

ETH_R <- estimate_ETH_R(incid)#,
                        #window=params[method, "window"],
                        #gt_mean=params[method, "gt_mean"],
                        #gt_sd=params[method, "gt_sd"])

Ilmenau_R <- estimate_Ilmenau_R(simple_incid, gt_type = "gamma",
                                window=params[method, "window"],
                                gt_mean=params[method, "gt_mean"],
                                gt_sd=params[method, "gt_sd"])

AGES_R <- estimate_AGES_R(simple_incid,
                          window=params[method, "window"],
                          gt_mean=params[method, "gt_mean"],
                          gt_sd=params[method, "gt_sd"])

SDSC_R <- estimate_SDSC_R(simple_incid, estimateOffsetting = 7,
                          window=params[method, "window"],
                          gt_mean=params[method, "gt_mean"],
                          gt_sd=params[method, "gt_sd"])

# merge estimates and plot for comparison
estimates <- RKI_R %>%
  full_join(ETH_R, by = "date") %>% 
  full_join(Ilmenau_R, by = "date") %>%
  full_join(AGES_R, by = "date") %>%
  full_join(SDSC_R, by = "date")
names(estimates) <- c("date", "RKI", "ETH", "Ilmenau", "AGES", "SDSC")

plot_multiple_estimates(estimates)

