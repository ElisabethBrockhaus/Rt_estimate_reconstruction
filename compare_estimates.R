# git repositories "reproductive_numbers" and "Rt_estimate_reconstruction"
# have to be located in the same directory
setwd("..")

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

##############
# RKI (AnDerHeiden2020)
##############

# load data
RKI_data <- load_RKI_data()

# estimation
RKI_est <- estimate_RKI_R(RKI_data) # only allows gt_sd = 0
RKI_EpiEstim_est <- estimate_RKI_R_EpiEstim(RKI_data, window = 7,
                                            gt_mean = 4, gt_sd = 0) # gt_sd > 0 possible

# plots for comparison
plot_published_vs_calculated(RKI_data, RKI_est, method_name="RKI")
plot_published_vs_calculated(RKI_data, RKI_EpiEstim_est, method_name="RKI with EpiEstim")



##############
# ETH (Huisman2021)
##############

# load data
ETH_data <- load_ETH_data()

# estimation
#ETH_est <- estimate_ETH_R(ETH_data)
ETH_est <- countryEstimates[countryEstimates$data_type=="Confirmed cases" &
                              countryEstimates$estimate_type=="Cori_slidingWindow",]
ETH_est <- data.frame(date=ETH_est$date,
                      R_calc=ETH_est$median_R_mean)

# plots for comparison
plot_published_vs_calculated(published=ETH_data[,c("date", "R_pub")],
                             calculated=ETH_est,
                             method_name="ETH")



##############
# Ilmenau (Hotz2020)
##############

# load data
Ilmenau_data <- load_Ilmenau_data()

# estimation
Ilmenau_est <- estimate_Ilmenau_R(Ilmenau_data, gt_type = "org")

# plots for comparison
plot_published_vs_calculated(Ilmenau_data, Ilmenau_est, method_name="Ilmenau")



##############
# AGES (Richter2020)
##############

# load data (Austria only)
AGES_data <- load_AGES_data()

# estimation with mean/sd used since 18th June 2021
AGES_est <- estimate_AGES_R(AGES_data, mean_si = 3.37, std_si = 1.83)

# plots for comparison
plot_published_vs_calculated(AGES_data, AGES_est, method_name="AGES")



##############
# epiforecasts (Abbot2020)
##############

# load data
#epiforecasts_data <- load_published_R_estimates("epiforecasts",
#                                                incid=data.frame(date=seq(as.Date("2021-04-06"), as.Date("2021-07-27"), by="day")))
estimates_published <- read_csv("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/master/national/cases/summary/rt.csv")
estimates_published <- estimates_published[estimates_published$country=="Germany", c("date", "median")]
names(estimates_published) <- c("date", "R_pub")

# estimation
#EpiNow2_est <- estimate_EpiNow2_R(epiforecasts_data)
EpiNow2_est <- qread("Rt_estimate_reconstruction/epiforecast_estimates/EpiNow2_est_horizon14.qs")
names(EpiNow2_est) <- c("date", "R_calc")

# plots for comparison
plot_published_vs_calculated(estimates_published, EpiNow2_est, method_name="EpiNow2")




###################

# compare estimates from different sources
estimates <- ETH_data[, c("date", "R_pub")] %>%
  full_join(Ilmenau_data[, c("date", "R_pub")], by = "date") %>% 
  full_join(RKI_data[, c("date", "R_pub")], by = "date") %>%
  full_join(estimates_published[, c("date", "R_pub")], by = "date")
names(estimates) <- c("date", "ETH", "Ilmenau", "RKI", "epiforecasts")

plot_multiple_estimates(estimates)

# compare estimates from different methods with same window size
# 7 day
RKI_est7 <- estimate_RKI_R(RKI_data, window = 7)
Ilmenau_est7 <- estimate_Ilmenau_R(Ilmenau_data, window = 7)

estimates <- RKI_est7[, c("date", "R_calc")] %>%
  full_join(Ilmenau_est7[, c("date", "R_calc")], by = "date") %>%
  full_join(EpiNow2_est[, c("date", "R_calc")], by = "date")
names(estimates) <- c("date", "RKI", "Ilmenau", "epiforecasts")

plot_multiple_estimates(estimates)

# 3 day
RKI_est3 <- estimate_RKI_R(RKI_data, window = 3)
Ilmenau_est3 <- estimate_Ilmenau_R(Ilmenau_data, window = 3)

estimates <- RKI_est3[, c("date", "R_calc")] %>%
  full_join(Ilmenau_est3[, c("date", "R_calc")], by = "date") %>%
  full_join(ETH_data[, c("date", "R_pub")], by = "date")
names(estimates) <- c("date", "RKI", "Ilmenau", "ETH")

plot_multiple_estimates(estimates)

Ilmenau_est3_shifted <- estimate_Ilmenau_R(Ilmenau_data, gt_type = "gamma")
Ilmenau_est3_shifted$date <- Ilmenau_est3_shifted$date+5
estimates <- RKI_est3[, c("date", "R_calc")] %>%
  full_join(Ilmenau_est3_shifted[, c("date", "R_calc")], by = "date")
names(estimates) <- c("date", "RKI", "Ilmenau")

plot_multiple_estimates(estimates)
