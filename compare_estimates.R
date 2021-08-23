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
RKI14_est <- estimate_RKI_R(RKI_data, window = 14)
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
ETH_est <- estimate_ETH_R(ETH_data)

# plots for comparison
plot_published_vs_calculated(published=data.frame(dates=ETH_data$date,
                                                  R_pub=ETH_data$R),
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
epiforecasts_data <- load_published_R_estimates("epiforecasts",
                                                incid=data.frame(date=seq(as.Date("2021-04-06"), as.Date("2021-07-27"), by="day")))

# check if data in git reproductive_numbers might be wrong (MV instead of national)
subnational_estimates <- read_csv("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/master/subnational/germany/cases/summary/rt.csv")
MV_estimates <- subnational_estimates[subnational_estimates$region == "Mecklenburg-Vorpommern", c("date", "median", "mean")]

plot(MV_estimates$date, MV_estimates$mean, type="l")
lines(epiforecasts_data$date, epiforecasts_data$R_pub, col="red")

# estimation
#EpiNow2_est <- estimate_EpiNow2_R(epiforecasts_data)
EpiNow2_est <- qread("EpiNow2_est_blue_line.qs")
names(EpiNow2_est) <- c("date", "R_calc")
EpiNow2_est <- full_join(data.frame(date=epiforecasts_data$date), EpiNow2_est, by="date")
EpiNow2_est <- EpiNow2_est[EpiNow2_est$date >= min(epiforecasts_data$date)
                           & EpiNow2_est$date <= max(epiforecasts_data$date),]

# plots for comparison
plot_published_vs_calculated(epiforecasts_data, EpiNow2_est, method_name="EpiNow2")

epiforecast_estimates <- read_csv("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/master/national/cases/summary/rt.csv")
epiforecast_estimates <- epiforecast_estimates[epiforecast_estimates$country=="Germany", c("date", "median", "mean")]
names(epiforecast_estimates) <- c("date", "R_pub", "R_mean")

estimates <- epiforecast_estimates[, c("date", "R_pub")] %>% 
  full_join(EpiNow2_est[, c("date", "R_calc")], by = "date")
names(estimates) <- c("date", "pub", "calc")

plot_multiple_estimates(estimates)

###################

# compare estimates from different sources
estimates <- epiforecasts_data[, c("date", "R_pub")] %>% 
  full_join(ETH_data[, c("date", "R_pub")], by = "date") %>%
  full_join(Ilmenau_data[, c("date", "R_pub")], by = "date") %>% 
  full_join(RKI_data[, c("date", "R_pub")], by = "date")
names(estimates) <- c("date", "epiforecasts", "ETH", "Ilmenau", "RKI")

plot_multiple_estimates(estimates)

