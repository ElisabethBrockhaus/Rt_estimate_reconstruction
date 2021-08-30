setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

##############
# RKI (AnDerHeiden2020)
##############

# load data
RKI_incid <- load_incidence_data(method = "RKI")
RKI_R_pub <- load_published_R_estimates(source = "RKI_7day")

# estimation
RKI_R_calc <- estimate_RKI_R(RKI_incid)

# plots for comparison
plot_published_vs_calculated(RKI_R_pub, RKI_R_calc, method_name="RKI")

# with slightly higher gt_sd
RKI_R_calc_EpiEstim <- estimate_RKI_R(RKI_incid, gt_sd = 1)
plot_published_vs_calculated(RKI_R_pub, RKI_R_calc_EpiEstim, method_name="RKI (gt_sd > 0)")



##############
# ETH (Huisman2021)
##############

# load data
ETH_countryData <-  load_incidence_data(method = "ETHZ_sliding_window")
ETH_R_pub <- load_published_R_estimates("ETHZ_sliding_window")

# estimation
ETH_R_calc <- estimate_ETH_R(ETH_countryData)

# plots for comparison
plot_published_vs_calculated(published=ETH_R_pub, calculated=ETH_R_calc, method_name="ETH")

###
# for Austria
# load data
ETH_countryData_AUT <-  load_incidence_data(method = "ETHZ_sliding_window", location = "AT")
ETH_R_pub_AUT <- load_published_R_estimates("ETHZ_sliding_window", location = "AU")

# estimation
ETH_R_calc_AUT <- estimate_ETH_R(ETH_countryData_AUT)

# plots for comparison
plot_published_vs_calculated(published=ETH_R_pub_AUT, calculated=ETH_R_calc_AUT, method_name="ETH (Austria)")


##############
# Ilmenau (Hotz2020)
##############

# load data
Ilmenau_incid <- load_incidence_data("ilmenau")
Ilmenau_R_pub <- load_published_R_estimates("ilmenau")

# estimation
Ilmenau_R_calc <- estimate_Ilmenau_R(Ilmenau_incid, gt_type = "org")

# plots for comparison
plot_published_vs_calculated(Ilmenau_R_pub, Ilmenau_R_calc, method_name="Ilmenau")



##############
# AGES (Richter2020)
##############

# load data (Austria only)
AGES_incid <- load_incidence_data("AGES", location = "AT")
AGES_R_pub <- load_published_R_estimates("AGES", location = "AT")

# estimation with mean/sd used since 18th June 2021
AGES_R_calc <- estimate_AGES_R(AGES_incid, mean_si = 3.37, std_si = 1.83)

# plots for comparison
plot_published_vs_calculated(AGES_R_pub, AGES_R_calc, method_name="AGES")



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
#EpiNow2_R_calc <- estimate_EpiNow2_R(epiforecasts_data)
EpiNow2_R_calc <- qread("Rt_estimate_reconstruction/epiforecasts/estimates/EpiNow2_est_horizon14.qs")
names(EpiNow2_R_calc) <- c("date", "R_calc")

# plots for comparison
plot_published_vs_calculated(estimates_published, EpiNow2_R_calc, method_name="EpiNow2")



##############
# SDSC
##############

# load data
SDSC_R_pub <- load_published_R_estimates("sdsc")
# set data_status to one day after publication day of R_pub
SDSC_incid <- load_incidence_data("sdsc", data_status = "2021-08-29")

# estimation
SDSC_R_calc <- estimate_SDSC_R(SDSC_incid, estimateOffsetting = 7)

# plots for comparison
plot_published_vs_calculated(SDSC_R_pub, SDSC_R_calc, method_name="SDSC")



###################

# compare estimates from different sources
estimates <- ETH_R_pub %>%
  full_join(Ilmenau_R_pub, by = "date") %>% 
  full_join(RKI_R_pub, by = "date") %>%
  full_join(estimates_published, by = "date") %>%
  full_join(SDSC_R_pub, by = "date")
names(estimates) <- c("date", "ETH", "Ilmenau", "RKI", "epiforecasts", "SDSC")

plot_multiple_estimates(estimates)

# compare estimates from different methods with same window size
# 7 day
RKI_est7 <- estimate_RKI_R(RKI_incid, window = 7)
Ilmenau_est7 <- estimate_Ilmenau_R(Ilmenau_incid, window = 7, gt_type="org")
SDSC_est7 <- estimate_SDSC_R(SDSC_incid, estimateOffsetting = 7, window=7)

estimates <- RKI_est7 %>%
  full_join(Ilmenau_est7, by = "date") %>%
  full_join(EpiNow2_R_calc, by = "date") %>%
  full_join(SDSC_est7, by = "date")
names(estimates) <- c("date", "RKI", "Ilmenau", "epiforecasts", "SDSC")

plot_multiple_estimates(estimates[estimates$date > "2020-03-08",])

# 3 day
RKI_est3 <- estimate_RKI_R(RKI_incid, window = 3)
Ilmenau_est3 <- estimate_Ilmenau_R(Ilmenau_incid, window = 3, gt_type = "org")
SDSC_est3 <- estimate_SDSC_R(SDSC_incid, estimateOffsetting = 7, window=3)

estimates <- RKI_est3 %>%
  full_join(Ilmenau_est3, by = "date") %>%
  full_join(ETH_R_pub, by = "date") %>%
  full_join(SDSC_est3, by = "date")
names(estimates) <- c("date", "RKI", "Ilmenau", "ETH", "SDSC")

plot_multiple_estimates(estimates[estimates$date > "2020-03-08",])

# 3 day Austria
AGES_est3 <- estimate_AGES_R(AGES_incid, window=3)

estimates <- AGES_est3 %>%
  full_join(ETH_R_calc_AUT, by = "date")
names(estimates) <- c("date", "AGES", "ETH (Austria)")

plot_multiple_estimates(estimates)

# use Ilmenau method with data and parameters from RKI
Ilmenau_est3_shifted <- estimate_Ilmenau_R(RKI_incid, window = 3, gt_type = "gamma", gt_mean=4, gt_sd=0.0001)
Ilmenau_est3_shifted$date <- Ilmenau_est3_shifted$date+6
estimates <- RKI_est3 %>%
  full_join(Ilmenau_est3_shifted, by = "date")
names(estimates) <- c("date", "RKI", "Ilmenau")

plot_multiple_estimates(estimates)

# use RKI method with mean of deconvoluted ETH data
ETH_incid <- ETH_countryData[, c("date", "value")]
ETH_incid <- aggregate(ETH_incid$value, by=list(ETH_incid$date), FUN=mean)
names(ETH_incid) <- c("date", "I")

RKI_est3_ETH <- estimate_RKI_R(ETH_incid, window = 3)
estimates <- RKI_est3_ETH %>% full_join(ETH_R_calc, by = "date")
names(estimates) <- c("date", "RKI with ETH data", "ETH")

plot_multiple_estimates(estimates)

