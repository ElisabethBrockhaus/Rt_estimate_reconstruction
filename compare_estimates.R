setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

#########################
# RKI (AnDerHeiden2020) #
#########################

# load data
RKI_incid <- load_incidence_data(method = "RKI")
RKI_R_pub <- load_published_R_estimates(source = "RKI_7day")

# estimation
RKI_R_calc <- estimate_RKI_R(RKI_incid)

RKI_R_calc_EpiEstim <- estimate_RKI_R(RKI_incid, method="EpiEstim")
names(RKI_R_calc_EpiEstim) <- c("date", "0.5", "0.025", "0.975")

RKI_R_calc_Ilmenau <- estimate_Ilmenau_R(RKI_incid, window = 7, shift = 6,
                                         gt_type = "gamma", gt_mean = 4, gt_sd = 0.0001)

# plots for comparison
plot_published_vs_calculated(RKI_R_pub, RKI_R_calc, method_name="RKI")
plot_published_vs_calculated_95CI(RKI_R_pub[RKI_R_pub$date > "2020-03-12",],
                                  RKI_R_calc_EpiEstim[RKI_R_calc_EpiEstim$date > "2020-03-12",],
                                  method_name="RKI")
plot_published_vs_calculated_95CI(RKI_R_pub[RKI_R_pub$date > "2020-03-12",],
                                  RKI_R_calc_Ilmenau[RKI_R_calc_Ilmenau$date > "2020-03-12",],
                                  method_name="RKI")

# with slightly higher gt_sd
RKI_R_calc_gamma <- estimate_RKI_R(RKI_incid, gt_type="gamma", gt_sd = 1)
plot_published_vs_calculated(RKI_R_pub, RKI_R_calc_gamma, method_name="RKI (gt_sd = 1)")


#####################
# ETH (Huisman2021) #
#####################

# load data
ETH_countryData <-  load_incidence_data(method = "ETHZ_sliding_window")#, force_deconvolution=T)
ETH_R_pub <- load_published_R_estimates("ETHZ_sliding_window")

# estimation
ETH_R_calc <- estimate_ETH_R(ETH_countryData)

# plots for comparison
plot_published_vs_calculated(published=ETH_R_pub, calculated=ETH_R_calc, method_name="ETH")


# for Austria
# load data
ETH_countryData_AUT <-  load_incidence_data(method = "ETHZ_sliding_window", location = "AT")
ETH_R_pub_AUT <- load_published_R_estimates("ETHZ_sliding_window", location = "AT")

# estimation
ETH_R_calc_AUT <- estimate_ETH_R(ETH_countryData_AUT)

# plots for comparison
plot_published_vs_calculated(published=ETH_R_pub_AUT, calculated=ETH_R_calc_AUT, method_name="ETH (Austria)")


######################
# Ilmenau (Hotz2020) #
######################

# load data
Ilmenau_incid <- load_incidence_data("ilmenau")
Ilmenau_R_pub <- load_published_R_estimates("ilmenau")

# estimation
Ilmenau_R_calc <- estimate_Ilmenau_R(Ilmenau_incid, gt_type = "ad hoc")

# plots for comparison
plot_published_vs_calculated_95CI(Ilmenau_R_pub, Ilmenau_R_calc, method_name="Ilmenau")

# compare confidence intervals to the ones from EpiEstim
Ilmenau_R_calc_EpiEstim <- estimate_RKI_R(Ilmenau_incid, method="EpiEstim", gt_type="ad hoc", window = 1, shift=-6)
names(Ilmenau_R_calc_EpiEstim) <- names(Ilmenau_R_pub)
Ilmenau_R_calc_EpiEstim <- Ilmenau_R_calc_EpiEstim[Ilmenau_R_calc_EpiEstim$date >= "2020-02-23",]
plot_published_vs_calculated_95CI(Ilmenau_R_calc_EpiEstim, Ilmenau_R_calc, method_name = "EpiEstim vs. Ilmenau confidence intervals")


######################
# AGES (Richter2020) #
######################

# load data (Austria only)
AGES_incid <- load_incidence_data("AGES", location = "AT")
AGES_R_pub <- load_published_R_estimates("AGES", location = "AT")

# estimation with mean/sd used since 18th June 2021
AGES_R_calc <- estimate_AGES_R(AGES_incid,
                               gt_mean = 3.37, gt_sd = 1.83,
                               shape = 3.38, scale = 1.00)

a <- 3.38
b <- 1.00
window <- 13 
gt_dist <- c(0, dgamma(1:1000, shape = a, scale = b))

for (t in window:nrow(AGES_incid)){
  
  # numerator
  new_infect <- a + sum(AGES_incid$I[t - 0:(window-1)])
  
  # nominator
  total_infect <- 1/b
  for (i in t - 0:(window-1)) {
    total_infect <- total_infect + (AGES_incid$I[i:1] %*% gt_dist[1:i])
  }
  
  # R_t
  AGES_R_calc[t, "R_calc"] <- new_infect / total_infect
}

plot_published_vs_calculated(AGES_R_pub, AGES_R_calc, method_name="AGES (Austria)")

incid <- AGES_incid
start_i <- min(which(!is.na(incid$I)))
end_i <- max(which(!is.na(incid$I)))
incid_wo_na <- incid[start_i:end_i,]

# start and end for estimations at each time point
start <- 2:(nrow(incid_wo_na) + 1 - window)
end <- start - 1 + window

# estimation with EpiEstim function
r_EpiEstim <- estimate_R(incid_wo_na$I,
                         method = "non_parametric_si",
                         config = make_config(list(
                           t_start=start, t_end=end,
                           si_distr=gt_dist,
                           mean_prior = 1,
                           std_prior = 5)))

len_est <- length(r_EpiEstim$R$`Mean(R)`)

# assign estimates to time points properly
estimate <- data.frame(date=incid$date, R_calc=c(rep(NA, (start_i+window-1)),
                                                 r_EpiEstim$R$`Mean(R)`))

# plots for comparison
plot_published_vs_calculated(AGES_R_pub, estimate, method_name="AGES (Austria)")



############################
# epiforecasts (Abbot2020) #
############################

# save new published estimates
#estimates_published <- read_csv("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/master/national/cases/summary/rt.csv")
#estimates_published <- estimates_published[estimates_published$country=="Germany", c("date", "type","median")]
#qsave(estimates_published, paste0("Rt_estimate_reconstruction/epiforecasts/estimates/R_pub_",
#                                  Sys.Date(), ".qs"))

# load estimates published and calculated from same date
data_status <- "2021-09-06"
epiforecasts_R_pub <- qread(paste0("Rt_estimate_reconstruction/epiforecasts/estimates/R_pub_",
                                   data_status, ".qs"))
epiforecasts_R_calc  <- qread(paste0("Rt_estimate_reconstruction/epiforecasts/estimates/R_calc_",
                                     data_status, ".qs"))

names(epiforecasts_R_pub) <- c("date", "type", "R_pub")
names(epiforecasts_R_calc) <- c("date", "type", "R_calc")

# plots for comparison
plot_published_vs_calculated(epiforecasts_R_pub, epiforecasts_R_calc, method_name="epiforecasts")



########
# SDSC #
########

# load data
SDSC_R_pub <- load_published_R_estimates("sdsc")
# set data_status to one day after publication day of R_pub
SDSC_incid <- load_incidence_data("sdsc", data_status = "2021-08-29")

# estimation
SDSC_R_calc <- estimate_SDSC_R(SDSC_incid, estimateOffsetting = 7)

# plots for comparison
plot_published_vs_calculated_95CI(SDSC_R_pub, SDSC_R_calc, method_name="SDSC")



############################
# globalrt (ArroyoMarioli) #
############################

# load published estimates
#globalrt_R_pub <- load_published_R_estimates("globalrt_7d")
globalrt_R_pub <- read_csv("https://raw.githubusercontent.com/crondonm/TrackingR/main/Estimates-Database/database.csv")
globalrt_R_pub <- globalrt_R_pub[globalrt_R_pub$`Country/Region` == "Germany" &
                                   globalrt_R_pub$days_infectious == 7, c("Date", "R", "ci_95_u", "ci_95_l")]
names(globalrt_R_pub) <- c("date", "0.5", "0.025", "0.975")
#names(globalrt_R_pub) <- c("date", "R_pub", "0.025", "0.975")

# load and filter calculated estimates
#globalrt_R_calc  <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/estimates/estimated_R.csv")
#globalrt_R_calc <- globalrt_R_calc[globalrt_R_calc$`Country/Region` == "Germany" &
#                                     globalrt_R_calc$days_infectious == 7, c("Date", "R", "ci_95_u", "ci_95_l")]
globalrt_R_calc  <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/estimates/bayesian_smoother.csv")
names(globalrt_R_calc) <- c("date", "0.5", "0.025", "0.975")
#globalrt_R_calc  <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/estimates/frequentist_estimates.csv")
#names(globalrt_R_calc) <- c("date", "R_filtered", "R_calc")

# plots for comparison
plot_published_vs_calculated_95CI(globalrt_R_pub, globalrt_R_calc, method_name="globalrt")





############################################
# compare estimates from different sources #
############################################

estimates <- ETH_R_pub %>%
  full_join(Ilmenau_R_pub, by = "date") %>% 
  full_join(RKI_R_pub, by = "date") %>%
  full_join(epiforecasts_R_pub, by = "date") %>%
  full_join(SDSC_R_pub, by = "date") %>%
  full_join(globalrt_R_pub, by = "date")

plot_multiple_estimates(estimates, methods = c("ETH", "Ilmenau", "RKI", "epiforecasts", "SDSC", "globalrt"))

# compare estimates from different methods with same window size
# 7 day
RKI_est7 <- estimate_RKI_R(RKI_incid, window = 7)
Ilmenau_est7 <- estimate_Ilmenau_R(Ilmenau_incid, window = 7, gt_type="ad hoc")
ETH_est7 <- estimate_ETH_R(ETH_countryData, window = 7)
SDSC_est7 <- estimate_SDSC_R(SDSC_incid, estimateOffsetting = 7, window=7)

estimates <- RKI_est7 %>%
  full_join(Ilmenau_est7, by = "date") %>%
  full_join(ETH_est7, by = "date") %>%
  full_join(SDSC_est7, by = "date")

plot_multiple_estimates(estimates[estimates$date > "2020-03-08",],
                        methods = c("date", "RKI", "Ilmenau", "ETH", "SDSC"))

# 3 day
RKI_est3 <- estimate_RKI_R(RKI_incid, window = 3)
Ilmenau_est3 <- estimate_Ilmenau_R(Ilmenau_incid, window = 3, gt_type = "ad hoc")
SDSC_est3 <- estimate_SDSC_R(SDSC_incid, estimateOffsetting = 7, window=3)

estimates <- RKI_est3 %>%
  full_join(Ilmenau_est3, by = "date") %>%
  full_join(ETH_R_pub[, c("date", "R_pub")], by = "date") %>%
  full_join(SDSC_est3, by = "date")

plot_multiple_estimates(estimates[estimates$date > "2020-03-08",],
                        methods = c("RKI", "Ilmenau", "ETH", "SDSC"))

# 3 day Austria
AGES_est3 <- estimate_AGES_R(AGES_incid, window=3)

estimates <- AGES_est3 %>%
  full_join(ETH_R_calc_AUT, by = "date")

plot_multiple_estimates(estimates, methods = c("AGES", "ETH (Austria)"))

# use Ilmenau method with data and parameters from RKI
Ilmenau_est3_shifted <- estimate_Ilmenau_R(RKI_incid, window = 3, gt_type = "gamma", gt_mean=4, gt_sd=0.0001)
Ilmenau_est3_shifted$date <- Ilmenau_est3_shifted$date+6
estimates <- RKI_est3 %>%
  full_join(Ilmenau_est3_shifted, by = "date")

plot_multiple_estimates(estimates, methods = c("RKI", "Ilmenau"))

# use RKI method with (mean of) deconvoluted ETH data and ETH parameters

ETH_incid <- ETH_countryData[, c("date", "value")]
ETH_incid <- aggregate(ETH_incid$value, by=list(ETH_incid$date), FUN=mean)
names(ETH_incid) <- c("date", "I")

RKI_est3_ETH <- estimate_RKI_R(ETH_incid, window = 3, gt_type = "gamma",
                               gt_mean = 4.8, gt_sd = 2.3)
estimates <- RKI_est3_ETH %>% full_join(ETH_R_calc, by = "date")

plot_multiple_estimates(estimates, methods = c("RKI (ETH data and parameters)", "ETH"))

# use ETH method with JHU data
ETH_countryDataJHU <-  load_incidence_data(method = "ETHZ_sliding_window",
                                           #force_deconvolution=T,
                                           source="JHU")
ETH_R_calc_JHU <- estimate_ETH_R(ETH_countryDataJHU)
estimates <- ETH_R_calc %>% full_join(ETH_R_calc_JHU, by = "date")

plot_multiple_estimates(as.data.frame(estimates), methods = c("ETH", "ETH with JHU data"))



# compare EpiEstim CI to ETH CI
EpiEstim_est3_ETH <- estimate_RKI_R(ETH_incid, window = 3, method = "EpiEstim",
                                    gt_type = "gamma", gt_mean = 4.8, gt_sd = 2.3)
estimates <- EpiEstim_est3_ETH %>% full_join(ETH_R_pub, by = "date")

plot_multiple_estimates(estimates[estimates$date > "2020-01-27",],
                        methods = c("EpiEstim.ETHpars", "ETH"),
                        include_CI = TRUE)

# compare RKI estimate to estimate when using unif[3,5] for generation time
RKI_R_calc_unif <- estimate_RKI_R(RKI_incid, gt_type = "uniform")
estimates <- RKI_R_calc %>% full_join(RKI_R_calc_unif, by = "date")

plot_multiple_estimates(estimates[estimates$date > "2020-01-27",],
                        methods = c("RKI", "RKI with unif[3,5]"),
                        include_CI = FALSE)
