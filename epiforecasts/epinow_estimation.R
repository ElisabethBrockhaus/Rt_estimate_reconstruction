library(data.table)
library(EpiNow2)
library(covidregionaldata)
library(dplyr)
library(lubridate)
library(qs)
library(readr)

wd <- "D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code/Rt_estimate_reconstruction/epiforecasts"
setwd(wd)

country <- "Germany"

# load published estimates
estimates_published <- read_csv("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/master/national/cases/summary/rt.csv")
estimates_published <- estimates_published[estimates_published$country==country,]

# load reported cases with function from covidregionaldata (complete history)
reported_cases <- get_national_data(countries = country)
reported_cases <- reported_cases[, c("date", "cases_new")]
names(reported_cases) <- c("date", "confirm")

# estimate only on 16 weeks of data due to computing times
reported_cases_latest <- reported_cases[reported_cases$date >= min(estimates_published$date) &
                                          reported_cases$date <= max(estimates_published[estimates_published$type != "forecast",]$date),]

# load distributions that were published at Github (https://github.com/seabbs/covid-climate-rt/tree/master/delays/data)
generation_time <- readRDS("distributions/generation_time.rds")
incubation_period <- readRDS("distributions/incubation_period.rds")
reporting_delay <- readRDS("distributions/onset_to_admission_delay.rds")

# do estimation
estimates <- epinow(reported_cases = reported_cases_latest,
                    generation_time = generation_time,
                    delays = delay_opts(incubation_period, reporting_delay),
                    rt = rt_opts(prior = list(mean = 1, sd = 0.2)),
                    stan = stan_opts(cores = 4),
                    horizon = 14,
                    CrIs = c(0.5, 0.95),
                    verbose = TRUE)

# plot estimates
plot(estimates)

# extract R estimates
epiforecasts_R_calc <- estimates$estimates$summarised[variable=="R", c("date", "type", "mean")]
qsave(epiforecasts_R_calc, paste0("R_calc_", Sys.Date(), ".qs"))

# compare published and calculated in plot
plot(epiforecasts_R_calc$date, epiforecasts_R_calc$mean, type="l")
lines(estimates_published$date, estimates_published$mean, col="red")

plot(epiforecasts_R_calc$date,
     epiforecasts_R_calc$mean - estimates_published[estimates_published$date %in% epiforecasts_R_calc$date,]$mean,
     type="l")
abline(h=0, col="grey", lty=2)



############################
# use different parameters #
############################

# load incidence data from RKI
incid <- read_csv("RKI_incid.csv")
incid_latest <- incid[incid$date >= as.Date("2021-06-05") &
                        incid$date <= max(estimates_published[estimates_published$type != "forecast",]$date),]
names(incid_latest) <- c("date", "confirm")

# set ETH parameters
generation_time_eth <- list(
  mean = 4.8, mean_sd = 0.1,
  sd = 2.3, sd_sd = 0.1,
  max = 30
)

incubation_period_eth <- list(
  mean = convert_to_logmean(5.3, 3.2), mean_sd = 0.1,
  sd = convert_to_logsd(5.3, 3.2), sd_sd = 0.1,
  max = 30
)

reporting_delay_eth <- list(
  mean = convert_to_logmean(5.5, 3.8), mean_sd = 0.1,
  sd = convert_to_logsd(5.5, 3.8), sd_sd = 0.1,
  max = 30
)

# do estimation
estimates_eth_params <- epinow(reported_cases = incid_latest,
                               generation_time = generation_time_eth,
                               delays = delay_opts(incubation_period_eth, reporting_delay_eth),
                               rt = rt_opts(prior = list(mean = 1, sd = 0.2)),
                               stan = stan_opts(cores = 4),
                               horizon = 14,
                               CrIs = c(0.5, 0.95),
                               verbose = TRUE)

plot(estimates_eth_params)

# compare published and calculated in plot
epiforecasts_R_calc_eth <- estimates_eth_params$estimates$summarised[variable=="R", c("date", "type", "mean")]
qsave(epiforecasts_R_calc_eth, paste0("R_calc_", Sys.Date(), "_ETHParams.qs"))

plot(epiforecasts_R_calc_eth$date, epiforecasts_R_calc_eth$mean, type="l")



# set rtlive parameters
generation_time_rtlive <- list(
  mean = 4.7, mean_sd = 0.1,
  sd = 2.9, sd_sd = 0.1,
  max = 30
)

incubation_period_rtlive <- list(
  mean = convert_to_logmean(5, 0), mean_sd = 0.1,
  sd = convert_to_logsd(5, 0), sd_sd = 0.1,
  max = 30
)

reporting_delay_rtlive <- list(
  mean = convert_to_logmean(7.1, 5.9), mean_sd = 0.1,
  sd = convert_to_logsd(7.1, 5.9), sd_sd = 0.1,
  max = 30
)

# do estimation
estimates_rtlive_params <- epinow(reported_cases = incid_latest,
                                  generation_time = generation_time_rtlive,
                                  delays = delay_opts(incubation_period_rtlive, reporting_delay_rtlive),
                                  rt = rt_opts(prior = list(mean = 1, sd = 0.2)),
                                  stan = stan_opts(cores = 4),
                                  horizon = 14,
                                  CrIs = c(0.5, 0.95),
                                  verbose = TRUE)

plot(estimates_rtlive_params)

# compare published and calculated in plot
epiforecasts_R_calc_rtlive <- estimates_rtlive_params$estimates$summarised[variable=="R", c("date", "type", "mean")]
qsave(epiforecasts_R_calc_rtlive, paste0("R_calc_", Sys.Date(), "_rtliveParams.qs"))

plot(epiforecasts_R_calc_rtlive$date, epiforecasts_R_calc_rtlive$mean, type="l")