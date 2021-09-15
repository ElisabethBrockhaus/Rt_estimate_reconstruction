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

### globalrt

# load reported cases with function from covidregionaldata (complete history)
# use source = "jhu" to be in line with globalrt
reported_cases_jhu <- get_national_data(countries = country, source = "jhu")
reported_cases_jhu <- reported_cases_jhu[, c("date", "cases_new")]
names(reported_cases_jhu) <- c("date", "confirm")

# estimate only on 16 weeks of data due to computing times
reported_cases_latest_jhu <- reported_cases_jhu[reported_cases_jhu$date >= min(estimates_published$date) &
                                                  reported_cases_jhu$date <= max(estimates_published[
                                                    estimates_published$type != "forecast",]$date),]


# compare case data visually
plot(reported_cases_latest_jhu, type="l", col="red")
lines(reported_cases_latest)
legend(x="topleft", legend = c("JHU (globalrt)", "WHO (epiforecasts)"), lty=1, col = c("red", "black"))

# adjust mean and sd of generation time to the values used in for the globart estimates
generation_time_globalrt <- generation_time
generation_time_globalrt$mean <- generation_time_globalrt$sd <- 7
#generation_time_globalrt$mean_sd <- generation_time_globalrt$sd_sd <- 0

# do estimation
estimates_globalrt_params <- epinow(reported_cases = reported_cases_latest_jhu,
                                    generation_time = generation_time_globalrt,
                                    delays = delay_opts(incubation_period, reporting_delay),
                                    rt = rt_opts(prior = list(mean = 1, sd = 0.2)),
                                    stan = stan_opts(cores = 4),
                                    horizon = 14,
                                    CrIs = c(0.5, 0.95),
                                    verbose = TRUE)

# extract and save R estimates
epiforecasts_R_calc_globalrt <- estimates_globalrt_params$estimates$summarised[variable=="R", c("date", "type", "mean")]
qsave(epiforecasts_R_calc_globalrt, paste0("R_calc_", Sys.Date(), "_globalrtParams.qs"))

# compare epiforecast estimates to the ones obtained when using JHU data and globalrt parameters
plot(epiforecasts_R_calc_globalrt$date, epiforecasts_R_calc_globalrt$mean, type="l", col="red")
lines(epiforecasts_R_calc$date, epiforecasts_R_calc$mean)
legend(x="topleft", legend = c("globalrt parameters", "epiforecasts parameters"), lty=1, col = c("red", "black"))



### ETH

# load reported cases with function from covidregionaldata (complete history)
# use source = "jhu" to be in line with globalrt
reported_cases_eth <- get_national_data(countries = country, source = "covid19")
reported_cases_eth <- reported_cases_eth[, c("date", "cases_new")]
reported_cases_eth$cases_new <- c(0, diff(reported_cases_eth$cases_new, lag = 1))
names(reported_cases_eth) <- c("date", "confirm")

# estimate only on 16 weeks of data due to computing times
reported_cases_latest_eth <- reported_cases_eth[reported_cases_eth$date >= min(estimates_published$date) &
                                                  reported_cases_eth$date <= max(estimates_published[
                                                    estimates_published$type != "forecast",]$date),]


# compare case data visually
plot(reported_cases_latest_eth, type="l", col="red")
lines(reported_cases_latest)
legend(x="topleft", legend = c("Covid19DataHub/RKI (ETH)", "WHO (epiforecasts)"), lty=1, col = c("red", "black"))

generation_time_eth <- generation_time
generation_time_eth$mean <- 4.8
generation_time_eth$sd <- 2.3

# do estimation
estimates_eth_params <- epinow(reported_cases = reported_cases_latest_eth,
                               generation_time = generation_time_eth,
                               delays = delay_opts(incubation_period, reporting_delay),
                               rt = rt_opts(prior = list(mean = 1, sd = 0.2)),
                               stan = stan_opts(cores = 4),
                               horizon = 14,
                               CrIs = c(0.5, 0.95),
                               verbose = TRUE)

# compare published and calculated in plot
epiforecasts_R_calc_eth <- estimates_eth_params$estimates$summarised[variable=="R", c("date", "type", "median")]
qsave(epiforecasts_R_calc_eth, paste0("R_calc_", Sys.Date(), "_ETHParams.qs"))

# compare epiforecast estimates to the ones obtained when using JHU data and globalrt parameters
plot(epiforecasts_R_calc_eth$date, epiforecasts_R_calc_eth$mean, type="l", col="red")
lines(epiforecasts_R_calc$date, epiforecasts_R_calc$mean)
legend(x="topleft", legend = c("ETH parameters", "epiforecasts parameters"), lty=1, col = c("red", "black"))

