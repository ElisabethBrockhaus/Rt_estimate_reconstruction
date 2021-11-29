library(data.table)
library(EpiNow2)
library(covidregionaldata)
library(dplyr)
library(lubridate)
library(qs)
library(readr)

#setwd("/home/brockhaus/reproductive_numbers")
setwd("D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code")
output_path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"


#######################
# load incidence data #
#######################

# load incidence data as used by rtlive (RKI line list aggregated)
incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_07_10.csv")
names(incid) <- c("date", "confirm")
incid <- incid[incid$date >= "2020-12-01",]

# set date for file names to the date, when the data was loaded
date_of_data <- "2021-07-10"

###########################
# function for estimation #
###########################

calculacte_and_save_estimates <- function(reported_cases,
                                          generation_time,
                                          incubation_period,
                                          reporting_delay,
                                          variation){
  start_time <- Sys.time()
  print(start_time)
  
  # do estimation
  estimates <- epinow(reported_cases = reported_cases,
                      generation_time = generation_time,
                      delays = delay_opts(incubation_period, reporting_delay),
                      rt = rt_opts(prior = list(mean = 1, sd = 0.2)),
                      stan = stan_opts(cores = 4),
                      horizon = 14,
                      CrIs = c(0.5, 0.95),
                      verbose = FALSE)
  
  end_time <- Sys.time()
  print(paste0("Estimation took ", round(difftime(end_time, start_time, units='mins'), digits=2), " minutes."))
  
  # compare published and calculated in plot
  result <- estimates$estimates$summarised[variable=="R",
                                           c("date", "type", "median", "mean", "sd",
                                             "lower_95", "lower_50", "upper_95", "upper_50")]
  
  qsave(result, paste0(output_path, "R_calc_", date_of_data, "_final_", variation, ".qs"))
  #plot(result$date, result$mean, type="l", main=method, xlab="date", ylab="Rt estimate")
  return(result)
}



############################
# use different parameters #
############################

# estimation with adjusted input data
path <- "Rt_estimate_reconstruction/epiforecasts/distributions/"
incubation_period <- readRDS(paste0(path, "incubation_period.rds"))
reporting_delay <- readRDS(paste0(path, "onset_to_admission_delay.rds"))
generation_time <- readRDS(paste0(path, "generation_time.rds"))

print("Start with estimation with adjusted input data")
R_epiforecasts_adjInput <- calculacte_and_save_estimates(incid, generation_time,
                                                         incubation_period, reporting_delay,
                                                         variation = "adjInput")

# no window size to adjust, save the same estimates for the second level of adjustment
qsave(R_epiforecasts_adjInput, paste0(output_path, "R_calc_", date_of_data, "_final_adjInputWindow.qs"))

# estimation with adjusted input data and generation time distribution
generation_time$mean <- 4
generation_time$sd <- 4

print("Start with estimation with adjusted input data and generation time distribution")
R_epiforecasts_adjInputWindowGTD <- calculacte_and_save_estimates(incid, generation_time,
                                                                  incubation_period, reporting_delay,
                                                                  variation = "adjInputWindowGTD")

# additionally adjust mean delay
# first calculate mean delay considered in the estimation
incubation <- function(x) dlnorm(x, meanlog = incubation_period$mean, sdlog = incubation_period$sd)
reporting <- function(y) dlnorm(y, meanlog = reporting_delay$mean, sdlog = reporting_delay$sd)
# convolution integral
delay <- function(z) integrate(function(x, z) reporting(z-x)*incubation(x),-Inf,Inf,z)$value
delay <- Vectorize(delay)
# mean of resulting convolution
z <- 0:10e4
mean <- z %*% delay(z)
mean
# shift estimates by mean delay
R_epiforecasts_final <- R_epiforecasts_adjInputWindowGTD
R_epiforecasts_final$date <- R_epiforecasts_final$date + round(mean)

# no window size to adjust, save the same estimates for the second level of adjustment
qsave(R_epiforecasts_final, paste0(output_path, "R_calc_", date_of_data, "_final_adjAll.qs"))

