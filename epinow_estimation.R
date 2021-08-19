library(data.table)
library(EpiNow2)
library(covidregionaldata)
library(dplyr)
library(lubridate)
library(qs)

wd <- "D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code"
setwd(wd)

# load reported cases (aggregated from RKI-COVID19.csv)
reported_cases <- qread("Rt_estimate_reconstruction/incidence_data/reported_cases_epiforecasts.qs")

# load distributions that were published at Github (https://github.com/seabbs/covid-climate-rt/tree/master/delays/data)
generation_time <- readRDS("Rt_estimate_reconstruction/distributions/generation_time.rds")
incubation_period <- readRDS("Rt_estimate_reconstruction/distributions/incubation_period.rds")
reporting_delay <- readRDS("Rt_estimate_reconstruction/distributions/onset_to_report_delay.rds")

# do estimation
estimates <- epinow(reported_cases = reported_cases[reported_cases$date >= "2021-04-01", ], 
                    generation_time = generation_time,
                    delays = delay_opts(incubation_period, reporting_delay),
                    rt = rt_opts(prior = list(mean = 1, sd = 1)),
                    stan = stan_opts(cores = 4),
                    verbose = TRUE)

# summary of the latest estimates
summary(estimates)

# plot estimates
plot(estimates)

# extract R estimates
EpiNow2_est <- summary(estimates, type = "parameters", params = "R")[,c("date", "median")]

