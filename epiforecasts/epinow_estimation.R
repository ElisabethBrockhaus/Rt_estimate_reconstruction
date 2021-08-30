library(data.table)
library(EpiNow2)
library(covidregionaldata)
library(dplyr)
library(lubridate)
library(qs)
library(readr)

wd <- "D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code"
setwd(wd)

country <- "Germany"

# load reported cases from github
reported_cases <- read_csv("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/master/national/cases/summary/reported_cases.csv")
reported_cases <- reported_cases[reported_cases$region==country,]
reported_cases <- na.omit(reported_cases)

# load distributions that were published at Github (https://github.com/seabbs/covid-climate-rt/tree/master/delays/data)
generation_time <- readRDS("Rt_estimate_reconstruction/distributions/generation_time.rds")
incubation_period <- readRDS("Rt_estimate_reconstruction/distributions/incubation_period.rds")
reporting_delay <- readRDS("Rt_estimate_reconstruction/distributions/onset_to_admission_delay.rds")

# do estimation
estimates <- regional_epinow(reported_cases = reported_cases,
                             generation_time = generation_time,
                             delays = delay_opts(incubation_period, reporting_delay),
                             rt = rt_opts(prior = list(mean = 1, sd = 0.2)),
                             stan = stan_opts(cores = 4),
                             horizon = 14,
                             verbose = TRUE)

# plot estimates
plot(estimates$regional[[country]])

# extract R estimates
EpiNow2_est <- summary(estimates$regional[[country]],
                       type = "parameters", params = "R")[,c("date", "median")]
qsave(EpiNow2_est, "EpiNow2_est.qs")

# load published estimates
estimates_published <- read_csv("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/master/national/cases/summary/rt.csv")
estimates_published <- estimates_published[estimates_published$country==country,]

# compare published and calculated in plot
plot(EpiNow2_est$date, EpiNow2_est$median, type="l")
lines(estimates_published$date, estimates_published$median, col="red")

