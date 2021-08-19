library(data.table)
library(EpiNow2)
library(covidregionaldata)
library(dplyr)
library(lubridate)
library(qs)

wd <- "D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code"
setwd(wd)

reported_cases <- qread("Rt_estimate_reconstruction/incidence_data/reported_cases.qs")

generation_time <- readRDS("Rt_estimate_reconstruction/distributions/generation_time.rds")
incubation_period <- readRDS("Rt_estimate_reconstruction/distributions/incubation_period.rds")
reporting_delay <- readRDS("Rt_estimate_reconstruction/distributions/onset_to_report_delay.rds")

estimates <- epinow(reported_cases = reported_cases, 
                    generation_time = generation_time,
                    delays = delay_opts(incubation_period, reporting_delay),
                    rt = rt_opts(prior = list(mean = 2, sd = 0.2)),
                    stan = stan_opts(cores = 4),
                    verbose = TRUE)