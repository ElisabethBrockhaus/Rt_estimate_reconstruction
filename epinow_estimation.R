library(data.table)
library(EpiNow2)
library(covidregionaldata)
library(dplyr)
library(lubridate)
library(qs)
library(readr)

wd <- "D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code"
setwd(wd)

# load reported cases (aggregated from RKI-COVID19.csv)
reported_cases <- qread("Rt_estimate_reconstruction/incidence_data/reported_cases_epiforecasts.qs")
subnational_data_DE <- read_csv("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/master/subnational/germany/cases/summary/reported_cases.csv")
data_aggregated <- aggregate(subnational_data_DE$confirm, by=list(subnational_data_DE$date), FUN=sum)
names(data_aggregated) <- c("date", "confirm")
qsave(data_aggregated, "Rt_estimate_reconstruction/incidence_data/reported_cases_aggregated.qs")


plot(reported_cases, type="l")
lines(data_aggregated, col="red")
plot(data_aggregated$date,
     reported_cases[reported_cases$date %in% data_aggregated$date, ]$confirm -
       data_aggregated$confirm,
     type="l")

cases_from_github <- read_csv("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/master/national/cases/summary/reported_cases.csv")
cases_from_github <- read_csv("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/master/national/cases/summary/cases_by_report.csv",
                              col_select = c("country", "date", "median", "mean"))
cases_from_github <- cases_from_github[cases_from_github$region=="Germany",]
cases_from_github <- na.omit(cases_from_github)

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

