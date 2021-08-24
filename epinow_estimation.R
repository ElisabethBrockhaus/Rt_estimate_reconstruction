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

# load reported cases preprocessed via https://github.com/epiforecasts/covidregionaldata/blob/master/R/Germany.R
#reported_cases_epiforecasts <- qread("reported_cases_epiforecasts.qs")
reported_cases_epiforecasts <- qread("reported_cases_aggregated.qs")
reported_cases_epiforecasts$region <- country

reported_cases <- read_csv("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/master/national/cases/summary/reported_cases.csv")
reported_cases <- reported_cases[reported_cases$region==country,]
reported_cases <- na.omit(reported_cases)

# load distributions that were published at Github (https://github.com/seabbs/covid-climate-rt/tree/master/delays/data)
generation_time <- readRDS("Rt_estimate_reconstruction/distributions/generation_time.rds")

incubation_period <- readRDS("Rt_estimate_reconstruction/distributions/incubation_period.rds")
#reporting_delay <- readRDS("distributions/onset_to_report_delay.rds") # wrong repo
reporting_delay <- readRDS("Rt_estimate_reconstruction/distributions/onset_to_admission_delay.rds")

# do estimation
estimates <- regional_epinow(#reported_cases = reported_cases_epiforecasts[reported_cases_epiforecasts$date >= "2021-04-01", ], 
                             reported_cases = reported_cases,
                             generation_time = generation_time,
                             delays = delay_opts(incubation_period, reporting_delay),
                             rt = rt_opts(prior = list(mean = 1, sd = 0.2)),
                             stan = stan_opts(cores = 4),
                             #horizon = 14,
                             verbose = TRUE)

# summary of the latest estimates
summary(estimates$regional$Germany)

# plot estimates
plot(estimates$regional$Germany)

# extract R estimates
EpiNow2_est <- summary(estimates$regional$Germany, type = "parameters", params = "R")[,c("date", "median")]
qsave(EpiNow2_est, "EpiNow2_est.qs")

# load published estimates
estimates_published <- read_csv("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/master/national/cases/summary/rt.csv")
estimates_published <- estimates_published[estimates_published$country=="Germany",]

# compare published and calculated in plot
plot(EpiNow2_est$date, EpiNow2_est$median, type="l")
lines(estimates_published$date, estimates_published$median, col="red")


#########
# different versions of the reported case data

# load reported cases (aggregated from RKI-COVID19.csv)
reported_cases <- qread("reported_cases.qs")

# precleaned ECDC data - alternative is to bring your own
reported_cases_ecdc <- covidregionaldata::get_national_data(country, source = "ECDC")
reported_cases_ecdc <- data.table::setDT(reported_cases_ecdc)
reported_cases_ecdc <- reported_cases_ecdc[, .(date, confirm = cases_new)]

# load reported cases from Github
reported_cases_all <- read_csv("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/master/national/cases/summary/reported_cases.csv")
reported_cases_DE <- na.omit(reported_cases_all[reported_cases_all$region==country, c("date", "confirm")])

plot(seq(as.Date("2020-03-01"), max(reported_cases_DE$date), by="day"), rep(1, 535), type="l", col="white",
     ylim=c(0,29000), xlab="Date", ylab="I")
lines(reported_cases$date, reported_cases$confirm)
lines(reported_cases_DE$date, reported_cases_DE$confirm, col = "blue")
lines(reported_cases_ecdc$date, reported_cases_ecdc$confirm, col = "red")
lines(reported_cases_epiforecasts$date, reported_cases_epiforecasts$confirm, col="darkgreen")

