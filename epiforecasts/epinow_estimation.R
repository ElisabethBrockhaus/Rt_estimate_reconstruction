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

#######################
# original estimation #
#######################

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
incid_latest <- incid[incid$date >= as.Date("2021-04-01") &
                        incid$date <= as.Date("2021-09-30"),]
names(incid_latest) <- c("date", "confirm")

# parameter combinations used in papers
mean_gt <-           c(4.8, 4, 5.6, 4.8, 5, 3.4, 3.6,  4.7, 7)
sd_gt <-             c(2.3, 0, 4.2, 2.3, 4, 1.8, 3.1,  2.9, 7)
o <- 0.1 # replacement for zeros that would lead to mathematiccal issues
m <- 4 # replacement for zeros in mean incubation period (shifted back manually after the estimation)
mean_incubation <-   c(5.3, m, m,   m,   m, m,   5.5,  5,   m)
sd_incubation <-     c(3.2, o, o,   o,   o, o,   2.4,  o,   o)
mean_report_delay <- c(5.5, 1, 7,   10,  o, o,   6.5,  7.1, o)
sd_report_delay <-   c(3.8, o, o,   o,   o, o,   17.1, 5.9, o)

params <- data.frame(gt_mean=mean_gt, gt_sd=sd_gt,
                     incubation_mean=mean_incubation, incubation_sd=sd_incubation,
                     report_delay_mean=mean_report_delay, report_delay_sd=sd_report_delay)
methods <- c("ETH", "RKI", "Ilmenau", "SDSC", "Zi", "AGES", "epiforecasts", "rtlive", "globalrt")
rownames(params) <- methods

for (method in methods){
  
  print(paste0("Start with estimation using parameters from ", method))
  
  # set parameters
  generation_time <- list(
    mean = params[method, "gt_mean"], mean_sd = 0.1,
    sd = params[method, "gt_sd"], sd_sd = 0.1,
    max = 30
  )
  
  incubation_period <- list(
    mean = convert_to_logmean(params[method, "incubation_mean"],
                              params[method, "incubation_sd"]),
    mean_sd = 0.1,
    sd = convert_to_logsd(params[method, "incubation_mean"],
                          params[method, "incubation_sd"]),
    sd_sd = 0.1,
    max = 30
  )
  
  reporting_delay <- list(
    mean = convert_to_logmean(params[method, "report_delay_mean"],
                              params[method, "report_delay_sd"]),
    mean_sd = 0.1,
    sd = convert_to_logsd(params[method, "report_delay_mean"],
                          params[method, "report_delay_sd"]),
    sd_sd = 0.1,
    max = 30
  )
  
  start_time <- Sys.time()
  
  # do estimation
  estimates_params <- epinow(reported_cases = incid_latest,
                             generation_time = generation_time,
                             delays = delay_opts(incubation_period, reporting_delay),
                             rt = rt_opts(prior = list(mean = 1, sd = 0.2)),
                             stan = stan_opts(cores = 4),
                             horizon = 14,
                             CrIs = c(0.5, 0.95),
                             verbose = TRUE)
  
  end_time <- Sys.time()
  print(paste0("Estimation took ", round(difftime(end_time, start_time, units='mins'), digits=2), " minutes."))
  
  plot(estimates_params)
  
  # compare published and calculated in plot
  epiforecasts_R_calc <- estimates_params$estimates$summarised[variable=="R",
                                                               c("date", "type", "median", "mean", "sd",
                                                                 "lower_95", "lower_50", "upper_95", "upper_50")]
  
  if (params[method, "incubation_mean"] == 4){
    # included 4 days too much delay
    # shift estimates forwar by 4 days
    epiforecasts_R_calc$date <- epiforecasts_R_calc$date + 4
  }
  
  qsave(epiforecasts_R_calc, paste0("R_calc_", Sys.Date(), method, "Params.qs"))
  
  plot(epiforecasts_R_calc$date, epiforecasts_R_calc$mean, type="l")
}
