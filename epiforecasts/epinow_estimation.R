library(data.table)
library(EpiNow2)
library(covidregionaldata)
library(dplyr)
library(lubridate)
library(qs)
library(readr)

setwd("/home/brockhaus/reproductive_numbers")

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
qsave(epiforecasts_R_calc, paste0("R_calc_", date_of_data, ".qs"))

# compare published and calculated in plot
plot(epiforecasts_R_calc$date, epiforecasts_R_calc$mean, type="l")
lines(estimates_published$date, estimates_published$mean, col="red")

plot(epiforecasts_R_calc$date,
     epiforecasts_R_calc$mean - estimates_published[estimates_published$date %in% epiforecasts_R_calc$date,]$mean,
     type="l")
abline(h=0, col="grey", lty=2)



###########################
# function for estimation #
###########################

calculacte_and_save_estimates <- function(reported_cases,
                                          generation_time,
                                          incubation_period,
                                          reporting_delay,
                                          method,
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
  
  # for originally constant delays correct for different mean used in epiforecasts estimation
  result$date <- result$date + params[method, "necessary_shift"]
  
  qsave(result, paste0("R_calc_", date_of_data, "_", method, variation, ".qs"))
  
  plot(result$date, result$mean, type="l", main=method, xlab="date", ylab="Rt estimate")
  
  # use result for other parameter sources if they differ only concerning a constant delay
  if ((method == "globalrt") & (variation == "_delays")){
    # reverse globalrt delay shift
    result$date <- result$date - params[method, "delay_shift"]

    for (method in c("RKI", "Ilmenau", "SDSC", "Zi", "AGES")){
      # shift as needed for this parameter sources
      result$date <- result$date + params[method, "delay_shift"]
      
      # save result
      qsave(result, paste0("R_calc_", date_of_data, "_", method, variation, ".qs"))
      
      # reverse previous shift for next method
      result$date <- result$date - params[method, "delay_shift"]
    }
  }
  
  # save result under delays/gtd only names if inherent parameters are used
  if (method == "epiforecasts"){
    qsave(result, paste0("R_calc_", date_of_data, "_", method, "_delays", ".qs"))
    qsave(result, paste0("R_calc_", date_of_data, "_", method, "_GTD", ".qs"))
  }
}



############################
# use different parameters #
############################

# load incidence data from RKI
#repo_incid <- "https://raw.githubusercontent.com/robert-koch-institut/SARS-CoV-2-Nowcasting_und_-R-Schaetzung/main/Nowcast_R_aktuell.csv"
#data <- read_csv(repo_incid)
#names(data) <- c("Datum", "NeuErkr", "lb_NeuErkr", "ub_NeuErkr",
#                 "NeuErkr_ma4", "lb_NeuErkr_ma4", "ub_NeuErkr_ma4",
#                 "R_7Tage", "lb_R_7Tage", "ub_R_7Tage")
#incid <- data.frame(date=data$Datum, confirm=data$NeuErkr)
#incid <- incid[incid$date < "2021-10-01",]
#rm(data)

# load incidence data as used by rtlive (RKI line list aggregated)
incid <- read_csv("rtlive_incid.csv")
names(incid) <- c("date", "confirm")
incid <- incid[incid$date < "2021-10-01",]
incid <- incid[incid$date >= "2021-01-01",]

# set date for file names to the date, when the data was loaded
date_of_data <- "2021-10-19"

# parameter combinations used in papers
o <- 0.1 # replacement for zeros that would lead to mathematical issues
oo <- 0.5
m <- 2 # mean incubation period and reporting delay in case of constant delay distributions (shifted back manually after the estimation)
mean_gt <-        c(4.8,  4,   5.6, 4.8, 5,   3.4, 3.6, 4.7,  7)
sd_gt <-          c(2.3,  0,   4.2, 2.3, 4,   1.8, 3.1, 2.9,  7)
mean_delay <-     c(10.8, 1,   7,   10,  0,   0,   12,  12.1, 0)
necessary_shift <- c(-1.2, -11, -5,  -2,  -12, -12, 0,   0.1,  0) # always using epiforecasts delay distributions with mean 12 -> correct through shift

params <- data.frame(gt_mean=mean_gt, gt_sd=sd_gt,
                     necessary_shift=round(necessary_shift))
methods <- c("ETH", "RKI", "Ilmenau", "SDSC", "Zi", "AGES", "epiforecasts", "rtlive", "globalrt")
rownames(params) <- methods

# estimation with adjusted input data using epiforecasts delay and generation time distributions
method <- "epiforecasts"
incubation_period <- readRDS("distributions/incubation_period.rds")
reporting_delay <- readRDS("distributions/onset_to_admission_delay.rds")
generation_time <- readRDS("distributions/generation_time.rds")

print("Start with estimation using original distributions")

calculacte_and_save_estimates(incid, generation_time,
                              incubation_period, reporting_delay,
                              method, variation = "")

# estimation with adjusted input data and generation time distribution using epiforecasts delay distributions
method <- "globalrt"
generation_time <- list(
  mean = 4, mean_sd = 0.5,
  sd = 4, sd_sd = 0.5,
  max = 30
)

print(paste0("Start with estimation using generation time from", method))

calculacte_and_save_estimates(incid, generation_time,
                              incubation_period, reporting_delay,
                              method, variation = "_GTD")



# estimation with adjusted delays and generation time
for (method in c("globalrt")){ #, "ETH", "epiforecasts", "rtlive", "SDSC", "Ilmenau", "Zi", "AGES")){
  
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
  
  calculacte_and_save_estimates(incid, generation_time,
                                incubation_period, reporting_delay,
                                method, variation = "Params")
}


# estimation with adjusted delays only, using epiforecast generation time distribution
generation_time <- readRDS("distributions/generation_time.rds")

for (method in c("globalrt")){ #, "ETH", "rtlive")){
  
  print(paste0("Start with estimation using delays from ", method))
  
  # set delay parameters
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
  
  calculacte_and_save_estimates(incid, generation_time,
                                incubation_period, reporting_delay,
                                method, variation = "_delays")
}


# estimation with adjusted generation time only, using epiforecast delay distributions
incubation_period <- readRDS("distributions/incubation_period.rds")
reporting_delay <- readRDS("distributions/onset_to_admission_delay.rds")

for (method in c("globalrt")){ #, "ETH", "rtlive", "SDSC", "Ilmenau", "Zi", "AGES")){
  
  print(paste0("Start with estimation using generation time from ", method))
  
  # set generation time
  generation_time <- list(
    mean = params[method, "gt_mean"], mean_sd = 0.1,
    sd = params[method, "gt_sd"], sd_sd = 0.1,
    max = 30
  )

  calculacte_and_save_estimates(incid, generation_time,
                                incubation_period, reporting_delay,
                                method, variation = "_GTD")
}
