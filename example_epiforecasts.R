# packages
#install.packages(c("data.table", "remotes", "EpiNow2"))		
#remotes::install_github("epiforecasts/covidregionaldata")
library(data.table)
library(EpiNow2)
library(covidregionaldata)
library(dplyr)
library(lubridate)

wd <- "D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code"
setwd(wd)

# target country (must be present in ECDC data)
country <- "germany"

# set number of cores to use fitting the model
# no benefit on runtime if cores > chains which is set to 4 by default
options(mc.cores = 2)

# literature distributions - please reach out if there are others you think should be supported
generation_time <- get_generation_time(disease = "SARS-CoV-2", source = "ganyani")
incubation_period <- get_incubation_period(disease = "SARS-CoV-2", source = "lauer")
incubation_period$mean <- 5.2
incubation_period$mean_sd <- 1.1
incubation_period$sd <- 1.52
incubation_period$sd_sd <- 1.1

# define reporting delay as lognormal with mean of 4 days and sd of 1 day in absence of
# evidence. If data on onset -> report then can use estimate_delay to estimate the delay
# EB: with original values
reporting_delay <- list(mean = convert_to_logmean(4, 1),
                        mean_sd = 0.1,
                        sd = convert_to_logsd(4, 1),
                        sd_sd = 0.1,
                        max = 15)
# EB: with values from Abbot2020
reporting_delay <- list(mean = 6.5,
                        mean_sd = 0.1,
                        sd = 5, # EB: paper says 17, but this leads to introduction of NAs
                        sd_sd = 0.1,
                        max = 30)

#delays <- delay_opts(incubation_period, reporting_delay)

# precleaned ECDC data - alternative is to bring your own
reported_cases <- covidregionaldata::get_national_data(country, source = "ECDC")
reported_cases <- data.table::setDT(reported_cases)
reported_cases <- reported_cases[, .(date, confirm = cases_new)]
# filter to the last 3 months of data (to speed up computation)
reported_cases <- reported_cases[date >= (max(date) - 12*7)]

# EB: use same data source as used for published estimates for germany
rki_data <- read.csv("Rt_estimate_reconstruction/incidence_data/RKI_COVID19.csv")
rki_data <- mutate(rki_data, date = as_date(Refdatum))
rki_agg <- aggregate(rki_data$AnzahlFall, by=list(rki_data$date), FUN=sum)
names(rki_agg) <- c("date", "I")
plot(rki_agg, type="l")
#lines(reported_cases$date, reported_cases$confirm, col="blue")

reported_cases <- data.table(date=rki_agg$date, confirm=rki_agg$I)
reported_cases <- reported_cases[date >= (max(date) - 24*7)]

rm(rki_data)
rm(rki_agg)

# estimate Rt and nowcast/forecast cases by date of infection
# on a 4 core computer this example should take between 2 ~ 5 minutes to run
# to see saved logs view the dated logs folder
# to see saved results see the dated results folder
# some data sets may produce divergent transition warnings
# this is not unexpected and is usually fine (if dts < 1% of samples)
# but is an area of research as we seek to optimise the underlying model.
# If you have some prior knowledge of the scaling between observations and 
# reports see ?obs_opts for options.
# If you have some prior knowledge on the truncation in your data or multiple
# snapshots of data see ?trunc_opts for options.
# Note that the default settings may not be appropriate for your use case.
# Example configurations are here: https://epiforecasts.io/EpiNow2/dev/reference/estimate_infections.html
out <- epinow(reported_cases = reported_cases, 
              generation_time = generation_time,
              delays = delay_opts(incubation_period, reporting_delay),
              rt = rt_opts(prior = list(mean = 1, sd = 1)), # EB: same as in the paper instead of mean = 1.5, sd = 0.5
              # here we define the quality of the gaussian process approximation
              # if the fit to data appears poor try increasing basis_prop and
              # potentially the boundary_scale (see ?gp_opts for details)
              # though this will likely increase runtimes.
              gp = gp_opts(basis_prop = 0.2),
              # in some instances stan chains can get stuck when estimating in 
              # these instances consider using the future fitting mode by passing 
              # `future = TRUE, max_execution_time = 60 * 30` to stan_opts and calling 
              # `future::plan("multiprocess")` prior to running epinow this will time out
              # chains after 30 minutes but still return results from completed chains
              stan = stan_opts(),
              horizon = 14, 
              target_folder = "results",
              logs = file.path("logs", Sys.Date()),
              return_output = TRUE, 
              verbose = TRUE)


# arguments used in example
delays = delay_opts(incubation_period, reporting_delay)
truncation = trunc_opts()
rt = rt_opts(prior = list(mean = 1, sd = 1))
backcalc = backcalc_opts()
gp = gp_opts()
obs = obs_opts()
stan = stan_opts()
horizon = 14 
CrIs = c(0.5)
return_output = TRUE
output = c("samples", "plots", "latest", "fit", "timing")
target_folder = "results"
target_date <- max(reported_cases$date, na.rm = TRUE)
forecast_args = NULL
logs = file.path("logs", Sys.Date())
id = "epinow"
verbose = TRUE

futile.logger::flog.threshold(futile.logger::DEBUG, 
                              name = "EpiNow2.epinow")

setup_default_logging(logs = logs, target_date = target_date, 
                      mirror_epinow = TRUE)

match_output_arguments <- function (input_args = c(), supported_args = c(), logger = NULL, 
                                    level = "info") 
{
  if (level %in% "info") {
    flog_fn <- futile.logger::flog.info
  }
  else if (level %in% "debug") {
    flog_fn <- futile.logger::flog.debug
  }
  output_args <- rep(FALSE, length(supported_args))
  names(output_args) <- supported_args
  found_args <- lapply(input_args, function(arg) {
    supported_args[grepl(arg, supported_args)]
  })
  found_args <- unlist(found_args)
  found_args <- unique(found_args)
  if (!is.null(logger)) {
    if (length(found_args) > 0) {
      flog_fn("Producing following optional outputs: %s", 
              paste(found_args, collapse = ", "), name = logger)
    }
    else {
      flog_fn("No optional output specified", name = logger)
    }
  }
  output_args[names(output_args) %in% found_args] <- TRUE
  return(output_args)
}

output <- match_output_arguments(output,
                                 supported_args = c("plots", "samples", "fit", "timing", "latest"),
                                 logger = "EpiNow2.epinow", 
                                 level = "debug")

target_folders <- setup_target_folder(target_folder, target_date)
target_folder <- target_folders$date
latest_folder <- target_folders$latest

start_time <- Sys.time()

# epinow internal step-by-step
futile.logger::flog.threshold(futile.logger::DEBUG, 
                              name = "EpiNow2.epinow")
reported_cases <- setup_dt(reported_cases)
save_input(reported_cases, target_folder)
horizon <- update_horizon(horizon, target_date, reported_cases)

# critical
estimates <- estimate_infections(reported_cases = reported_cases, 
                                 generation_time = generation_time, delays = delays, 
                                 truncation = truncation, rt = rt, backcalc = backcalc, 
                                 gp = gp, obs = obs, stan = stan, CrIs = CrIs, horizon = horizon, 
                                 verbose = verbose, id = id)
if (!output["fit"]) {
  estimates$fit <- NULL
  estimates$args <- NULL
}
save_estimate_infections(estimates, target_folder, samples = output["samples"], 
                         return_fit = output["fit"])
if (!is.null(forecast_args)) {
  forecast <- do.call(forecast_infections, c(list(infections = estimates$summarised[variable == 
                                                                                      "infections"][type != "forecast"][, `:=`(type, 
                                                                                                                               NULL)], rts = estimates$summarised[variable == 
                                                                                                                                                                    "R"][type != "forecast"][, `:=`(type, NULL)], 
                                                  gt_mean = estimates$summarised[variable == "gt_mean"]$mean, 
                                                  gt_sd = estimates$summarised[variable == "gt_sd"]$mean, 
                                                  gt_max = generation_time$max, horizon = horizon, 
                                                  CrIs = CrIs), forecast_args))
  save_forecast_infections(forecast, target_folder, 
                           samples = output["samples"])
}
else {
  forecast <- NULL
}
estimated_reported_cases <- estimates_by_report_date(estimates, 
                                                     forecast, delays = delays, target_folder = target_folder, 
                                                     samples = output["samples"], CrIs = CrIs)
summary <- summary.estimate_infections(estimates, return_numeric = TRUE, 
                                       target_folder = target_folder)
if (output["plots"]) {
  plots <- plot.estimate_infections(estimates, type = "all", 
                                    target_folder = target_folder)
}
else {
  plots <- NULL
}
if (return_output) {
  out <- construct_output(estimates, forecast, estimated_reported_cases, 
                          plots = plots, summary, samples = output["samples"])
  return(out)
}
else {
  return(invisible(NULL))
}
end_time <- Sys.time()

if (!is.null(out$error)) {
  out$trace <- rlang::trace_back()
}
if (!is.null(target_folder) & !is.null(out$error)) {
  saveRDS(out$error, paste0(target_folder, "/error.rds"))
  saveRDS(out$trace, paste0(target_folder, "/trace.rds"))
}
if (output["timing"]) {
  out$timing <- round(as.numeric(end_time - start_time), 
                      1)
  if (!is.null(target_folder)) {
    saveRDS(out$timing, paste0(target_folder, "/runtime.rds"))
  }
}
if (output["latest"]) {
  copy_results_to_latest(target_folder, latest_folder)
}
if (return_output) {
  class(out) <- c("epinow", class(out))
  return(out)
}
else {
  return(invisible(NULL))
}


estimates <- estimate_infections(reported_cases = rc, 
                                 generation_time = generation_time,
                                 delays = delay_opts(incubation_period, reporting_delay), 
                                 truncation = trunc_opts(),
                                 rt = rt_opts(prior = list(mean = 1, sd = 1)),
                                 backcalc = backcalc_opts(), gp = gp_opts(),
                                 obs = obs_opts(), stan = stan_opts(),
                                 CrIs = c(0.5), horizon = horizon, 
                                 verbose = TRUE, id = "epinow")

reported_cases <- data.table::rbindlist(list(data.table::data.table(date = seq(min(reported_cases$date) - 
                                                                                 delays$seeding_time - backcalc$prior_window, min(reported_cases$date) - 
                                                                                 1, by = "days"), confirm = 0, breakpoint = 0), reported_cases))


# summary of the latest estimates
summary(out)
# plot estimates
plot(out)
# summary of R estimates
summary(out, type = "parameters", params = "R")
