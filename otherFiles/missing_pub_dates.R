library(lubridate)
library(dplyr)

setwd("../..")
getwd()

source("Rt_estimate_reconstruction/load_data.R")

path_estimates <- "reproductive_numbers/data-processed/"
methods <- c("RKI_7day", "ETHZ_sliding_window", "ilmenau", "sdsc", "globalrt_7d", "epiforecasts", "rtlive", "Braunschweig")

start_date <- as_date("2020-11-16")
end_date <- as_date("2021-05-01")
target_dates <- seq(start_date, end_date, by = "day")

df_missing_dates <- data.frame(matrix(rep(NA, length(methods)*length(target_dates)),
                                      ncol = length(methods)),
                               row.names = target_dates) %>%
  setNames(methods)

for (method in methods){
  files <- list.files(paste0(path_estimates, method))
  dates <- files %>% substr(1,10) %>% as_date
  dates <- dates[dates >= start_date & dates <= end_date]
  dates_missing <- target_dates %>%
    lubridate::setdiff(dates)
}
