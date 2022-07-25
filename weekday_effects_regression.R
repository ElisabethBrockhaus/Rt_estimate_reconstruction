library(lubridate)
library(dplyr)
library(splines)

Sys.setlocale("LC_TIME", "English")

setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
getwd()



################################
# get functions and parameters #
################################

source("Rt_estimate_reconstruction/load_data.R")

pub_delays <- read.csv("Rt_estimate_reconstruction/otherFiles/pub_delays.csv", row.names = 1)

methods <- c("Braunschweig", "epiforecasts", "ETHZ_sliding_window",
             "globalrt_7d", "ilmenau", "RKI_7day", "rtlive", "SDSC")

method <- methods[6]
country <- "DE"
start_date <- as_date("2020-11-16")
end_date <- as_date("2021-05-01")
conf_level <- "95"
path_estimates <- "reproductive_numbers/data-processed/"

min_lag <- pub_delays[method, country]
max_lag <- min_lag + 6

pub_dates <- list.files(paste0(path_estimates, method),
                        pattern = "\\d{4}-\\d{2}-\\d{2}",
                        full.names = F) %>% substr(1, 10)
pub_dates <- pub_dates[which(as_date(pub_dates) <= end_date &
                               as_date(pub_dates) >= start_date - max_lag)]

R_est_ts <<- data.frame(date = seq(min(as_date(pub_dates)) - max_lag,
                                   max(as_date(pub_dates)) - min_lag,
                                   by = "day"))

# read 7 most recent R values for all estimation dates
for (pub_date in pub_dates){
  tryCatch(
    {
      R_est <- load_published_R_estimates(method,
                                          start = as_date(pub_date) - max_lag,
                                          end = as_date(pub_date) - min_lag,
                                          pub_date = pub_date,
                                          location = country,
                                          conf_level = conf_level,
                                          verbose = F) %>%
        dplyr::select(date, R_pub)
      
      R_est <- R_est %>%
        rename(R = R_pub) %>%
        rename_with(~ paste0(.x, "_estimated_", pub_date), !date)
      
    },
    error = function(e) {R_est <<- data.frame(date = seq(as_date(pub_date) - max_lag,
                                                         as_date(pub_date) - min_lag,
                                                         by = "day"))}
  )
  R_est_ts <- R_est_ts %>% full_join(R_est, by="date")
}
  
R_est_ts <- R_est_ts %>% rename(target_date = date)
num_est_dates <- dim(R_est_ts)[2]-1
num_target_dates <- dim(R_est_ts)[1]

# create filter for existing values in R_est_ts
non_na_entries <- !(R_est_ts %>% dplyr::select(-target_date) %>% is.na())
  
#### first regressor: splines for target date
# create splines of target_date
target_dates <- R_est_ts %>%
  transmute(across(!target_date, list(td = ~ {target_date})))

splines_td <- bs(as_date(target_dates[non_na_entries]))

#### second regressor: weekday of target date
weekdays_td <- data.frame(lapply(target_dates, weekdays))[non_na_entries]
weekdays_td <- relevel(as.factor(weekdays_td), ref = "Monday")

#### third regressor: splines for difference between estimation date and target date
time_diffs <- R_est_ts %>%
  # estimation date - target date
  transmute(across(!target_date, list(lag = ~ {as.numeric(difftime(as_date(substr(cur_column(), 13, 22)), target_date), units = "days")})))
splines_time_diff <- bs(time_diffs[non_na_entries])

#### Dependent variable: log(R(target_date, estimation_date))
R <- (R_est_ts %>% dplyr::select(-target_date))[non_na_entries]
log_R <- log(R)


####################
# fit linear model #
####################

model <- lm(log_R ~ splines_td + weekdays_td + splines_time_diff)
summary(model)
  
weekday_effects <- c(0, coef(model)[c("weekdays_tdTuesday", "weekdays_tdWednesday",
                                      "weekdays_tdThursday", "weekdays_tdFriday",
                                      "weekdays_tdSaturday", "weekdays_tdSunday")])
names(weekday_effects) <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
plot(weekday_effects, type="l")




