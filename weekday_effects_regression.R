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

wds <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
weekday_effects_target_date <- data.frame(matrix(rep(NA, length(methods)*7), nrow=7))
colnames(weekday_effects_target_date) <- methods
row.names(weekday_effects_target_date) <- wds
weekday_effects_estimation_date <- weekday_effects_target_date

country <- "DE"
start_date <- as_date("2020-11-16")
end_date <- as_date("2021-05-01")
path_estimates <- "reproductive_numbers/data-processed/"

# choose whether to encode (estimation date - target date) as splines or dummies
splines_for_difftime <- F

for (method in methods){
  print(method)
  
  min_lag <- pub_delays[method, country]
  max_lag <- min_lag + 6
  
  pub_dates <- list.files(paste0(path_estimates, method),
                          pattern = "\\d{4}-\\d{2}-\\d{2}",
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which(as_date(pub_dates) <= end_date &
                                 as_date(pub_dates) >= start_date - max_lag)]
  
  R_est_ts <- data.frame(date = seq(min(as_date(pub_dates)) - max_lag,
                                    max(as_date(pub_dates)) - min_lag,
                                    by = "day"))
  
  # read 7 most recent R values for all estimation dates (ed)
  for (pub_date in pub_dates){
    tryCatch(
      {
        R_est <- load_published_R_estimates(method,
                                            start = as_date(pub_date) - max_lag,
                                            end = as_date(pub_date) - min_lag,
                                            pub_date = pub_date,
                                            location = country,
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
  
  #### first regressor: splines for target date (td)
  # create splines of target_date
  target_dates <- R_est_ts %>%
    transmute(across(!target_date, list(td = ~ {target_date})))
  
  splines_td <- bs(as_date(target_dates[non_na_entries]))
  
  #### second regressor: weekday of target date
  weekdays_td <- data.frame(lapply(target_dates, weekdays))[non_na_entries]
  weekdays_td <- relevel(as.factor(weekdays_td), ref = "Monday") # define baseline category
  
  #### third regressor: splines or dummies for difference between estimation date and target date
  time_diffs <- R_est_ts %>%
    # estimation date - target date
    transmute(across(!target_date, list(lag = ~ {as.numeric(difftime(as_date(substr(cur_column(), 13, 22)), target_date), units = "days")})))
  if(splines_for_difftime) {
    time_diff <- bs(time_diffs[non_na_entries])
  } else {
    time_diff <- as.factor(time_diffs[non_na_entries])
  }
  
  #### optional fourth regressor: weekday of estimation date
  estimation_dates <- R_est_ts %>%
    transmute(across(!target_date, list(td = ~ {as_date(substr(cur_column(), 13, 22))})))
  weekdays_ed <- data.frame(lapply(estimation_dates, weekdays))[non_na_entries]
  weekdays_ed <- relevel(as.factor(weekdays_ed), ref = "Monday")
  
  #### Dependent variable: log(R(target_date, estimation_date))
  R <- (R_est_ts %>% dplyr::select(-target_date))[non_na_entries]
  log_R <- log(R)
  
  
  ####################
  # fit linear model #
  ####################
  
  model <- lm(log_R ~ splines_td + weekdays_td + time_diff + weekdays_ed)
  print(summary(model))
  
  wd_td_regressors <- paste0("weekdays_td", wds[2:7])
  wd_ed_regressors <- paste0("weekdays_ed", wds[2:7])
  
  weekday_effects_target_date[,method] <- c(0, coef(model)[wd_td_regressors])
  weekday_effects_estimation_date[,method] <- c(0, coef(model)[wd_ed_regressors])
}

write.csv(weekday_effects_target_date, "Rt_estimate_reconstruction/otherFiles/weekday_effects_target_date.csv")
write.csv(weekday_effects_estimation_date, "Rt_estimate_reconstruction/otherFiles/weekday_effects_estimation_date.csv")

source("Rt_estimate_reconstruction/prepared_plots.R")

for (date_type in c("target", "estimation")){
  data <- read.csv(paste0("Rt_estimate_reconstruction/otherFiles/weekday_effects_", date_type, "_date.csv"), row.names = 1)
  
  plot_data <- data %>%
    rownames_to_column("weekday") %>%
    gather("method", "value", 2:(length(methods)+1)) %>%
    mutate(method = plyr::mapvalues(method,
                                    c("Braunschweig", "ETHZ_sliding_window", "globalrt_7d", "ilmenau", "RKI_7day"),
                                    c("HZI",          "ETH",                 "globalrt",    "Ilmenau", "RKI"))) %>%
    arrange(method) %>%
    mutate(exp_value = exp(value))
  
  methods_legend <- unique(plot_data$method)
  col_values <- get_colors(methods = methods_legend, palette = "methods")
  
  plot_raw <- ggplot() +
    theme_minimal() +
    theme(
      plot.margin = unit(c(3,14,2,3), "mm"),
      plot.title = element_text(size=18),
      axis.text = element_text(size=16),
      axis.title.y = element_text(size=18),
      axis.title.x = element_blank(),
      legend.text = element_text(size=16),
      legend.title = element_text(size=18),
      legend.position = "bottom",
      panel.border = element_rect(fill = "transparent", size = 0.5),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major = element_line(),
      panel.grid.minor = element_blank()
    )
  
  plot <- plot_raw +
    geom_line(data=plot_data,
              aes(x = factor(weekday, levels = wds),
                  y = value,
                  group = method,
                  color = method),
              size = .8, na.rm = T) +
    scale_x_discrete() + 
    scale_color_manual(values = col_values, name = "method") +
    ylab(paste("coefficient for", date_type, "date weekday"))
  
  print(plot)
  
  plot <- plot_raw +
    geom_line(data=plot_data,
              aes(x = factor(weekday, levels = wds),
                  y = exp_value,
                  group = method,
                  color = method),
              size = .8, na.rm = T) +
    scale_x_discrete() + 
    scale_color_manual(values = col_values, name = "method") +
    ylab(paste("exp(coefficient) for", date_type, "date weekday"))
  
  print(plot)
}



