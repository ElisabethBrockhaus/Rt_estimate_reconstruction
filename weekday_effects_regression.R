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

path_estimates <- "reproductive_numbers/data-processed/"

country <- "DE"
start_date <- as_date("2020-11-16")
end_date <- as_date("2021-05-01")
number_days <- 28

wds <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
weekday_effects_target_date <- data.frame(matrix(rep(NA, length(methods)*3*7), nrow=7))
colnames(weekday_effects_target_date) <- c(outer(methods, c("_coef", "_l", "_u"), FUN=paste0))
row.names(weekday_effects_target_date) <- wds
weekday_effects_estimation_date <- weekday_effects_target_date
time_diff_effects <- data.frame(matrix(rep(NA, 3*length(methods)*
                                             (max(pub_delays[,"DE"], na.rm=T) +
                                                number_days + 1 -
                                                min(pub_delays[,"DE"], na.rm=T))),
                                       ncol = length(methods)*3))
colnames(time_diff_effects) <- c(outer(methods, c("_coef", "_l", "_u"), FUN=paste0))
row.names(time_diff_effects) <- min(pub_delays[,"DE"], na.rm=T) : (max(pub_delays[,"DE"], na.rm=T) + number_days)

for (method in methods){
  print(method)
  
  min_lag <- pub_delays[method, country]
  max_lag <- min_lag + number_days - 1
  
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
  
  # check if splines are smooth
  plot(as_date(target_dates[non_na_entries]), splines_td[,1], type="l",
       ylab=NA, xlab=NA, main=paste0("Splines for target date (", method, ")"),
       ylim=c(0,1))
  lines(as_date(target_dates[non_na_entries]), splines_td[,2], col="blue")
  lines(as_date(target_dates[non_na_entries]), splines_td[,3], col="darkred")

  #### second regressor: weekday of target date
  weekdays_td <- data.frame(lapply(target_dates, weekdays))[non_na_entries]
  weekdays_td <- relevel(as.factor(weekdays_td), ref = "Monday") # define baseline category
  
  #### third regressor: splines or dummies for difference between estimation date and target date
  time_diffs <- R_est_ts %>%
    # estimation date - target date
    transmute(across(!target_date,
                     list(lag = ~ {as.numeric(difftime(as_date(substr(cur_column(), 13, 22)),
                                                       target_date),
                                              units = "days")})))
  time_diff <- as.factor(time_diffs[non_na_entries])
  
  #### fourth regressor: weekday of estimation date
  estimation_dates <- R_est_ts %>%
    transmute(across(!target_date,
                     list(td = ~ {as_date(substr(cur_column(), 13, 22))})))
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
  time_diff_regressors <- paste0("time_diff", (min_lag+1):max_lag)
  
  weekday_effects_target_date[,paste0(method, "_coef")] <- c(0, coef(model)[wd_td_regressors])
  weekday_effects_estimation_date[,paste0(method, "_coef")] <- c(0, coef(model)[wd_ed_regressors])
  time_diff_effects[as.character(min_lag:max_lag),paste0(method, "_coef")] <- c(0, coef(model)[time_diff_regressors])
  
  weekday_effects_target_date[wds[2:7], paste0(method, "_l")] <- confint(model)[wd_td_regressors, "2.5 %"]
  weekday_effects_target_date[wds[2:7], paste0(method, "_u")] <- confint(model)[wd_td_regressors, "97.5 %"]
  weekday_effects_estimation_date[wds[2:7], paste0(method, "_l")] <- confint(model)[wd_ed_regressors, "2.5 %"]
  weekday_effects_estimation_date[wds[2:7], paste0(method, "_u")] <- confint(model)[wd_ed_regressors, "97.5 %"]
  time_diff_effects[as.character((min_lag+1):max_lag),
                    paste0(method, "_l")] <- confint(model)[time_diff_regressors, "2.5 %"]
  time_diff_effects[as.character((min_lag+1):max_lag),
                    paste0(method, "_u")] <- confint(model)[time_diff_regressors, "97.5 %"]
}

write.csv(weekday_effects_target_date,
          paste0("Rt_estimate_reconstruction/otherFiles/weekday_effects_target_date_", number_days, ".csv"))
write.csv(weekday_effects_estimation_date,
          paste0("Rt_estimate_reconstruction/otherFiles/weekday_effects_estimation_date_", number_days, ".csv"))
write.csv(time_diff_effects,
          paste0("Rt_estimate_reconstruction/otherFiles/time_diff_effects_", number_days, ".csv"))

source("Rt_estimate_reconstruction/prepared_plots.R")

for (date_type in c("target", "estimation")){
  for (number_days in c(7,14,28)){
    data <- read.csv(paste0("Rt_estimate_reconstruction/otherFiles/weekday_effects_",
                            date_type, "_date_", number_days, ".csv"), row.names = 1)
    
    plot_data <- data %>%
      rownames_to_column("weekday") %>%
      gather(variable, value, -weekday) %>%
      separate(variable,
               into = c("method", "type"),
               sep = "_(?!.*_)",
               extra = "merge",
               remove = TRUE) %>%
      mutate(method = plyr::mapvalues(method,
                                      c("Braunschweig", "ETHZ_sliding_window", "globalrt_7d", "ilmenau", "RKI_7day"),
                                      c("HZI",          "ETH",                 "globalrt",    "Ilmenau", "RKI"))) %>%
      arrange(method) %>%
      pivot_wider(names_from = type, values_from = value) %>%
      mutate(exp_coef = exp(coef), exp_l = exp(l), exp_u = exp(u))
    
    methods_legend <- unique(plot_data$method)
    col_values <- get_colors(methods = methods_legend, palette = "methods")
    
    plot_raw <- ggplot(data=plot_data,
                       aes(x = factor(weekday, levels = wds),
                           group = method,
                           color = method)) +
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
      geom_line(aes(y = coef),
                size = .8, na.rm = T) +
      scale_x_discrete() + 
      scale_color_manual(values = col_values, name = "method") +
      ylab(paste0("coefficient for weekday of ",
                  date_type, " date (", number_days, " days in regression)"))
    
    print(plot)
    
    ylim <- c(0.9,1.15)
    
    plot <- plot_raw +
      geom_line(aes(y = exp_coef),
                size = .8, na.rm = T) +
      geom_text(data=subset(plot_data, exp_coef < ylim[1]),
                aes(label = round(exp_coef, 2)),
                y = ylim[1]) +
      geom_text(data=subset(plot_data, exp_coef > ylim[2]),
                aes(label = round(exp_coef, 2)),
                y = ylim[2]) +
      scale_color_manual(values = col_values, name = "method") +
      geom_ribbon(aes(ymin = exp_l, ymax = exp_u, fill = method), alpha = .15, colour = NA) +
      scale_fill_manual(values=col_values) +
      scale_x_discrete() + 
      ylab(paste0("exp(coefficient) for weekday of ",
                  date_type, " date (", number_days, " days in regression)")) +
      coord_cartesian(ylim=ylim) +
      geom_hline(yintercept = 1, size=.8)
    print(plot)
    
    ggsave(plot, filename = paste0("Figures/weekday_effects/influence_weekday_",
                                   date_type, "_date_", number_days, ".pdf"),
           bg = "transparent", width = 8, height = 5.8)
  }
}



