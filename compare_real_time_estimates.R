library(data.table)
library(ggpubr)

setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers"
getwd()



################################
# get functions and parameters #
################################

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

# path of published estimates
path_estimates <- "reproductive_numbers/data-processed/"

# sources of published real-time estimates
methods <- c("Braunschweig", "ETHZ_sliding_window", "RKI_7day",
             "rtlive", "SDSC", "epiforecasts", "ilmenau", "globalrt_7d")

available_countries <- read.csv("Rt_estimate_reconstruction/otherFiles/available_countries.csv", row.names = 1)
pub_delays <- read.csv("Rt_estimate_reconstruction/otherFiles/pub_delays.csv", row.names = 1)

# period for which RKI was criticized to correct always upwards
start_default <- "2020-10-01"
start_globalrt <- "2021-02-15"
start_ilmenau <- "2020-11-16"


#################################
# plot estimates as time series #
#################################

plots_CI <- list()

for (method in methods){
  print(method)
  
  if (method == "globalrt_7d"){
    start <- start_globalrt
  } else if (method == "ilmenau"){
    start <- start_ilmenau
  } else {
    start <- start_default
  }
  end <- as.character(as_date(start) + weeks(10))
  
  pub_dates <- list.files(paste0(path_estimates, method),
                          pattern = "\\d{4}-\\d{2}-\\d{2}",
                          full.names = F) %>% substr(1, 10)
  
  final_version <- as.character(as_date(end) %m+% months(6))
  if (method == "epiforecasts") final_version <- as.character(as.Date(start) + 105)
  if (method == "globalrt_7d") final_version <- "2021-10-28"
  if (method == "Braunschweig") final_version <- "2021-06-09"
  
  pub_dates_available <- pub_dates[which(pub_dates <= end &
                                 pub_dates >= start)]
  start_date <- as_date(start)
  end_date <- as_date(end)
  
  pub_dates_wanted <- seq(from=start_date+weeks(1), to=end_date, by="week")
  
  pub_dates <- c()
  for (pd in pub_dates_wanted){
    pub_dates <- c(pub_dates,
                   pub_dates_available[which.min(abs(as_date(pub_dates_available) - as_date(pd)))])
  }
  
  pub_dates <- c(unique(pub_dates), final_version)
  
  for (country in c("DE", "AT", "CH")[1]){
    print(country)
    if (available_countries[method, country]) {
      if (exists("R_est_ts")) rm(R_est_ts)
      if (exists("R_est")) rm(R_est)
      for (pub_date in pub_dates){
        tryCatch(
          {
            R_est <- load_published_R_estimates(method,
                                                end = as_date(pub_date) - pub_delays[method, country],
                                                pub_date = pub_date,
                                                location = country,
                                                include_label = TRUE,
                                                verbose = F) %>%
              dplyr::select(any_of(c("date", "label", "R_pub", "lower", "upper")))
            if ((pub_date != final_version) & (dim(R_est)[2] > 1)){
              last <- max(R_est[rowSums(!is.na(R_est))>1, "date"])
              R_est <- R_est %>% dplyr::filter(date <= last, date > last - 21)
            }
            names <- c("date",
                       paste0("label.", pub_date),
                       paste0("R.", pub_date),
                       paste0("l.", pub_date),
                       paste0("u.", pub_date))
            names(R_est) <- names[1:dim(R_est)[2]]
          },
          error = function(e) {R_est <<- data.frame(date = seq(as_date("2019-12-28"),
                                                               as_date(final_version),
                                                               by = "day"))}
        )
        if (!exists("R_est_ts")){
          R_est_ts <- R_est
        } else{
          R_est_ts <- R_est_ts %>% full_join(R_est, by="date")
        }
      }
      if (exists("R_est_ts")){
        last_date <- max(R_est_ts[rowSums(!is.na(R_est_ts))>1, "date"])
        
        folder_ending <- paste0("_realtime_raw_CI/", end, "/")
        folder <- paste0("Figures/estimates", folder_ending)
        if (!dir.exists(folder)) {
          dir.create(folder)
        }
        
        if (method == "ilmenau"){
          ylim <- c(0.15, 1.75)
        } else if(method == "Braunschweig"){
          ylim <- c(0.9, 1.8)
        } else {
          ylim <- c(0.9, 1.5)
        }
        
        title <- ifelse(method == "Braunschweig", "HZI",
                        ifelse(method == "ETHZ_sliding_window", "ETH",
                               ifelse(method == "globalrt_7d", "globalrt",
                                      ifelse(method == "ilmenau", "Ilmenau",
                                             ifelse(method == "RKI_7day", "RKI",
                                                    method)))))
        
        plots_CI[[method]] <- plot_real_time_estimates_with_CI(R_est_ts,
                                                               start_date = start,
                                                               end_date = end,
                                                               plot_title = title,
                                                               name_consensus = final_version,
                                                               filenames = paste0(folder_ending,
                                                                                  method, "_", country, ".pdf"),
                                                               ylim_l = ylim[1], ylim_u = ylim[2])
      }
    } else {
      print(paste("No estimates from", method, "for", country, "."))
    }
  }
}

plots_arranged <- ggarrange(plots_CI[[methods[1]]], plots_CI[[methods[2]]], plots_CI[[methods[3]]],
                            plots_CI[[methods[4]]], plots_CI[[methods[5]]], plots_CI[[methods[6]]],
                            plots_CI[[methods[7]]], plots_CI[[methods[8]]],
                            ncol=1, nrow=length(plots_CI),
                            #common.legend = T, legend="bottom",
                            #legend.grob = get_legend(plots_CI[["RKI_7day"]]))
                            legend = "none")
ggsave(plots_arranged, filename = paste0("Figures/estimates_realtime_raw_CI/",
                                         as.character(as_date(start_default) + weeks(10)),
                                         "/all_methods.pdf"),
       bg = "transparent", width = 17, height = 21)


###########
# plot epiforecast history in more detail

method <- "epiforecasts"
country <- "DE"
start <- start_default
end <- as.character(as.Date(start) + 105)
final_version <- end

pub_dates <- list.files(paste0(path_estimates, method),
                        pattern = "\\d{4}-\\d{2}-\\d{2}",
                        full.names = F) %>% substr(1, 10)

pub_dates <- pub_dates[which(pub_dates <= end &
                               pub_dates >= start)]
start_date <- as_date(min(pub_dates))
end_date <- as_date(max(pub_dates))

if (exists("R_est_ts")) rm(R_est_ts)
if (exists("R_est")) rm(R_est)
for (pub_date in pub_dates){
  tryCatch(
    {
      R_est <- load_published_R_estimates(method,
                                          end = as_date(pub_date) - pub_delays[method, country],
                                          pub_date = pub_date,
                                          location = country,
                                          verbose = F) %>%
        dplyr::select(any_of(c("date", "R_pub")))
      if ((pub_date != final_version) & (dim(R_est)[2] > 1)){
        last <- max(R_est[rowSums(!is.na(R_est))>1, "date"])
        R_est <- R_est %>% dplyr::filter(date <= last)
      }
      names(R_est) <- c("date", paste0("R.", pub_date))
    },
    error = function(e) {R_est <<- data.frame(date = seq(as_date("2019-12-28"),
                                                         as_date(final_version),
                                                         by = "day"))}
  )
  if (!exists("R_est_ts")){
    R_est_ts <- R_est
  } else{
    R_est_ts <- R_est_ts %>% full_join(R_est, by="date")
  }
}
if (exists("R_est_ts")){
  last_date <- max(R_est_ts[rowSums(!is.na(R_est_ts))>1, "date"])
  
  folder_ending <- paste0("_realtime_raw/", end, "/")
  folder <- paste0("Figures/estimates", folder_ending)
  if (!dir.exists(folder)) {
    dir.create(folder)
  }

  ylim <- c(0.9, 1.5)
  
  plot_real_time_estimates(R_est_ts,
                           start_date = start,
                           end_date = end,
                           plot_title = method,
                           name_consensus = final_version,
                           filenames = paste0(folder_ending,
                                              method, "_", country, ".pdf"),
                           ylim_l = ylim[1], ylim_u = ylim[2],
                           show.legend = FALSE)
}





#################################################
# plot different estimates for same target date #
#################################################

plots_lags <- list()

for (method in methods){
  print(method)
  
  if (method == "globalrt_7d"){
    start <- start_globalrt
  } else if (method == "ilmenau"){
    start <- start_ilmenau
  } else {
    start <- start_default
  }
  start_date <- as_date(start)
  end <- as.character(start_date + weeks(10))
  end_date <- as_date(end)
  
  pub_dates_available <- list.files(paste0(path_estimates, method),
                              pattern = "\\d{4}-\\d{2}-\\d{2}",
                              full.names = F) %>% substr(1, 10)
  
  pub_dates <- seq(start_date, end_date, by = "day")
  
  for (country in c("DE", "AT", "CH")[1]){
    if (available_countries[method, country]) {
      if (exists("R_est_ts")) rm(R_est_ts)
      for (pub_date in pub_dates){
        if (exists("R_est")) rm(R_est)
        for (lag in c(1,14,50,70)){
          if (exists("R_est_lag")) rm(R_est_lag)
          if ((as.character(as_date(pub_date) + days(lag)) %in% pub_dates_available) &
              (pub_delays[method, country] <= lag)){
              tryCatch(
                {
                  R_est_lag <- load_published_R_estimates(method,
                                                          start = as_date(pub_date),
                                                          end = as_date(pub_date),
                                                          pub_date = as.character(as_date(pub_date) + days(lag)),
                                                          location = country,
                                                          include_label = TRUE,
                                                          verbose = F) %>%
                    dplyr::select(any_of(c("date", "label", "R_pub")))
                  names(R_est_lag) <- c("date",
                                    paste0("label+", lag),
                                    paste0("R+", lag))
                },
                error = function(e) {R_est_lag <<- data.frame(date = seq(start_date, end_date, by = "day"))}
              )
              if (!exists("R_est")){
                R_est <- R_est_lag
              } else{
                R_est <- R_est %>% full_join(R_est_lag, by="date")
              }
          }
        }
        if(exists("R_est")){
          if (!exists("R_est_ts")){
            R_est_ts <- R_est
          } else{
            R_est_ts <- bind_rows(R_est_ts, R_est)
          }
        }
      }
      if (exists("R_est_ts")){
        folder_ending <- "_realtime_compare_lags/"
        folder <- paste0("Figures/estimates", folder_ending)
        if (!dir.exists(folder)) {
          dir.create(folder)
        }
        
        if (method == "ilmenau"){
          ylim <- c(0.2, 1.75)
        } else if(method == "Braunschweig"){
          ylim <- c(0.9, 1.8)
        } else {
          ylim <- c(0.9, 1.5)
        }
        
        title <- ifelse(method == "Braunschweig", "HZI",
                        ifelse(method == "ETHZ_sliding_window", "ETH",
                               ifelse(method == "globalrt_7d", "globalrt",
                                      ifelse(method == "ilmenau", "Ilmenau",
                                             ifelse(method == "RKI_7day", "RKI",
                                                    method)))))
        
        plots_lags[[method]] <- plot_estimates_with_different_lags(R_est_ts,
                                                                 start_date = start,
                                                                 end_date = end,
                                                                 plot_title = title,
                                                                 name_consensus = "70",
                                                                 filenames = paste0(folder_ending,
                                                                                    method, "_", country, ".pdf"),
                                                                 ylim_l = ylim[1], ylim_u = ylim[2])
      }
    } else {
      print(paste("No estimates from", method, "for", country, "."))
    }
  }
}

plots_arranged <- ggarrange(plots_lags[[methods[1]]], plots_lags[[methods[2]]], plots_lags[[methods[3]]],
                            plots_lags[[methods[4]]], plots_lags[[methods[5]]], plots_lags[[methods[6]]],
                            plots_lags[[methods[7]]], plots_lags[[methods[8]]],
                            ncol=1, nrow=length(plots_lags),
                            common.legend = T, legend="bottom",
                            legend.grob = get_legend(plots_lags[["SDSC"]]))
ggsave(plots_arranged, filename = "Figures/estimates_realtime_compare_lags/all_methods.pdf",
       bg = "transparent", width = 17, height = 21)


#####################################################################################
# plot estimates over days between estimation and target day averaged over weekdays #
#####################################################################################

# dates over which the mean is calculated
start_date <- as_date("2020-11-16")
end_date <- as_date("2021-05-01")
target_dates <- seq(start_date, end_date, by = "day")
wds <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")

min_lag <- 1
max_lag <- 30

for (method in methods){
  print(method)
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which((as_date(pub_dates) <= end_date) &
                                 (as_date(pub_dates) >= start_date))]
  for (country in c("DE", "AT", "CH")[1]){
    print(country)
    if (available_countries[method, country]) {
      if (exists("R_est_ts")) rm(R_est_ts)
      if (exists("R_est")) rm(R_est)
      for (pub_date in pub_dates){
        tryCatch(
          {
            R_est <- load_published_R_estimates(method,
                                                start = max(as_date(min(pub_dates)) - pub_delays[method, country],
                                                            as_date(pub_date) - max_lag),
                                                end = min(as_date(max(pub_dates)) - pub_delays[method, country] - 6,
                                                          as_date(pub_date) - min_lag),
                                                pub_date = pub_date,
                                                location = country,
                                                verbose = F) %>%
              dplyr::select("date", "R_pub") %>%
              mutate(estimated_after = as_date(pub_date) - date) %>%
              reshape(idvar = "estimated_after", timevar = "date", direction = "wide") %>%
              rename_with(~ gsub(".", "_", .x, fixed = TRUE)) %>%
              dplyr::select("estimated_after",
                            which((colnames(.)[-1] %>%
                                     substr(7, 16) %>%
                                     as_date()) %in% target_dates) + 1)
          },
          error = function(e) {
            R_est <<- data.frame(estimated_after = make_difftime(day = seq(min_lag, max_lag),
                                                                 units = "day"))
            }
        )
        
        if (!exists("R_est_ts")){
          R_est_ts <- R_est
        } else {
          R_est_ts <- R_est_ts %>% merge(R_est, all = T) %>%
            group_by(estimated_after) %>%
            summarise_all(list(~ mean(., na.rm = T)))
          R_est_ts <- R_est_ts[order(colnames(R_est_ts))]
        }
      }
      if (exists("R_est_ts")){
        
        min_lag_plot <- pub_delays[method, country]
        max_lag_plot <- min(min_lag_plot + 6, max_lag)
        
        R_est_ts <- R_est_ts %>%
          dplyr::filter((estimated_after >= min_lag_plot) &
                          (estimated_after <= max_lag_plot))
        
        # extract weekdays from target dates
        org_names <- colnames(R_est_ts)[2:(dim(R_est_ts)[2])]
        colnames(R_est_ts)[2:(dim(R_est_ts)[2])] <- weekdays(org_names %>%
                                                                 substr(7, 16) %>%
                                                                 as_date()) %>%
          paste0(org_names)
        
        # calculate mean over weekdays
        R_est_means <- R_est_ts["estimated_after"]
        for (wd in wds) {
          R_est_means[wd] <- R_est_ts %>%
            dplyr::select(starts_with(wd)) %>%
            apply(1, function(x) exp(mean(log(x), na.rm = TRUE)))
        }

        R_est_plot <- R_est_means %>% rename(date = estimated_after)
        
        plot_for_comparison(R_est_plot,
                            comp_methods = wds,
                            start_date = make_difftime(day = min_lag_plot),
                            end_date = make_difftime(day = max_lag_plot),
                            ylim_l=0.5, ylim_u=1.5,
                            legend_name = "mean estimate for",
                            plot_title = paste(method, country),
                            col_palette = "Spectral",
                            filenames = paste0("_realtime_corrections_per_weekday/",
                                               method, "_", country, ".pdf"),
                            verbose = F)
      }
    } else {
      print(paste("No estimates from", method, "for", country, "."))
    }
  }
}



############################
# weekday effects averaged #
############################

# pub_dates
start_date <- as_date("2020-11-16")
end_date <- as_date("2021-05-01")
wds <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")

for (method in methods){
  print(method)
  
  if (method == "Braunschweig") {
    print("Method handled manually due to irregular pub dates.")
    next
  }
  
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which((as_date(pub_dates) <= end_date) &
                                 (as_date(pub_dates) >= start_date))]
    
  for (country in c("DE", "AT", "CH")[1]){
    print(country)
    
    min_lag <- pub_delays[method, country]
    max_lag <- min_lag + 6
    
    if (available_countries[method, country]) {
      if (exists("R_est_list")) rm(R_est_list)
      if (exists("R_est")) rm(R_est)
      
      R_est_empty <- data.frame(target_weekday = wds,
                                Monday = rep(NA, 7),
                                Tuesday = rep(NA, 7),
                                Wednesday = rep(NA, 7),
                                Thursday = rep(NA, 7),
                                Friday = rep(NA, 7),
                                Saturday = rep(NA, 7),
                                Sunday = rep(NA, 7))
      
      for (pub_date in pub_dates){
        tryCatch(
          {
            R_est <- load_published_R_estimates(method,
                                                start = max(as_date(min(pub_dates)) - min_lag,
                                                            as_date(pub_date) - max_lag),
                                                end = min(as_date(max(pub_dates)) - max_lag,
                                                          as_date(pub_date) - min_lag),
                                                pub_date = pub_date,
                                                location = country,
                                                verbose = F) %>%
              mutate(target_weekday = weekdays(date)) %>%
              dplyr::select("target_weekday", "R_pub")
            R_est <- R_est[match(wds, R_est$target_weekday),]
          },
          error = function(e) {R_est <<- data.frame(target_weekday = wds,
                                                    R_pub = rep(NA, 7))}
        )
        
        weekday <- weekdays(as_date(pub_date))
        R_est <- R_est_empty %>% mutate(!!weekday := R_est$R_pub)
          
        if (!exists("R_est_list")){
          R_est_list <- R_est
        } else {
          R_est_list <- rbindlist(list(R_est_list, R_est))
        }
      }
      if (exists("R_est_list")){
        
        R_est_mean <- R_est_list[, lapply(.SD, function(x) exp(mean(log(x), na.rm = TRUE))),
                                 by = .(target_weekday)]
        
        plot_weekday_effects(R_est_mean,
                             ylim_l=0.75, ylim_u=1.25,
                             plot_title = paste0("Mean over latest 7 estimates per pub date (",
                                                 method, " ", country, ")"),
                             min_lag = min_lag, max_lag = max_lag,
                             filenames = paste0(method, "_", country, ".pdf"),
                             verbose = F)
      }
    } else {
      print(paste("No estimates from", method, "for", country, "."))
    }
  }
}

# Braunschweig: mean over pub weekdays, otherwise estimates not published often enough
{
  method <- "Braunschweig"
  country <- "DE"
  print(method)
  print(country)
    
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which((as_date(pub_dates) <= end_date) &
                                 (as_date(pub_dates) >= start_date))]
  
  
    
  min_lag <- pub_delays[method, country]
  max_lag <- min_lag + 6
  
  if (exists("R_est_list")) rm(R_est_list)
  if (exists("R_est")) rm(R_est)
  
  R_est_empty <- data.frame(target_weekday = wds,
                            R_pub = rep(NA, 7))
  
  for (pub_date in pub_dates){
    tryCatch(
      {
        R_est <- load_published_R_estimates(method,
                                            start = max(as_date(min(pub_dates)) - min_lag,
                                                        as_date(pub_date) - max_lag),
                                            end = min(as_date(max(pub_dates)) - max_lag,
                                                      as_date(pub_date) - min_lag),
                                            pub_date = pub_date,
                                            location = country,
                                            verbose = F) %>%
          mutate(target_weekday = weekdays(date)) %>%
          dplyr::select("target_weekday", "R_pub")
        R_est <- R_est[match(wds, R_est$target_weekday),]
      },
      error = function(e) {R_est <<- data.frame(target_weekday = wds,
                                                R_pub = rep(NA, 7))}
    )
    
    if (!exists("R_est_list")){
      R_est_list <- R_est
    } else {
      R_est_list <- rbindlist(list(R_est_list, R_est))
    }
  }
  if (exists("R_est_list")){
    
    R_est_list <- R_est_list %>% na.omit()
    
    R_est_mean <- R_est_list[, lapply(.SD, function(x) exp(mean(log(x), na.rm = TRUE))),
                             by = .(target_weekday)]
    
    plot_title <- paste0("Mean over latest 7 estimates per pub date (",
                         method, " ", country, ")")
    filenames <- paste0(method, "_", country, ".pdf")
    ylim_l <- 0.75
    ylim_u <- 1.25

    R_est <- R_est_mean %>%
      rename(mean_R = R_pub)
    
    wds <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
    
    # plot
    R_plot <- ggplot(data = R_est, aes(x = factor(target_weekday,
                                                  levels = wds),
                                       y = mean_R,
                                       group = 1)) +
      geom_line() +
      geom_hline(aes(yintercept = 1))
    
    R_plot <- R_plot +
      labs(x = "weekday (target date)", y = "Rt estimate",
           subtitle = paste(min_lag, "-", max_lag, "days previous to pub date")) +
      scale_x_discrete()
    
    R_plot <- R_plot +
      theme_minimal() +
      theme(
        plot.title = element_text(size=18),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        axis.line = element_line(),
        axis.line.y.right = element_line(),
        axis.line.x.top = element_line(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent")
      ) +
      ggtitle(plot_title)
    
    R_plot <- R_plot +
      coord_cartesian(ylim = c(ylim_l, ylim_u), expand = FALSE)
    
    ggsave(R_plot, filename = paste0("Figures/estimates_realtime_target_date_influence/", filenames),  bg = "transparent",
           width = 13.1, height = 5.8)
    print(R_plot)
  }
}
