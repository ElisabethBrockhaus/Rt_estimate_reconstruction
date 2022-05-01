library(data.table)

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
methods <- list.dirs(path_estimates, full.names = F)
methods <- methods[!methods %in% c("", "AW_7day", "AW_WVday", "owid", "ETHZ_step",
                                   "ETHZ_sliding_window_deaths", "ETHZ_step_deaths")]
methods

available_countries <- read.csv("Rt_estimate_reconstruction/otherFiles/available_countries.csv", row.names = 1)
pub_delays <- read.csv("Rt_estimate_reconstruction/otherFiles/pub_delays.csv", row.names = 1)


#################################
# plot estimates as time series #
#################################

#start <- "2021-04-01"
#end <- "2021-05-01"

start <- "2021-02-16"
end <- "2021-03-18"

for (method in methods){
  print(method)
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)

  final_version <- "2021-07-16"
  if (method == "Braunschweig") final_version <- "2021-07-18"
  if (method == "epiforecasts") final_version <- as.character(as.Date(start) + 106)
  
  pub_dates <- pub_dates[which(pub_dates <= end &
                                 pub_dates >= start)]
  end_date <- as_date(max(pub_dates))
  pub_dates <- c(pub_dates, final_version)
  for (country in c("DE", "AT", "CH")[1]){
    print(country)
    if (available_countries[method, country]) {
      if (exists("R_est_ts")) rm(R_est_ts)
      if (exists("R_est")) rm(R_est)
      for (pub_date in pub_dates){
        tryCatch(
          {
            R_est <- load_published_R_estimates(method,
                                                end = end_date,
                                                pub_date = pub_date,
                                                location = country,
                                                verbose = F)
            if (pub_date != final_version){
              last <- max(R_est[rowSums(!is.na(R_est))>1, "date"])
              R_est <- R_est %>% dplyr::filter(date <= last, date > last - 7)
            }
            names(R_est) <- c("date", paste0("R.", pub_date),
                              paste0("lower.", pub_date), paste0("upper.", pub_date))
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
        if (!file.exists(folder)) {
          dir.create(folder)
        }
        
        plot_real_time_estimates(R_est_ts,
                                 start_date = last_date - 30,
                                 end_date = last_date,
                                 plot_title = paste(method, country),
                                 name_consensus = final_version,
                                 filenames = paste0(folder_ending,
                                                    method, "_", country, ".png"))
      }
    } else {
      print(paste("No estimates from", method, "for", country, "."))
    }
  }
}



##############################################################
# plot estimates over days between estimation and target day #
##############################################################

start_date <- as_date("2020-11-16")
end_date <- as_date("2021-07-16")

# plot monthly
target_dates <- seq(start_date, end_date, by = "month")

# plot 4-weekly
#target_dates <- seq(start_date, end_date, by = "week")[1+4*(0:11)]

min_lag <- 0
max_lag <- 30

for (method in methods){
  print(method)
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which((as_date(pub_dates) <= end_date) &
                                 (as_date(pub_dates) >= start_date))]
  end_date <- as_date(max(pub_dates))
  for (country in c("DE", "AT", "CH")[1]){
    print(country)
    if (available_countries[method, country]) {
      if (exists("R_est_ts")) rm(R_est_ts)
      if (exists("R_est")) rm(R_est)
      for (pub_date in pub_dates){
        tryCatch(
          {
            R_est <- load_published_R_estimates(method,
                                                start = as_date(pub_date) - max_lag,
                                                end = as_date(pub_date) - min_lag,
                                                pub_date = pub_date,
                                                location = country,
                                                verbose = F) %>%
              dplyr::select("date", "R_pub") %>%
              mutate(estimated_after = as_date(pub_date) - date) %>%
              reshape(idvar = "estimated_after", timevar = "date", direction = "wide") %>%
              rename_with(~ gsub(".", "_", .x, fixed = TRUE)) %>%
              dplyr::select("estimated_after",
                            which((colnames(.)[2+(0:(max_lag-min_lag))] %>%
                                     substr(7, 16) %>%
                                     as_date()) %in% target_dates) + 1)
          },
          error = function(e) {R_est <<- data.frame(estimated_after = make_difftime(day = seq(min_lag, max_lag), units = "day"))}
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
        
        R_est_plot <- R_est_ts %>%
          dplyr::select("estimated_after", which(colMeans(is.na(.)) < 0.67))       
        
        R_est_plot <- R_est_plot %>% rename(date = estimated_after)
        
        tryCatch(
          {
            resulted_target_dates <- intersect(target_dates,
                                               colnames(R_est_plot[2:dim(R_est_plot)[2]]) %>%
                                                 substr(7, 16) %>%
                                                 as_date()) %>% as_date()
            min_lag_plot <- Inf
            for (col in 2:dim(R_est_plot)[2]) {
              min_lag_plot <- min(min_lag_plot, which(!is.na(R_est_plot[,col]))[1] - 1)
            }
            max_lag_plot <- min(min_lag_plot + 6, max_lag)
            
            plot_for_comparison(R_est_plot,
                                comp_methods = resulted_target_dates,
                                start_date = make_difftime(day = min_lag_plot),
                                end_date = make_difftime(day = max_lag_plot),
                                ylim_l=0.5, ylim_u=1.5,
                                legend_name = "estimated R for",
                                plot_title = paste(method, country),
                                col_palette = "Spectral",
                                filenames = paste0("_realtime_corrections_over_time/", method, "_", country, ".png"),
                                verbose = F)
          },
          error = function(c) {cat("Too many estimates missing. \n Minimal proportion of NA values:",
                                   min(colMeans(is.na(R_est_ts[,2:dim(R_est_ts)[2]]))), "\n")}
        )
      }
    } else {
      print(paste("No estimates from", method, "for", country, "."))
    }
  }
}






#########
# same plot with averages for weekdays
#########

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
                                                start = as_date(pub_date) - max_lag,
                                                end = as_date(pub_date) - min_lag,
                                                pub_date = pub_date,
                                                location = country,
                                                verbose = F) %>%
              dplyr::select("date", "R_pub") %>%
              mutate(estimated_after = as_date(pub_date) - date) %>%
              reshape(idvar = "estimated_after", timevar = "date", direction = "wide") %>%
              rename_with(~ gsub(".", "_", .x, fixed = TRUE)) %>%
              dplyr::select("estimated_after",
                            which((colnames(.)[2+(0:(max_lag-min_lag))] %>%
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
        
        #min_lag_plot <- Inf
        #for (col in 2:dim(R_est_ts)[2]) {
        #  min_lag_plot <- min(min_lag_plot,
        #                      as.numeric(R_est_ts[which(!is.na(R_est_ts[,col]))[1], "estimated_after"]),
        #                      na.rm = TRUE)
        #}
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
                                               method, "_", country, ".png"),
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
                                                start = as_date(pub_date) - max_lag,
                                                end = as_date(pub_date) - min_lag,
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
                             filenames = paste0(method, "_", country, ".png"),
                             verbose = F)
      }
    } else {
      print(paste("No estimates from", method, "for", country, "."))
    }
  }
}



