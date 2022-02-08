setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
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
methods <- methods[!methods %in% c("", "AW_7day", "AW_WVday",
                                   "ETHZ_sliding_window_deaths", "ETHZ_step_deaths")]
methods



#################################
# plot estimates as time series #
#################################

for (method in methods){
  print(method)
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which(pub_dates <= "2021-12-31")]
  pub_dates <- pub_dates[length(pub_dates)-27:0]
  end_date <- as_date(max(pub_dates))
  for (country in c("DE", "AT", "CH")){
    if (exists("R_est_ts")) rm(R_est_ts)
    if (exists("R_est")) rm(R_est)
    print(country)
    tryCatch(
      {
        for (pub_date in pub_dates){
          R_est <- load_published_R_estimates(method,
                                              end = end_date,
                                              pub_date = pub_date,
                                              location = country,
                                              verbose = F)
          if (pub_date != end_date){
            last <- max(R_est[rowSums(!is.na(R_est))>1, "date"])
            R_est <- R_est %>% dplyr::filter(date <= last, date > last - 7)
          }
          names(R_est) <- c("date", paste0("R_pub_", pub_date),
                            paste0("lower_", pub_date), paste0("upper_", pub_date))
          if (!exists("R_est_ts")){
            R_est_ts <- R_est
          } else{
            R_est_ts <- R_est_ts %>% full_join(R_est, by="date")
          }
        }
      },
      error = function(c) {print(paste("No estimates from", method,
                                       "for", country, "."))}
    )
    if (exists("R_est_ts")){
      last_date <- max(R_est_ts[rowSums(!is.na(R_est_ts))>1, "date"])
      plot_for_comparison(R_est_ts,
                          comp_methods = pub_dates,
                          start_date = last_date - 30,
                          end_date = last_date,
                          name_consensus = end_date,
                          legend_name = "published on",
                          plot_title = paste(method, country),
                          col_palette = "Spectral",
                          filenames = paste0("_realtime_", method, "_", country, ".png"),
                          verbose = F)
    }
  }
}



##############################################################
# plot estimates over days between estimation and target day #
##############################################################
source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

for (method in methods){
  print(method)
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which(pub_dates <= "2021-12-31")]
  pub_dates <- pub_dates[length(pub_dates)-27:0]
  end_date <- as_date(max(pub_dates))
  for (country in c("DE", "AT", "CH")[1]){
    print(country)
    if (exists("R_est_ts")) rm(R_est_ts)
    if (exists("R_est")) rm(R_est)
    tryCatch(
      {
        for (pub_date in pub_dates){
          R_est <- load_published_R_estimates(method,
                                              start = as_date(pub_date) - 40,
                                              end = as_date(pub_date),
                                              pub_date = pub_date,
                                              location = country,
                                              verbose = F) %>%
            dplyr::select("date", "R_pub") %>%
            mutate(estimated_after = as_date(pub_date) - date) %>%
            reshape(idvar = "estimated_after", timevar = "date", direction = "wide") %>%
            rename_with(~ gsub(".", "_", .x, fixed = TRUE))
            
          if (!exists("R_est_ts")){
            R_est_ts <- R_est
          } else{
            R_est_ts <- R_est_ts %>% merge(R_est, all = T) %>%
              group_by(estimated_after) %>%
              summarise_all(list(~ mean(., na.rm = T)))
            R_est_ts <- R_est_ts[order(colnames(R_est_ts))]
          }
        }
      },
      error = function(c) {print(paste("No estimates from", method,
                                       "for", country, "."))}
    )
    if (exists("R_est_ts")){
      R_est_ts <- R_est_ts[, c(TRUE, rep(FALSE, 19),
                               (seq_len(ncol(R_est_ts)) %% 7 == 2)[21:(ncol(R_est_ts)-14)],
                               rep(FALSE, 14))]
      target_dates <- colnames(R_est_ts)[2:dim(R_est_ts)[2]] %>% substr(7, 16)
      R_est_ts <- R_est_ts %>% rename(date = estimated_after)
      plot_for_comparison(R_est_ts,
                          comp_methods = target_dates,
                          start_date = make_difftime(day = 0),
                          end_date = make_difftime(day = 30),
                          legend_name = "estimated R for",
                          plot_title = paste(method, country),
                          col_palette = "Spectral",
                          filenames = paste0("_realtime_", method, "_", country, "_corrections_over_time.png"),
                          verbose = F)
    }
  }
}
source("Rt_estimate_reconstruction/prepared_plots.R")

# TODO: Solve get_colors in prepared_plots.R line 43 ff.
