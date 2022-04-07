setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()



################################
# get functions and parameters #
################################

source("Rt_estimate_reconstruction/load_data.R")

# sources of published real-time estimates
methods <- list.dirs(path_estimates, full.names = F)
methods <- methods[!methods %in% c("", "AW_7day", "AW_WVday", "owid", "ETHZ_step",
                                   "ETHZ_sliding_window_deaths", "ETHZ_step_deaths",
                                   "zidatalab")]
methods

available_countries <- read.csv("Rt_estimate_reconstruction/otherFiles/available_countries.csv", row.names = 1)
pub_delays <- read.csv("Rt_estimate_reconstruction/otherFiles/pub_delays_mode.csv", row.names = 1)


#############################################
# plot coverage of 95% confidence intervals #
#############################################

calc_CI_coverages <- function(methods,
                              country = "DE",
                              start_date = as_date("2020-11-16"),
                              end_date = as_date("2021-05-01"),
                              conf_level = "95",
                              path_estimates = "reproductive_numbers/data-processed/",
                              save_csv = F) {
  
  n <- length(methods)
  CI_coverage <- data.frame(matrix(rep(NA, n*24), nrow = n), row.names = methods)
  colnames(CI_coverage) <- c(0:20, "num_CIs", "min_pub_date", "max_pub_date")

  for (method in methods){
    print(method)
    pub_dates <- list.files(paste0(path_estimates, method),
                            full.names = F) %>% substr(1, 10)
    pub_dates <- pub_dates[which(as_date(pub_dates) <= end_date &
                                   as_date(pub_dates) >= start_date)]
    
    final_version <- "2021-07-16"
    if (method == "Braunschweig") final_version <- "2021-07-18"
    
    if (available_countries[method, country]) {
      if (exists("R_est_ts")) rm(R_est_ts)
      if (exists("R_est")) rm(R_est)
      
      min_lag <- pub_delays[method, country]
      max_lag <- min_lag + 6
      
      CI_not_available <- tryCatch(
        {
          R_est_ts <- load_published_R_estimates(method,
                                                 start = start_date - max_lag,
                                                 end = as_date(final_version),
                                                 pub_date = final_version,
                                                 location = country,
                                                 conf_level = conf_level,
                                                 verbose = F) %>%
            dplyr::select("date", "R_pub")
          names(R_est_ts)[2] <- "R_final"
        },
        error = function(e) e
      )
      
      if (inherits(CI_not_available, "error")){
        print(paste(conf_level, "%-CI not available for method:", method))
        
      } else {
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
                dplyr::select("date", "lower", "upper")
              
              names(R_est) <- c("date",
                                paste0("02.5q_", pub_date),
                                paste0("97.5q_", pub_date))
            },
            error = function(e) {R_est <<- data.frame(date = seq(as_date(pub_date) - max_lag,
                                                                 as_date(pub_date),
                                                                 by = "day"))}
          )
          R_est_ts <- R_est_ts %>% full_join(R_est, by="date")
        }
        if (exists("R_est_ts")){
          available_pub_dates <- colnames(R_est_ts[3:dim(R_est_ts)[2]]) %>%
            substr(7, 16) %>%
            unique()
          
          R_covered <- data.frame(date = R_est_ts$date)
          R_covered_difftime <- data.frame(estimated_after = make_difftime(day = seq(min_lag, max_lag),
                                                                           units = "day"))
          
          for (pd in available_pub_dates){
            R_covered[pd] <- ifelse(is.na(R_est_ts[, paste0("02.5q_", pd)]), NA,
                                    ifelse((R_est_ts$R_final >= R_est_ts[, paste0("02.5q_", pd)]) &
                                             (R_est_ts$R_final <= R_est_ts[, paste0("97.5q_", pd)]),
                                           TRUE,
                                           FALSE))
            if (dim(na.omit(R_covered[pd]))[1] == 7) {
              R_covered_difftime[pd] <- na.omit(R_covered[pd])
            }
          }
          
          ea <- R_covered_difftime$estimated_after
          CI_coverage[method, as.character(ea)] <- rowMeans(R_covered_difftime[,-1])
          CI_coverage[method, "num_CIs"] <- dim(R_covered_difftime)[2]
          CI_coverage[method, "min_pub_date"] <- min(colnames(R_covered_difftime[,-1]))
          CI_coverage[method, "max_pub_date"] <- max(colnames(R_covered_difftime[,-1]))
    
        }
      }
    } else {
      print(paste("No estimates from", method, "for", country))
    }
  }
  
  View(CI_coverage)
  return(CI_coverage)
}

CI_coverage_95 <- calc_CI_coverages(methods)
write.csv(CI_coverage_95, "Rt_estimate_reconstruction/otherFiles/95_CI_coverage.csv")

CI_coverage_50 <- calc_CI_coverages(methods[c(2,3,4,8,11,12,13,14,15)], conf_level = 50)
write.csv(CI_coverage_50, "Rt_estimate_reconstruction/otherFiles/50_CI_coverage.csv")

source("Rt_estimate_reconstruction/prepared_plots.R")

plot_CI_coverage_rates()
plot_CI_coverage_rates("50")
