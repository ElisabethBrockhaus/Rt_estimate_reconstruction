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
methods <- methods[!methods %in% c("", "AW_7day", "AW_WVday", "owid", "ETHZ_step",
                                   "ETHZ_sliding_window_deaths", "ETHZ_step_deaths")]
methods

available_countries <- read.csv("Rt_estimate_reconstruction/otherFiles/available_countries.csv", row.names = 1)



#############################################
# plot coverage of 95% confidence intervals #
#############################################

for (method in methods[3]){
  print(method)
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which(pub_dates <= "2021-05-01" &
                                 pub_dates >= "2021-04-01")]
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
                                                start = as_date(pub_date) - 7,
                                                end = as_date(pub_date),
                                                pub_date = pub_date,
                                                location = country,
                                                verbose = F)
            names(R_est) <- c("date",
                              paste0("50.0q_", pub_date),
                              paste0("02.5q_", pub_date),
                              paste0("97.5q_", pub_date))
          },
          error = function(e) {R_est <<- data.frame(date = seq(as_date("2019-12-28"),
                                                               as_date(end_date),
                                                               by = "day"))}
        )
        if (!exists("R_est_ts")){
          R_est_ts <- R_est
        } else{
          R_est_ts <- R_est_ts %>% full_join(R_est, by="date")
        }
      }
      if (exists("R_est_ts")){
        available_median_dates <- colnames(R_est_ts[2:dim(R_est_ts)[2]] %>%
                                             dplyr::select(starts_with("50.0q"))) %>%
          substr(7, 16) %>%
          as_date()
        
        
        
        

      }
    } else {
      print(paste("No estimates from", method, "for", country, "."))
    }
  }
}