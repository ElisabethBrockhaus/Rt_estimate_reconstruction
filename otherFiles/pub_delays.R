setwd("../..")
getwd()

source("Rt_estimate_reconstruction/load_data.R")

# path of published estimates
path_estimates <- "reproductive_numbers/data-processed/"

# sources of published real-time estimates
methods <- list.dirs(path_estimates, full.names = F)
methods <- methods[!methods %in% c("", "AW_7day", "AW_WVday", "owid", "ETHZ_step",
                                   "ETHZ_sliding_window_deaths", "ETHZ_step_deaths")]
methods

n <- length(methods)
init <- rep(NA, n)
pub_delays <- data.frame(DE = init, AT = init, CH = init, row.names = methods)

for (method in methods){
  print(method)
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which((pub_dates <= "2021-07-18") & (pub_dates >= "2020-11-16"))]

  for (country in c("DE", "AT", "CH")){
    print(country)
    min_lag <- 0
    
    if (available_countries[method, country]) {
      if (exists("R_est")) rm(R_est)
      for (pub_date in pub_dates){
        tryCatch(
          {
            R_est <- load_published_R_estimates(method,
                                                start = as_date(pub_date) - 30,
                                                end = as_date(pub_date) - min_lag,
                                                pub_date = pub_date,
                                                location = country,
                                                verbose = F) %>%
              dplyr::select("date", "R_pub") %>%
              mutate(estimated_after = as_date(pub_date) - date) %>%
              dplyr::select("estimated_after", "R_pub") %>%
              arrange(estimated_after)
            
            min_lag <- as.numeric(R_est[which(!is.na(R_est[,"R_pub"]))[1], "estimated_after"])

          },
          error = function(e) {}
        )
      }
      pub_delays[method, country] <- min_lag
    }
  }
}

pub_delays
write.csv(pub_delays, "Rt_estimate_reconstruction/otherFiles/pub_delays.csv")



pub_delays_list <- list()
for (method in methods) {
  pub_delays_list[[method]] <- data.frame(DE = rep(0, 17),
                                     AT = rep(0, 17),
                                     CH = rep(0, 17),
                                     row.names = 0:16)
}

min_lag <- 0
max_lag <- 16

for (method in methods){
  print(method)
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which((pub_dates <= "2021-07-18") & (pub_dates >= "2020-11-16"))]

  for (country in c("DE", "AT", "CH")){
    print(country)
    
    if (available_countries[method, country]) {
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
              dplyr::select("estimated_after", "R_pub") %>%
              arrange(estimated_after)
            
            min_delay <- as.numeric(R_est[which(!is.na(R_est[,"R_pub"]))[1], "estimated_after"])
            
            pub_delays_list[[method]][as.character(min_delay), country] <- 1 + pub_delays_list[[method]][as.character(min_delay), country]
            
          },
          error = function(e) {}
        )
      }
      pub_delays[method, country] <- min_lag
    }
  }
}

pub_delays_list
saveRDS(pub_delays_list, file="Rt_estimate_reconstruction/otherFiles/pub_delays_list.RData")
min_pub_delays_list <- readRDS("Rt_estimate_reconstruction/otherFiles/pub_delays_list.RData")
          