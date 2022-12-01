setwd("../..")
getwd()

source("Rt_estimate_reconstruction/load_data.R")

# path of published estimates
path_estimates <- "reproductive_numbers/data-processed/"
available_countries <- read.csv("Rt_estimate_reconstruction/otherFiles/available_countries.csv", row.names = 1)

# sources of published real-time estimates
methods <- c("Braunschweig", "ETHZ_sliding_window", "RKI_7day",
             "rtlive", "SDSC", "epiforecasts", "ilmenau", "globalrt_7d")

n <- length(methods)
init <- rep(NA, n)

min_lag <- 0
max_lag <- 20

pub_delays_list <- list()
for (method in methods) {
  pub_delays_list[[method]] <- data.frame(DE = rep(0, max_lag+1),
                                          AT = rep(0, max_lag+1),
                                          CH = rep(0, max_lag+1),
                                          row.names = min_lag:max_lag)
}

start_default <- "2020-10-01"
start_globalrt <- "2021-02-18"
start_ilmenau <- "2020-11-19"

for (method in methods){
  print(method)
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  
  if (method == "globalrt_7d"){
    start <- start_globalrt
  } else if (method == "ilmenau"){
    start <- start_ilmenau
  } else {
    start <- start_default
  }
  end <- as.character(as_date(as_date(start) + weeks(10)) %m+% months(5))
  
  pub_dates <- pub_dates[which((pub_dates <= end) & (pub_dates >= start))]
  print(length(pub_dates))

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
          error = function(e) {if (country == "DE") print(paste("Error for:", method, pub_date))}
        )
      }
    }
  }
}

pub_delays_list
saveRDS(pub_delays_list, file="Rt_estimate_reconstruction/otherFiles/pub_delays_list.RData")

# summarize pub_delays: one value per method and country
min_pub_delays_list <- readRDS("Rt_estimate_reconstruction/otherFiles/pub_delays_list.RData")
min_pub_delays_list

delays_given_min_share <- data.frame(total_num = init,
                                     d0.8 = init, d0.85 = init, d0.9 = init, d0.95 = init,
                                     row.names = methods)
country <- "DE"

for (method in methods){
  num.est <- sum(min_pub_delays_list[[method]][, country])
  delays_given_min_share[method, "total_num"] <- num.est
  cum.num.est <- cumsum(min_pub_delays_list[[method]][, country])
  
  for (share in seq(0.8,0.95,0.05)){
    min.index <- which.max(cum.num.est/num.est >= share)
    delays_given_min_share[method, paste0("d", share)] <- row.names(min_pub_delays_list[[method]])[min.index]
  }
}
delays_given_min_share # d0.x is the min_pub_delay given at least x% of the estimates are available

pub_delays <- data.frame(DE = init, AT = init, CH = init, row.names = methods)

for (method in methods){
  for (country in c("DE", "AT", "CH")){
    num.est <- sum(min_pub_delays_list[[method]][, country])
    print(paste("total number of estimates", num.est))
    if(num.est > 0){
      cum.num.est <- cumsum(min_pub_delays_list[[method]][, country])
      min.index <- which.max(cum.num.est/num.est >= 0.9)
      pub_delays[method, country] <- row.names(min_pub_delays_list[[method]])[min.index]
    }
  }
}

pub_delays
write.csv(pub_delays, "Rt_estimate_reconstruction/otherFiles/pub_delays.csv")


