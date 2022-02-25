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
init <- rep(NA, 16)

available_countries <- data.frame(DE = init, AT = init, CH = init, row.names = methods)

for (method in methods){
  print(method)
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which((pub_dates <= "2021-12-31") & (pub_dates >= "2021-01-01"))]
  end_date <- as_date(max(pub_dates))
  for (country in c("DE", "AT", "CH")){
    print(country)
    for (pub_date in pub_dates) {
      tryCatch(
        {
          temp <- load_published_R_estimates(method,
                                             start = as_date(pub_date) - 30,
                                             end = as_date(pub_date),
                                             pub_date = pub_date,
                                             location = country,
                                             verbose = F)
          available_countries[method, country] <- TRUE
          break
        },
        error = function(e) {
          if (is.na(available_countries[method, country])) {
            available_countries[method, country] <<- FALSE
          }
        }
      )
    }
  }
}

available_countries
write.csv(available_countries, "Rt_estimate_reconstruction/otherFiles/available_countries.csv")
          
          