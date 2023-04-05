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
ggsave(plots_arranged, filename = "Figures/estimates_real_time_all_methods.pdf",
       bg = "transparent", width = 17, height = 21)


