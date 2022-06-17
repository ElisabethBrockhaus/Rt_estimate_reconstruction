setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
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
pub_delays <- read.csv("Rt_estimate_reconstruction/otherFiles/pub_delays.csv", row.names = 1)


#############################################
# plot coverage of 95% confidence intervals #
#############################################

calc_CI_coverages <- function(methods,
                              country = "DE",
                              start_date_min = as_date("2020-11-16"),
                              end_date = as_date("2021-05-01"),
                              conf_level = "95",
                              path_estimates = "reproductive_numbers/data-processed/") {
  
  n <- length(methods)
  CI_coverage <- data.frame(matrix(rep(NA, n*25), nrow = n), row.names = methods)
  colnames(CI_coverage) <- c(0:20, "num_CIs", "min_lag", "min_pub_date", "max_pub_date")
  CI_width_mean <- data.frame(matrix(rep(NA, n*25), nrow = n), row.names = methods)
  colnames(CI_width_mean) <- c(0:20, "num_CIs", "min_lag", "min_pub_date", "max_pub_date")
  diff_first_mean <- data.frame(matrix(rep(NA, n*25), nrow = n), row.names = methods)
  colnames(diff_first_mean) <- c(0:20, "num_datapoints", "min_lag", "min_pub_date", "max_pub_date")
  diff_prev_mean <- data.frame(matrix(rep(NA, n*25), nrow = n), row.names = methods)
  colnames(diff_prev_mean) <- c(0:20, "num_datapoints", "min_lag", "min_pub_date", "max_pub_date")
  
  for (method in methods){
    print(method)
    
    max_lag <- 20
    
    if (method == "epiforecasts") {
      start_date <- max(start_date_min, as_date("2021-04-14"))
    } else if (method == "globalrt_7d") {
      start_date <- max(start_date_min, as_date("2021-02-15"))
    } else {
      start_date <- start_date_min
    }
    
    pub_dates <- list.files(paste0(path_estimates, method),
                            full.names = F) %>% substr(1, 10)
    pub_dates <- pub_dates[which(as_date(pub_dates) <= end_date &
                                   as_date(pub_dates) >= start_date - max_lag)]
    
    final_version <- "2021-07-16"
    if (method == "Braunschweig") final_version <- "2021-07-20"
    
    if (available_countries[method, country]) {
      if (exists("R_est_ts")) {rm(R_est_ts)}
      if (exists("R_est")) {rm(R_est)}
      
      min_lag <- pub_delays[method, country]
      
      CI_not_available <- tryCatch(
        {
          R_est_ts <- load_published_R_estimates(method,
                                                 start = min(as_date(pub_dates)) - max_lag,
                                                 end = as_date(final_version),
                                                 pub_date = final_version,
                                                 location = country,
                                                 conf_level = conf_level,
                                                 verbose = F) %>%
            dplyr::select("date", "R_pub") %>%
            rename(R_final = R_pub)
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
                dplyr::select("date", "R_pub", "lower", "upper") %>%
                mutate(width = upper - lower)
              
              names(R_est) <- c("date",
                                paste0("R_pub_", pub_date),
                                paste0("lower_", pub_date),
                                paste0("upper_", pub_date),
                                paste0("width_", pub_date))
              
            },
            error = function(e) {R_est <<- data.frame(date = seq(as_date(pub_date) - max_lag,
                                                                 as_date(pub_date) - min_lag,
                                                                 by = "day"))}
          )
          R_est_ts <- R_est_ts %>% full_join(R_est, by="date")
        }
        if (exists("R_est_ts")){
          available_pub_dates_long <- colnames(R_est_ts[3:dim(R_est_ts)[2]]) %>%
            substr(7, 16) %>%
            unique()
          available_pub_dates <- available_pub_dates_long[as_date(available_pub_dates_long) >= start_date]

          R_covered <- data.frame(date = R_est_ts$date)
          R_covered_difftime <-
            CI_width <-
            abs_diff_first <-
            abs_diff_prev <-
            data.frame(estimated_after = make_difftime(day = seq(max_lag, min_lag),
                                                       units = "day"))
          
          for (i in seq_along(available_pub_dates)){
            pd <- available_pub_dates[i]
            R_covered[pd] <- ifelse(is.na(R_est_ts[, paste0("lower_", pd)]), NA,
                                    ifelse((R_est_ts$R_final >= R_est_ts[, paste0("lower_", pd)]) &
                                             (R_est_ts$R_final <= R_est_ts[, paste0("upper_", pd)]),
                                           TRUE,
                                           FALSE))
            
            if (dim(na.omit(R_covered[pd]))[1] == (max_lag - min_lag + 1)) {
              
              ea <- difftime(as.Date(pd), R_covered[,1], units = "day")
              indices <- which(!is.na(R_covered[pd]))
              
              R_covered_difftime[R_covered_difftime$estimated_after
                                 %in% ea[indices], pd] <- na.omit(R_covered[pd])
              
              CI_width[CI_width$estimated_after
                       %in% ea[indices], pd] <- na.omit(R_est_ts[, paste0("width_", pd)])
              
              j <- as_date(pd) - min(as_date(available_pub_dates)) + 1
              
              for (k in 0:(max_lag-min_lag)) {
                first <- as_date(pd) - k
                
                if (first %in% as_date(available_pub_dates_long)) {
                  abs_diff_first[max_lag-min_lag+1-k, pd] <-
                  abs(R_est_ts[j+max_lag-min_lag-k, paste0("R_pub_", pd)] -
                        R_est_ts[j+max_lag-min_lag-k, paste0("R_pub_", first)])
                }
              }
              
              prev <- as_date(pd) - 1
              if (prev %in% as_date(available_pub_dates_long)) {
                abs_diff_prev[abs_diff_prev$estimated_after
                              %in% ea[indices], pd] <-
                  abs(R_est_ts[j+0:(max_lag-min_lag), paste0("R_pub_", pd)] -
                        R_est_ts[j+0:(max_lag-min_lag), paste0("R_pub_", prev)])
              }
            }
          }
          
          R_covered_difftime <- R_covered_difftime %>%
            column_to_rownames(var = "estimated_after")
          CI_width <- CI_width %>%
            column_to_rownames(var = "estimated_after")
          abs_diff_first <- abs_diff_first %>%
            column_to_rownames(var = "estimated_after")
          abs_diff_prev <- abs_diff_prev %>%
            column_to_rownames(var = "estimated_after")
            
          CI_coverage[method, rownames(R_covered_difftime)] <- rowMeans(R_covered_difftime)
          CI_coverage[method, "min_lag"] <- min_lag
          CI_coverage[method, "num_CIs"] <- dim(R_covered_difftime)[2]
          CI_coverage[method, "min_pub_date"] <- min(colnames(R_covered_difftime))
          CI_coverage[method, "max_pub_date"] <- max(colnames(R_covered_difftime))
          
          CI_width_mean[method, rownames(CI_width)] <- rowMeans(CI_width)
          CI_width_mean[method, "min_lag"] <- min_lag
          CI_width_mean[method, "num_CIs"] <- dim(CI_width)[2]
          CI_width_mean[method, "min_pub_date"] <- min(colnames(CI_width))
          CI_width_mean[method, "max_pub_date"] <- max(colnames(CI_width))
          
          diff_first_mean[method, rownames(abs_diff_first)] <- rowMeans(abs_diff_first, na.rm = TRUE)
          diff_first_mean[method, "min_lag"] <- min_lag
          diff_first_mean[method, "num_datapoints"] <- dim(abs_diff_first)[2]
          diff_first_mean[method, "min_pub_date"] <- min(colnames(abs_diff_first))
          diff_first_mean[method, "max_pub_date"] <- max(colnames(abs_diff_first))
          
          diff_prev_mean[method, rownames(abs_diff_prev)] <- rowMeans(abs_diff_prev, na.rm = TRUE)
          diff_prev_mean[method, "min_lag"] <- min_lag
          diff_prev_mean[method, "num_datapoints"] <- dim(abs_diff_prev)[2]
          diff_prev_mean[method, "min_pub_date"] <- min(colnames(abs_diff_prev))
          diff_prev_mean[method, "max_pub_date"] <- max(colnames(abs_diff_prev))
        }
      }
    } else {
      print(paste("No estimates from", method, "for", country))
    }
  }
  
  return(list(CI_coverage, CI_width_mean, diff_first_mean, diff_prev_mean))
}

CI_eval <- calc_CI_coverages(c("Braunschweig", "epiforecasts", "ETHZ_sliding_window",
                               "globalrt_7d", "ilmenau", "RKI_7day", "rtlive", "SDSC"))
write.csv(CI_eval[[1]], "Rt_estimate_reconstruction/otherFiles/95_CI_coverage.csv")
write.csv(CI_eval[[2]], "Rt_estimate_reconstruction/otherFiles/95_CI_width.csv")
write.csv(CI_eval[[3]], "Rt_estimate_reconstruction/otherFiles/diff_to_first.csv")
write.csv(CI_eval[[4]], "Rt_estimate_reconstruction/otherFiles/diff_to_prev.csv")

source("Rt_estimate_reconstruction/prepared_plots.R")

plot_CI_coverage_rates()
plot_CI_widths()
plot_abs_diff_first()
plot_abs_diff_prev()

CI_coverage_50 <- calc_CI_coverages(c("Braunschweig", "epiforecasts", "rtlive"),
                                    conf_level = 50)
write.csv(CI_coverage_50[[1]], "Rt_estimate_reconstruction/otherFiles/50_CI_coverage.csv")
write.csv(CI_coverage_50[[2]], "Rt_estimate_reconstruction/otherFiles/50_CI_width.csv")

plot_CI_coverage_rates("50")
plot_CI_widths("50")

