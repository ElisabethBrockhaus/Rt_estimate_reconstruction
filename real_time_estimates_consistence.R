setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
getwd()



################################
# get functions and parameters #
################################

source("Rt_estimate_reconstruction/load_data.R")

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
  
  methods_CI <- methods[methods!="Braunschweig"]
  n_CI <- length(methods_CI)
  n <- length(methods)
  CI_coverage <- data.frame(matrix(rep(NA, n_CI*25), nrow = n_CI), row.names = methods_CI)
  colnames(CI_coverage) <- c(0:20, "num_CIs", "min_lag", "min_pub_date", "max_pub_date")
  CI_width_mean <- data.frame(matrix(rep(NA, n_CI*25), nrow = n_CI), row.names = methods_CI)
  colnames(CI_width_mean) <- c(0:20, "num_CIs", "min_lag", "min_pub_date", "max_pub_date")
  diff_first_mean <- data.frame(matrix(rep(NA, n*21), nrow = n), row.names = methods)
  colnames(diff_first_mean) <- c(0:20)
  diff_prev_mean <- data.frame(matrix(rep(NA, n*21), nrow = n), row.names = methods)
  colnames(diff_prev_mean) <- c(0:20)
  
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
              cols <- c(date = NA, R_pub = NA, lower = NA, upper = NA)
                
              R_est <- load_published_R_estimates(method,
                                                  start = as_date(pub_date) - max_lag,
                                                  end = as_date(pub_date) - min_lag,
                                                  pub_date = pub_date,
                                                  location = country,
                                                  conf_level = conf_level,
                                                  verbose = F) %>%
                dplyr::select(any_of(names(cols)))
              
              R_est <- R_est %>%
                add_column(!!!cols[setdiff(names(cols), names(R_est))]) %>%
                mutate(width = upper - lower) %>%
                rename_with(~ paste0(.x, "_", pub_date), !date)
              
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
          
          
          ######
          # CI coverage rates and width...
          
          # if CIs are available calculate mean coverage rates and widths
          if (dim(R_est_ts %>% dplyr::select(starts_with("lower")) %>% na.omit)[1] > 0){
            R_covered <- data.frame(date = R_est_ts$date)
            R_covered_difftime <-
              CI_width <-
              data.frame(estimated_after = make_difftime(day = seq(max_lag, min_lag),
                                                         units = "day"))
  
            for (pd in available_pub_dates){
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
              }
            }
            
            R_covered_difftime <- R_covered_difftime %>%
              column_to_rownames(var = "estimated_after")
            CI_width <- CI_width %>%
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
          }
          
          
          ######
          # MAD to first estimate...

          # for each target date extract date and value of first estimate
          R_pub <- R_est_ts %>%
            dplyr::select(starts_with("R_pub_") | date) %>%
            column_to_rownames("date")
          first_est <- R_est_ts %>%
            dplyr::select(date) %>%
            rename(target_date = date)
          first_est$first_date <- apply(R_pub, 1,
                                        function(x) names(which(!is.na(x)))[1]) %>%
            substr(7, 16)
          first_est$first_R <- apply(R_pub, 1,
                                     function(x) x[which(!is.na(x))][1])
          first_est <- first_est %>% na.omit()
          
          abs_diff_first <- R_est_ts %>%
            dplyr::select(date | starts_with("R_pub_")) %>%
            # only use rows for which the first published R value is in "first_est"
            dplyr::filter(date %in% first_est$target_date) %>%
            # calculate absolute difference to first published estimate for the target_date
            mutate(across(!date, function(x) abs(x - first_est$first_R))) %>%
            # add columns with difference between pub_dates
            mutate(across(!date, list(lag = ~ {as_date(substr(cur_column(), 7, 16)) -
                as_date(first_est[first_est$target_date == date, "first_date"])} ))) %>%
            # order columns alphabetically
            dplyr::select(colnames(.)[order(colnames(.))]) %>%
            column_to_rownames("date") %>%
            # select columns corresponding to pub_dates after start_date
            dplyr::select(colnames(.)[as_date(substr(colnames(.), 7, 16)) >= start_date])
          
          # reshape wide to long
          cols <- colnames(dplyr::select(abs_diff_first, !ends_with("_lag")))
          abs_diff_first_long <- data.frame()
          for (col in cols){
            df <- abs_diff_first[,c(col, paste0(col, "_lag"))] %>% setNames(c("diff", "lag"))
            abs_diff_first_long <- bind_rows(df, abs_diff_first_long)
          }
          abs_diff_first_long <- na.omit(abs_diff_first_long)
          abs_diff_first_mean <- abs_diff_first_long %>%
            group_by(lag) %>%
            summarise(diff = mean(diff))
          
          diff_first_mean[method, as.character(abs_diff_first_mean$lag)] <- abs_diff_first_mean$diff
      
          
          ######
          # MAD to previous estimate...
          
          R_pub <- R_pub %>%
            rename_with(~ substr(.x, 7, 16))
          abs_diff_prev <- data.frame(estimated_after = make_difftime(day = seq(max_lag, min_lag),
                                                                      units = "day")) %>%
            column_to_rownames(var = "estimated_after")
          
          for (c in colnames(R_pub)){
            if ((as_date(c) - 1) %in% as_date(colnames(R_pub))){
              j <- which(as_date(colnames(R_pub)) == as_date(c) - 1)
              abs_diff <- abs(R_pub[,c] - R_pub[,j])
              abs_diff_prev[as.character(as_date(c) -
                              as_date(rownames(R_pub[which(!is.na(abs_diff)),]))),
                            c] <- abs_diff %>% na.omit
            }
          }
          
          diff_prev_mean[method, rownames(abs_diff_prev)] <- rowMeans(abs_diff_prev, na.rm = TRUE)
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

