library(ggpubr)

setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
getwd()



################################
# get functions and parameters #
################################

source("Rt_estimate_reconstruction/load_data.R")

available_countries <- read.csv("Rt_estimate_reconstruction/otherFiles/available_countries.csv", row.names = 1)
pub_delays <- read.csv("Rt_estimate_reconstruction/otherFiles/pub_delays.csv", row.names = 1)


methods <- c("Braunschweig", "epiforecasts", "ETHZ_sliding_window",
             "globalrt_7d", "ilmenau", "RKI_7day", "rtlive", "SDSC")
method <- methods[2]
country = "DE"
conf_level = "95"
path_estimates = "reproductive_numbers/data-processed/"

#############################################
# plot coverage of 95% confidence intervals #
#############################################

start_default <- "2020-10-01"
start_globalrt <- "2021-02-18"
start_ilmenau <- "2020-11-19"

calc_consistence_metrics <- function(methods,
                                     country = "DE",
                                     conf_level = "95",
                                     path_estimates = "reproductive_numbers/data-processed/") {
  
  methods_CI <- methods[methods!="Braunschweig"]
  n_CI <- length(methods_CI)
  n <- length(methods)
  
  df_CI_coverage <- data.frame(matrix(rep(NA, n_CI*25), nrow = n_CI), row.names = methods_CI)
  colnames(df_CI_coverage) <- c(0:20, "num_CIs", "min_lag", "min_pub_date", "max_pub_date")
  
  df_CI_width <- data.frame(matrix(rep(NA, n_CI*25), nrow = n_CI), row.names = methods_CI)
  colnames(df_CI_width) <- c(0:20, "num_CIs", "min_lag", "min_pub_date", "max_pub_date")
  
  df_diff_to_first <- data.frame(matrix(rep(NA, n*21), nrow = n), row.names = methods)
  colnames(df_diff_to_first) <- c(0:20)
  
  df_diff_to_prev <- data.frame(matrix(rep(NA, n*21), nrow = n), row.names = methods)
  colnames(df_diff_to_prev) <- c(0:20)
  
  df_abs_diff_to_prev <- data.frame(matrix(rep(NA, n*21), nrow = n), row.names = methods)
  colnames(df_abs_diff_to_prev) <- c(0:20)
  
  df_diff_to_final <- data.frame(matrix(rep(NA, n*21), nrow = n), row.names = methods)
  colnames(df_diff_to_final) <- c(0:20)
  
  df_abs_diff_to_final <- data.frame(matrix(rep(NA, n*21), nrow = n), row.names = methods)
  colnames(df_abs_diff_to_final) <- c(0:20)
  
  df_labels <- data.frame(matrix(rep(NA, n*21), nrow = n), row.names = methods)
  colnames(df_labels) <- c(0:20)
  
  for (method in methods){
    print(method)
    
    max_lag <- 20
    days_until_final <- 70
    
    if (method == "globalrt_7d"){
      start_date <- start_globalrt
    } else if (method == "ilmenau"){
      start_date <- start_ilmenau
    } else {
      start_date <- start_default
    }

    end_date <- as.character(as_date(as_date(start_date) + weeks(10)) %m+% months(5))
    final_version <- as.character(as.Date(end_date) + days(days_until_final))
    if (method == "epiforecasts") final_version <- "2021-07-16" # "2021-07-19" is not available
    if (method == "Braunschweig") final_version <- "2021-07-20" # "2021-07-19" is not available
    
    print(paste0("end_date: ", end_date, ", final_version: ", final_version))

    pub_dates <- list.files(paste0(path_estimates, method),
                            pattern = "\\d{4}-\\d{2}-\\d{2}",
                            full.names = F) %>% substr(1, 10)
    pub_dates <- pub_dates[which(as_date(pub_dates) <= as_date(final_version) &
                                   as_date(pub_dates) >= as_date(start_date) - days(max_lag))]
    
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
                                                 include_label = TRUE,
                                                 verbose = F) %>%
            dplyr::select("date", "label", "R_pub") # %>%
            # rename(R_final = R_pub)
        },
        error = function(e) e
      )
      
      if (inherits(CI_not_available, "error")){
        print(paste(conf_level, "%-CI not available for method:", method))
        
      } else {
        for (pub_date in pub_dates[1:length(pub_dates)-1]){
          tryCatch(
            {
              cols <- c(date = NA, label = NA, R_pub = NA, lower = NA, upper = NA)
                
              R_est <- load_published_R_estimates(method,
                                                  start = min(as_date(pub_dates)) - max_lag,
                                                  end = as_date(pub_date) - min_lag,
                                                  pub_date = pub_date,
                                                  location = country,
                                                  conf_level = conf_level,
                                                  include_label = TRUE,
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
          available_pub_dates_long <- colnames(R_est_ts[4:dim(R_est_ts)[2]]) %>%
            substr(7, 16) %>%
            unique()
          available_pub_dates <- available_pub_dates_long[(as_date(available_pub_dates_long) >= start_date) &
                                                            (as_date(available_pub_dates_long) <= end_date)]
          
          ######
          # CI coverage rates and width...
          
          # if CIs are available calculate mean coverage rates and widths
          if (sum(!is.na(R_est_ts %>% dplyr::select(starts_with("lower")))) > 0){
            R_covered <- data.frame(date = R_est_ts$date)
            R_covered_difftime <-
              CI_width <-
              data.frame(estimated_after = make_difftime(day = seq(max_lag, min_lag),
                                                         units = "day"))
  
            for (pd in available_pub_dates){
              if ((as_date(pd) + days_until_final) %in% as_date(available_pub_dates)){
                i_final <- paste0("R_pub_", as_date(pd) + days_until_final)
                R_covered[pd] <- ifelse(is.na(R_est_ts[, paste0("lower_", pd)]), NA,
                                        ifelse((R_est_ts[,i_final] >= R_est_ts[, paste0("lower_", pd)]) &
                                                 (R_est_ts[,i_final] <= R_est_ts[, paste0("upper_", pd)]),
                                               TRUE,
                                               FALSE))
                
                if (dim(na.omit(R_covered[pd]))[1] >= (max_lag - min_lag + 1)) {
                  ea <- difftime(as.Date(pd), R_covered[,1], units = "day")
                  ea_wanted <- ea[ea <= max_lag & ea >= min_lag]
                  
                  R_covered_difftime[R_covered_difftime$estimated_after
                                     %in% ea_wanted, pd] <- R_covered[which(ea <= max_lag & ea >= min_lag), pd]
                  
                  CI_width[CI_width$estimated_after
                           %in% ea_wanted, pd] <- R_est_ts[R_est_ts$date %in%
                                                             R_covered[which(ea <= max_lag & ea >= min_lag), "date"],
                                                           paste0("width_", pd)]
                }
              }
            }
            
            R_covered_difftime <- R_covered_difftime %>%
              column_to_rownames(var = "estimated_after")
            CI_width <- CI_width %>%
              column_to_rownames(var = "estimated_after")
              
            df_CI_coverage[method, rownames(R_covered_difftime)] <- rowMeans(R_covered_difftime, na.rm = TRUE)
            df_CI_coverage[method, "min_lag"] <- min_lag
            df_CI_coverage[method, "num_CIs"] <- dim(R_covered_difftime)[2]
            df_CI_coverage[method, "min_pub_date"] <- min(colnames(R_covered_difftime))
            df_CI_coverage[method, "max_pub_date"] <- max(colnames(R_covered_difftime))
            
            df_CI_width[method, rownames(CI_width)] <- rowMeans(CI_width, na.rm = TRUE)
            df_CI_width[method, "min_lag"] <- min_lag
            df_CI_width[method, "num_CIs"] <- dim(CI_width)[2]
            df_CI_width[method, "min_pub_date"] <- min(colnames(CI_width))
            df_CI_width[method, "max_pub_date"] <- max(colnames(CI_width))
          }
          
          
          ######
          # MAD to first estimate...

          # for each target date extract date and value of first estimate
          {R_pub <- R_est_ts %>%
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
          first_est <- first_est %>%
            na.omit() %>%
            dplyr::filter(target_date <= end_date)
          
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
          
          idx <- as.character(intersect(0:max_lag, abs_diff_first_mean$lag))
          df_diff_to_first[method, idx] <- abs_diff_first_mean %>%
            dplyr::filter(as.character(lag) %in% idx) %>%
            arrange(lag) %>%
            dplyr::select(diff) %>% pull()}
          
          ######
          # M(A)D to previous estimate...
          
          {R_pub <- R_est_ts %>%
            dplyr::select(starts_with("R_pub_") | date) %>%
            column_to_rownames("date") %>%
            rename_with(~ substr(.x, 7, 16)) %>%
            dplyr::select(all_of(available_pub_dates))
          diff_prev <- data.frame(estimated_after = make_difftime(day = seq(max_lag, min_lag),
                                                                  units = "day")) %>%
            column_to_rownames(var = "estimated_after")
          
          for (c in colnames(R_pub)){
            if ((as_date(c) - 1) %in% as_date(colnames(R_pub))){
              j <- which(as_date(colnames(R_pub)) == as_date(c) - 1)
              diff <- R_pub[,c] - R_pub[,j]
              diff_prev[as.character(as_date(c) -
                                       as_date(rownames(R_pub[which(!is.na(diff)),]))),
                        c] <- diff %>% na.omit
            }
          }
          
          idx <- intersect(rownames(diff_prev), as.character(0:max_lag))
          df_diff_to_prev[method, idx] <- rowMeans(diff_prev, na.rm = TRUE)[idx]
          df_abs_diff_to_prev[method, idx] <- rowMeans(abs(diff_prev), na.rm = TRUE)[idx]}
          
          ######
          # M(A)D to final estimate (constant lag)...
          
          {R_pub <- R_est_ts %>%
            dplyr::select(starts_with("R_pub_") | date) %>%
            column_to_rownames("date") %>%
            rename_with(~ substr(.x, 7, 16)) %>%
            dplyr::select(all_of(available_pub_dates))
          diff_final <- data.frame(estimated_after = make_difftime(day = seq(0, max_lag),
                                                                   units = "day")) %>%
            column_to_rownames(var = "estimated_after")
          
          for (i_realtime in colnames(R_pub)[as_date(colnames(R_pub)) <= end_date]){
            if ((as_date(i_realtime) + days_until_final) %in% as_date(colnames(R_pub))){
              i_final <- as.character(as_date(i_realtime) + days_until_final)
              diff <- R_pub[,i_final] - R_pub[,i_realtime]
              i_diff <- as.character(as_date(i_realtime) - as_date(rownames(R_pub[which(!is.na(diff)),])))
              diff_final[i_diff, i_realtime] <- diff %>% na.omit()
            }
          }
          
          idx <- intersect(rownames(diff_final), as.character(0:max_lag))
          df_diff_to_final[method, idx] <- rowMeans(diff_final, na.rm = TRUE)[idx]
          df_abs_diff_to_final[method, idx] <- rowMeans(abs(diff_final), na.rm = TRUE)[idx]}
          
          ######
          # deprecated (nicht mehr verwendet und noch nicht angepasst an neue available pub dates)
          # M(A)D to final estimate (constant day of final estimate)... 
          {
          # diff_final <- R_est_ts %>%
          #   dplyr::select(date | R_final | starts_with("R_pub_")) %>%
          #   # calculate difference to "final" estimate
          #   mutate(across(!c(date,R_final), function(x) (R_final - x))) %>%
          #   dplyr::select(!R_final) %>%
          #   # add columns with difference between pub_date and target_date
          #   mutate(across(!date, list(lag = ~ {as_date(substr(cur_column(), 7, 16)) - date}))) %>%
          #   # order columns alphabetically
          #   dplyr::select(colnames(.)[order(colnames(.))]) %>%
          #   column_to_rownames("date") %>%          
          #   # select columns corresponding to pub_dates after start_date
          #   dplyr::select(colnames(.)[as_date(substr(colnames(.), 7, 16)) >= start_date])
          # 
          # # reshape wide to long
          # cols <- colnames(dplyr::select(diff_final, !ends_with("_lag")))
          # diff_final_long <- data.frame()
          # for (col in cols){
          #   df <- diff_final[,c(col, paste0(col, "_lag"))] %>% setNames(c("diff", "lag"))
          #   diff_final_long <- bind_rows(df, diff_final_long)
          # }
          # diff_final_long <- na.omit(diff_final_long) %>% group_by(lag)
          # diff_final_mean <- diff_final_long %>%
          #   summarise(diff = mean(diff))
          # abs_diff_final_mean <- diff_final_long %>%
          #   summarise(diff = mean(abs(diff)))
          # 
          # idx <- as.character(intersect(0:max_lag, diff_final_long$lag))
          # df_diff_to_final[method, idx] <- diff_final_mean %>%
          #   dplyr::filter(as.character(lag) %in% idx) %>%
          #   arrange(lag) %>%
          #   dplyr::select(diff) %>% pull()
          # df_abs_diff_to_final[method, idx] <- abs_diff_final_mean %>%
          #   dplyr::filter(as.character(lag) %in% idx) %>%
          #   arrange(lag) %>%
          #   dplyr::select(diff) %>% pull()
          }
          
          ######
          # extract labels...
          
          {labels <- R_est_ts %>%
            dplyr::select(date | starts_with("label_")) %>%
            # add columns with difference between pub_date and target_date
            mutate(across(!date, list(lag = ~ {as_date(substr(cur_column(), 7, 16)) - date}))) %>%
            # order columns alphabetically
            dplyr::select(colnames(.)[order(colnames(.))]) %>%
            column_to_rownames("date") %>%          
            # select columns corresponding to pub_dates after start_date
            dplyr::select(colnames(.)[(as_date(substr(colnames(.), 7, 16)) >= start_date) &
                                        (as_date(substr(colnames(.), 7, 16)) <= end_date)])
          
          # reshape wide to long
          cols <- colnames(dplyr::select(labels, !ends_with("_lag")))
          labels_long <- data.frame()
          for (col in cols){
            df <- labels[,c(col, paste0(col, "_lag"))] %>% setNames(c("label", "lag"))
            labels_long <- bind_rows(df, labels_long)
          }
          labels_long <- na.omit(labels_long) %>% group_by(lag)
          
          getmode <- function(x) {
            ux <- unique(x)
            tab <- tabulate(match(x, ux))
            return(ux[tab == max(tab)])
          }
          
          labels_mode <- labels_long %>%
            summarise(label = getmode(label))

          idx <- as.character(intersect(0:max_lag, labels_mode$lag))
          df_labels[method, idx] <- labels_mode %>%
            dplyr::filter(as.character(lag) %in% idx) %>%
            arrange(lag) %>%
            dplyr::select(label) %>% pull()}
        }
      }
    } else {
      print(paste("No estimates from", method, "for", country))
    }
  }
  
  return(list(df_CI_coverage, df_CI_width,
              df_abs_diff_to_prev, df_abs_diff_to_final,
              df_diff_to_prev, df_diff_to_final,
              df_diff_to_first, df_labels))
}

methods <- c("Braunschweig", "epiforecasts", "ETHZ_sliding_window",
             "globalrt_7d", "ilmenau", "RKI_7day", "rtlive", "SDSC")

CI_eval <- calc_consistence_metrics(methods)
View(CI_eval[[2]])
{write.csv(CI_eval[[1]], "Rt_estimate_reconstruction/otherFiles/95_CI_coverage.csv")
write.csv(CI_eval[[2]], "Rt_estimate_reconstruction/otherFiles/95_CI_width.csv")
write.csv(CI_eval[[3]], "Rt_estimate_reconstruction/otherFiles/abs_diff_to_prev.csv")
write.csv(CI_eval[[4]], "Rt_estimate_reconstruction/otherFiles/abs_diff_to_final.csv")
write.csv(CI_eval[[5]], "Rt_estimate_reconstruction/otherFiles/diff_to_prev.csv")
write.csv(CI_eval[[6]], "Rt_estimate_reconstruction/otherFiles/diff_to_final.csv")
write.csv(CI_eval[[7]], "Rt_estimate_reconstruction/otherFiles/diff_to_first.csv")
write.csv(CI_eval[[8]], "Rt_estimate_reconstruction/otherFiles/estimate_labels.csv")}

source("Rt_estimate_reconstruction/prepared_plots.R")

plot1 <- plot_CI_coverage_rates(); plot1
plot2 <- plot_CI_widths(); plot2
plot3 <- plot_diff_final(diff_type = "abs_diff"); plot3
plot4 <- plot_diff_final(diff_type = "diff", ylim = c(-0.02, 0.037)); plot4

consistence_plot <- ggarrange(plot1, plot2, plot3, plot4, ncol=2, nrow=2,
                              labels = list("A", "B", "C", "D"),
                              font.label = list(size = 18, face = "bold"),
                              common.legend = T, legend="bottom",
                              legend.grob = get_legend(plot3))
print(consistence_plot)
ggsave(consistence_plot, filename = paste0("Figures/CI/consistence_plots.pdf"),
       bg = "transparent", width = 16, height = 11.6)

# appendix
plot_diff_prev_A <- plot_diff_prev(diff_type = "abs_diff", ylim = c(-0.01, 0.115)); plot_diff_prev_A
plot_diff_prev_B <- plot_diff_prev(diff_type = "diff", ylim = c(-0.006, 0.021)); plot_diff_prev_B
plot_diff_prev <- ggarrange(plot_diff_prev_A, plot_diff_prev_B, ncol=2, nrow=1,
                            labels = list("A", "B"),
                            font.label = list(size = 18, face = "bold"),
                            common.legend = T, legend="bottom",
                            legend.grob = get_legend(plot_diff_prev_A))
print(plot_diff_prev)
ggsave(plot_diff_prev, filename = paste0("Figures/CI/diff_prev_plots.pdf"),
       bg = "transparent", width = 16, height = 5)

plot_abs_diff_first(max_lag=20)


