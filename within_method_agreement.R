library(ggpubr)

setwd("/home/johannes/Documents/Projects/R_t")
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
method <- methods[5]
country = "DE"
conf_level = "95"
days_until_final = 70
path_estimates = "reproductive_numbers/data-processed/"

#############################################
# plot coverage of 95% confidence intervals #
#############################################

start_default <- "2020-10-01"
start_globalrt <- "2021-02-15"
start_ilmenau <- "2020-11-16"

calc_consistence_metrics <- function(methods,
                                     country = "DE",
                                     conf_level = "95",
                                     days_until_final = 70,
                                     path_estimates = "reproductive_numbers/data-processed/") {
  
  methods_CI <- methods[methods!="Braunschweig"]
  n_CI <- length(methods_CI)
  n <- length(methods)
  
  df_CI_coverage <- data.frame(matrix(rep(NA, (n_CI+2)*44), nrow = n_CI+2), row.names = c(methods_CI, "ETH_old", "ETH_new"))
  colnames(df_CI_coverage) <- c(0:39, "num_CIs", "min_lag", "min_pub_date", "max_pub_date")
  
  df_CI_width <- data.frame(matrix(rep(NA, (n_CI+2)*44), nrow = n_CI+2), row.names = c(methods_CI, "ETH_old", "ETH_new"))
  colnames(df_CI_width) <- c(0:39, "num_CIs", "min_lag", "min_pub_date", "max_pub_date")
  
  df_diff_to_final <- data.frame(matrix(rep(NA, (n+2)*41), nrow = n+2), row.names = c(methods, "ETH_old", "ETH_new"))
  colnames(df_diff_to_final) <- c(0:39, "num_est")
  
  df_abs_diff_to_final <- data.frame(matrix(rep(NA, (n+2)*41), nrow = n+2), row.names = c(methods, "ETH_old", "ETH_new"))
  colnames(df_abs_diff_to_final) <- c(0:39, "num_est")
  
  df_switch_final <- df_switch2_final <- data.frame(matrix(rep(NA, (n+2)*41), nrow = n+2), row.names = c(methods, "ETH_old", "ETH_new"))
  colnames(df_switch_final) <- colnames(df_switch2_final) <- c(0:39, "num_est")
  
  df_labels <- data.frame(matrix(rep(NA, n*41), nrow = n), row.names = methods)
  colnames(df_labels) <- c(0:39, "num_est")
  
  for (method in methods){
    print(method)
    
    max_lag <- 39

    if (method == "globalrt_7d"){
      start_date <- start_globalrt
    } else if (method == "ilmenau"){
      start_date <- start_ilmenau
    } else {
      start_date <- start_default
    }

    end_date <- "2021-09-30"
    if (method == "rtlive") end_date <- as.character(as_date("2021-07-31") - days(days_until_final))
    
    pub_dates <- list.files(paste0(path_estimates, method),
                            pattern = "\\d{4}-\\d{2}-\\d{2}",
                            full.names = F) %>% substr(1, 10)
    
    time_diffs <- as_date(pub_dates) - as_date(end_date)
    time_diffs[time_diffs<days_until_final] <- 1000 # "big-M" in order to only consider dates at least days_until_final days after end_date 
    final_version <- pub_dates[which.min(time_diffs)]

    pub_dates <- pub_dates[which(as_date(pub_dates) <= as_date(final_version) &
                                   as_date(pub_dates) >= as_date(start_date) - days(max_lag))]
    
    print(paste0("end_date: ", end_date, ", final_version: ", final_version))
    
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
          # M(A)D to final estimate (constant lag)...
          
          {R_pub <- R_est_ts %>%
            dplyr::select(starts_with("R_pub_") | date) %>%
            column_to_rownames("date") %>%
            rename_with(~ substr(.x, 7, 16)) %>%
            dplyr::select(all_of(available_pub_dates))
          diff_final <- switch_final <- switch2_final <- data.frame(estimated_after = make_difftime(day = seq(0, max_lag),
                                                                   units = "day")) %>%
            column_to_rownames(var = "estimated_after")
          
          for (i_realtime in colnames(R_pub)[as_date(colnames(R_pub)) <= as_date(end_date) - days_until_final]){
            time_diffs <- as_date(available_pub_dates) - as_date(i_realtime)
            time_diffs[time_diffs<days_until_final] <- 1000 # "big-M" in order to only consider dates at least days_until_final days later 
            i_final <- available_pub_dates[which.min(time_diffs)]
            
            diff <- R_pub[,i_final] - R_pub[,i_realtime]
            i_diff <- as.character(as_date(i_realtime) - as_date(rownames(R_pub[which(!is.na(diff)),])))
            diff_final[i_diff, i_realtime] <- diff %>% na.omit()
            
            switch <- as.numeric((R_pub[,i_final] > 1) != (R_pub[,i_realtime] > 1))
            i_switch <- as.character(as_date(i_realtime) - as_date(rownames(R_pub[which(!is.na(switch)),])))
            switch_final[i_switch, i_realtime] <- switch %>% na.omit()
            
            switch2 <- as.numeric(((R_pub[,i_final] > 1.03) &(R_pub[,i_realtime] < 0.97)) |
                                   ((R_pub[,i_final] < 0.97) & (R_pub[,i_realtime] > 1.03)))
            i_switch2 <- as.character(as_date(i_realtime) - as_date(rownames(R_pub[which(!is.na(switch2)),])))
            switch2_final[i_switch2, i_realtime] <- switch2 %>% na.omit()
            
          }
          
          # extract relevant rows
          idx <- intersect(rownames(diff_final), as.character(min_lag:max_lag))
          diff_final <- diff_final[idx,]
          # drop columns containing any NA
          diff_final <- diff_final[, colSums(is.na(diff_final)) == 0]
          
          df_diff_to_final[method, idx] <- rowMeans(diff_final)
          df_abs_diff_to_final[method, idx] <- rowMeans(abs(diff_final))
          df_abs_diff_to_final[method, "num_est"] <- df_diff_to_final[method, "num_est"] <- ncol(diff_final)
          
          # extract relevant rows for switch
          idx_switch <- intersect(rownames(switch_final), as.character(min_lag:max_lag))
          switch_final <- switch_final[idx_switch,]
          # drop columns containing any NA
          switch_final <- switch_final[, colSums(is.na(switch_final)) == 0]
          
          df_switch_final[method, idx_switch] <- rowMeans(switch_final)
          df_switch_final[method, "num_est"] <- ncol(switch_final)
          
          
          # extract relevant rows for switch2
          idx_switch2 <- intersect(rownames(switch2_final), as.character(min_lag:max_lag))
          switch2_final <- switch2_final[idx_switch2,]
          # drop columns containing any NA
          switch2_final <- switch2_final[, colSums(is.na(switch2_final)) == 0]
          
          df_switch2_final[method, idx_switch2] <- rowMeans(switch2_final)
          df_switch2_final[method, "num_est"] <- ncol(switch2_final)
          }
          
          if(method=="ETHZ_sliding_window"){
            diff_final_old <- diff_final %>% dplyr::select(1:"2021-01-25")
            diff_final_new <- diff_final %>% dplyr::select("2021-01-26":dim(diff_final)[2])
            df_diff_to_final["ETH_old", idx] <- rowMeans(diff_final_old)
            df_diff_to_final["ETH_new", idx] <- rowMeans(diff_final_new)
            df_abs_diff_to_final["ETH_old", idx] <- rowMeans(abs(diff_final_old))
            df_abs_diff_to_final["ETH_new", idx] <- rowMeans(abs(diff_final_new))
            df_abs_diff_to_final["ETH_old", "num_est"] <- df_diff_to_final["ETH_old", "num_est"] <- ncol(diff_final_old)
            df_abs_diff_to_final["ETH_new", "num_est"] <- df_diff_to_final["ETH_new", "num_est"] <- ncol(diff_final_new)
          }
          
          dates_used_for_diff_to_final <- colnames(diff_final)

        
                    
          ######
          # CI coverage rates and width...
          
          # if CIs are available calculate mean coverage rates and widths
          if (sum(!is.na(R_est_ts %>% dplyr::select(starts_with("lower")))) > 0){
            R_covered <- data.frame(date = R_est_ts$date)
            R_covered_difftime <-
              CI_width <-
              data.frame(estimated_after = make_difftime(day = seq(max_lag, min_lag),
                                                         units = "day"))
  
            for (pd in dates_used_for_diff_to_final){
              time_diffs <- as_date(available_pub_dates) - as_date(pd)
              time_diffs[time_diffs<days_until_final] <- 1000 # "big-M" in order to only consider dates at least days_until_final days later 
              date_final <- available_pub_dates[which.min(time_diffs)]
              i_final <- paste0("R_pub_", as_date(date_final))
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
            
            R_covered_difftime <- R_covered_difftime %>%
              column_to_rownames(var = "estimated_after")
            CI_width <- CI_width %>%
              column_to_rownames(var = "estimated_after")
              
            df_CI_coverage[method, rownames(R_covered_difftime)] <- rowMeans(R_covered_difftime)
            df_CI_coverage[method, "min_lag"] <- min_lag
            df_CI_coverage[method, "num_CIs"] <- dim(R_covered_difftime)[2]
            df_CI_coverage[method, "min_pub_date"] <- min(colnames(R_covered_difftime))
            df_CI_coverage[method, "max_pub_date"] <- max(colnames(R_covered_difftime))
            
            df_CI_width[method, rownames(CI_width)] <- rowMeans(CI_width)
            df_CI_width[method, "min_lag"] <- min_lag
            df_CI_width[method, "num_CIs"] <- dim(CI_width)[2]
            df_CI_width[method, "min_pub_date"] <- min(colnames(CI_width))
            df_CI_width[method, "max_pub_date"] <- max(colnames(CI_width))
          }
          
          if(method=="ETHZ_sliding_window"){
            R_covered_difftime_old <- R_covered_difftime %>% dplyr::select(1:"2021-01-25")
            R_covered_difftime_new <- R_covered_difftime %>% dplyr::select("2021-01-26":dim(diff_final)[2])
            df_CI_coverage["ETH_old", rownames(R_covered_difftime_old)] <- rowMeans(R_covered_difftime_old)
            df_CI_coverage["ETH_new", rownames(R_covered_difftime_new)] <- rowMeans(R_covered_difftime_new)
            df_CI_coverage["ETH_old", "num_est"] <- ncol(R_covered_difftime_old)
            df_CI_coverage["ETH_new", "num_est"] <- ncol(R_covered_difftime_new)
            
            CI_width_old <- CI_width %>% dplyr::select(1:"2021-01-25")
            CI_width_new <- CI_width %>% dplyr::select("2021-01-26":dim(diff_final)[2])
            df_CI_width["ETH_old", idx] <- rowMeans(CI_width_old)
            df_CI_width["ETH_new", idx] <- rowMeans(CI_width_new)
            df_CI_width["ETH_old", "num_est"] <- ncol(diff_final_old)
            df_CI_width["ETH_new", "num_est"] <- ncol(diff_final_new)
          }
          
          
          ######
          # extract labels...
          
          {labels <- R_est_ts %>%
            dplyr::select(date | any_of(paste0("label_", dates_used_for_diff_to_final))) %>%
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
            summarise(label = getmode(label), num_est = n())

          idx <- as.character(intersect(0:max_lag, labels_mode$lag))
          df_labels[method, idx] <- labels_mode %>%
            dplyr::filter(as.character(lag) %in% idx) %>%
            arrange(lag) %>%
            dplyr::select(label) %>% pull()
          df_labels[method, "num_est"] <- min((labels_mode %>% dplyr::filter(as.character(lag) %in% idx))[, "num_est"])
          }
        }
      }
    } else {
      print(paste("No estimates from", method, "for", country))
    }
  }
  
  return(list(CI_coverage = df_CI_coverage, 
              CI_width = df_CI_width,
              abs_diff_to_final = df_abs_diff_to_final, 
              diff_to_final = df_diff_to_final,
              switch_final = df_switch_final,
              switch2_final = df_switch2_final,
              labels = df_labels))
}

methods <- c("Braunschweig", "epiforecasts", "ETHZ_sliding_window",
             "globalrt_7d", "ilmenau", "RKI_7day", "rtlive", "SDSC")

# methods <- c("rtlive")

for (d in c(50,70,80)[2]){
  CI_eval <- calc_consistence_metrics(methods, days_until_final=d)
  path <- paste0("Rt_estimate_reconstruction/otherFiles/consistence_measures/", d, "/")
  write.csv(CI_eval$CI_coverage, paste0(path, "95_CI_coverage.csv"))
  write.csv(CI_eval$CI_width, paste0(path, "95_CI_width.csv"))
  write.csv(CI_eval$abs_diff_to_final, paste0(path, "abs_diff_to_final.csv"))
  write.csv(CI_eval$diff_to_final, paste0(path, "diff_to_final.csv"))
  write.csv(CI_eval$switch_to_final, paste0(path, "switch_to_final.csv"))
  write.csv(CI_eval$switch2_to_final, paste0(path, "switch2_to_final.csv"))
  write.csv(CI_eval$labels, paste0(path, "estimate_labels.csv"))
}

source("Rt_estimate_reconstruction/prepared_plots.R")

d <- 70
plot1 <- plot_CI_coverage_rates(days_until_final = d); plot1
plot2 <- plot_CI_widths(days_until_final = d); plot2
plot3 <- plot_diff_final(diff_type = "abs_diff", days_until_final = d); plot3
plot4 <- plot_diff_final(diff_type = "diff", days_until_final = d, ylim = c(-0.012,0.06)); plot4
plot5 <- plot_diff_final(diff_type = "switch", days_until_final = d, ylim = c(-0.012,0.4)); plot5
plot6 <- plot_diff_final(diff_type = "switch2", days_until_final = d, ylim = c(-0.012,0.4)); plot6

# consistence plot for main manuscript:
consistence_plot <- ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol=2, nrow=3,
                              labels = list("A", "B", "C", "D", "E", "F"),
                              font.label = list(size = 18, face = "bold", family = "serif"),
                              common.legend = T, legend="bottom",
                              legend.grob = get_legend(plot3))
print(consistence_plot)
ggsave(consistence_plot, filename = paste0("Figures/CI/consistence_plots_", d, ".pdf"),
       bg = "transparent", width = 16, height = 13.6)


# # supplementary figure on switching above / below 0 or [0.97, 1.03]
# switch_plot <- ggarrange(plot5, plot6, ncol=2, nrow=1,
#                               labels = list("A", "B"),
#                               font.label = list(size = 18, face = "bold", family = "serif"),
#                               common.legend = T, legend="bottom",
#                               legend.grob = get_legend(plot3))
# print(switch_plot)
# ggsave(switch_plot, filename = paste0("Figures/CI/switch_plots_", d, ".pdf"),
#        bg = "transparent", width = 16, height = 6)


