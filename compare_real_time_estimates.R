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
# methods <- list.dirs(path_estimates, full.names = F, recursive = F)
# methods <- methods[!methods %in% c("", "AW_7day", "AW_WVday", "owid", "ETHZ_step",
#                                    "ETHZ_sliding_window_deaths", "ETHZ_step_deaths",
#                                    "zidatalab")]
methods <- c("Braunschweig", "epiforecasts", "ETHZ_sliding_window",
             "globalrt_7d", "ilmenau", "RKI_7day", "rtlive", "SDSC")
methods

available_countries <- read.csv("Rt_estimate_reconstruction/otherFiles/available_countries.csv", row.names = 1)
pub_delays <- read.csv("Rt_estimate_reconstruction/otherFiles/pub_delays.csv", row.names = 1)


#################################
# plot estimates as time series #
#################################

#start <- "2021-04-01"
#end <- "2021-05-01"

#start <- "2021-02-16"
#end <- "2021-03-18"

# period for which RKI was criticized to correct always upwards
start <- "2020-09-28"
end <- "2020-12-07"

plots <- list()

for (method in methods[methods!="globalrt_7d"]){
  print(method)
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)

  final_version <- "2021-07-16"
  if (method == "Braunschweig") final_version <- "2021-07-20"
  if (method == "epiforecasts") final_version <- as.character(as.Date(start) + 106)
  
  pub_dates <- pub_dates[which(pub_dates <= end &
                                 pub_dates >= start)]
  end_date <- as_date(max(pub_dates))
  pub_dates <- c(pub_dates, final_version)
  for (country in c("DE", "AT", "CH")[1]){
    print(country)
    if (available_countries[method, country]) {
      if (exists("R_est_ts")) rm(R_est_ts)
      if (exists("R_est")) rm(R_est)
      for (pub_date in pub_dates){
        tryCatch(
          {
            R_est <- load_published_R_estimates(method,
                                                end = end_date,
                                                pub_date = pub_date,
                                                location = country,
                                                verbose = F) %>%
              dplyr::select(date, R_pub)
            if (pub_date != final_version){
              last <- max(R_est[rowSums(!is.na(R_est))>1, "date"])
              R_est <- R_est %>% dplyr::filter(date <= last, date > last - 7)
            }
            names(R_est) <- c("date", paste0("R.", pub_date))
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
        
        folder_ending <- paste0("_realtime_raw/", end, "/")
        folder <- paste0("Figures/estimates", folder_ending)
        if (!file.exists(folder)) {
          dir.create(folder)
        }
        
        if (method == "ilmenau"){
          ylim <- c(0.3, 2.15)
        } else if(method == "Braunschweig"){
          ylim <- c(0.9, 1.75)
        } else {
          ylim <- c(0.9, 1.5)
        }
        
        plots[[method]] <- plot_real_time_estimates(R_est_ts,
                                                    start_date = start,
                                                    end_date = end,
                                                    plot_title = method,
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

plots_arranged <- ggarrange(plots[[methods[1]]], plots[[methods[2]]], plots[[methods[3]]],
                            plots[[methods[5]]], plots[[methods[6]]],
                            plots[[methods[7]]], plots[[methods[8]]],
                            ncol=1, nrow=7,
                            common.legend = T, legend="bottom", legend.grob = get_legend(plots[["epiforecasts"]]))
print(plots_arranged)
ggsave(plots_arranged, filename = paste0("Figures/estimates", folder_ending,
                                         "all_methods.pdf"),
       bg = "transparent", width = 21, height = 29.7)



#####################################################################################
# plot estimates over days between estimation and target day averaged over weekdays #
#####################################################################################

# dates over which the mean is calculated
start_date <- as_date("2020-11-16")
end_date <- as_date("2021-05-01")
target_dates <- seq(start_date, end_date, by = "day")
wds <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")

min_lag <- 1
max_lag <- 30

for (method in methods){
  print(method)
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which((as_date(pub_dates) <= end_date) &
                                 (as_date(pub_dates) >= start_date))]
  for (country in c("DE", "AT", "CH")[1]){
    print(country)
    if (available_countries[method, country]) {
      if (exists("R_est_ts")) rm(R_est_ts)
      if (exists("R_est")) rm(R_est)
      for (pub_date in pub_dates){
        tryCatch(
          {
            R_est <- load_published_R_estimates(method,
                                                start = max(as_date(min(pub_dates)) - pub_delays[method, country],
                                                            as_date(pub_date) - max_lag),
                                                end = min(as_date(max(pub_dates)) - pub_delays[method, country] - 6,
                                                          as_date(pub_date) - min_lag),
                                                pub_date = pub_date,
                                                location = country,
                                                verbose = F) %>%
              dplyr::select("date", "R_pub") %>%
              mutate(estimated_after = as_date(pub_date) - date) %>%
              reshape(idvar = "estimated_after", timevar = "date", direction = "wide") %>%
              rename_with(~ gsub(".", "_", .x, fixed = TRUE)) %>%
              dplyr::select("estimated_after",
                            which((colnames(.)[-1] %>%
                                     substr(7, 16) %>%
                                     as_date()) %in% target_dates) + 1)
          },
          error = function(e) {
            R_est <<- data.frame(estimated_after = make_difftime(day = seq(min_lag, max_lag),
                                                                 units = "day"))
            }
        )
        
        if (!exists("R_est_ts")){
          R_est_ts <- R_est
        } else {
          R_est_ts <- R_est_ts %>% merge(R_est, all = T) %>%
            group_by(estimated_after) %>%
            summarise_all(list(~ mean(., na.rm = T)))
          R_est_ts <- R_est_ts[order(colnames(R_est_ts))]
        }
      }
      if (exists("R_est_ts")){
        
        min_lag_plot <- pub_delays[method, country]
        max_lag_plot <- min(min_lag_plot + 6, max_lag)
        
        R_est_ts <- R_est_ts %>%
          dplyr::filter((estimated_after >= min_lag_plot) &
                          (estimated_after <= max_lag_plot))
        
        # extract weekdays from target dates
        org_names <- colnames(R_est_ts)[2:(dim(R_est_ts)[2])]
        colnames(R_est_ts)[2:(dim(R_est_ts)[2])] <- weekdays(org_names %>%
                                                                 substr(7, 16) %>%
                                                                 as_date()) %>%
          paste0(org_names)
        
        # calculate mean over weekdays
        R_est_means <- R_est_ts["estimated_after"]
        for (wd in wds) {
          R_est_means[wd] <- R_est_ts %>%
            dplyr::select(starts_with(wd)) %>%
            apply(1, function(x) exp(mean(log(x), na.rm = TRUE)))
        }

        R_est_plot <- R_est_means %>% rename(date = estimated_after)
        
        plot_for_comparison(R_est_plot,
                            comp_methods = wds,
                            start_date = make_difftime(day = min_lag_plot),
                            end_date = make_difftime(day = max_lag_plot),
                            ylim_l=0.5, ylim_u=1.5,
                            legend_name = "mean estimate for",
                            plot_title = paste(method, country),
                            col_palette = "Spectral",
                            filenames = paste0("_realtime_corrections_per_weekday/",
                                               method, "_", country, ".pdf"),
                            verbose = F)
      }
    } else {
      print(paste("No estimates from", method, "for", country, "."))
    }
  }
}



############################
# weekday effects averaged #
############################

# pub_dates
start_date <- as_date("2020-11-16")
end_date <- as_date("2021-05-01")
wds <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")

for (method in methods){
  print(method)
  
  if (method == "Braunschweig") {
    print("Method handled manually due to irregular pub dates.")
    next
  }
  
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which((as_date(pub_dates) <= end_date) &
                                 (as_date(pub_dates) >= start_date))]
    
  for (country in c("DE", "AT", "CH")[1]){
    print(country)
    
    min_lag <- pub_delays[method, country]
    max_lag <- min_lag + 6
    
    if (available_countries[method, country]) {
      if (exists("R_est_list")) rm(R_est_list)
      if (exists("R_est")) rm(R_est)
      
      R_est_empty <- data.frame(target_weekday = wds,
                                Monday = rep(NA, 7),
                                Tuesday = rep(NA, 7),
                                Wednesday = rep(NA, 7),
                                Thursday = rep(NA, 7),
                                Friday = rep(NA, 7),
                                Saturday = rep(NA, 7),
                                Sunday = rep(NA, 7))
      
      for (pub_date in pub_dates){
        tryCatch(
          {
            R_est <- load_published_R_estimates(method,
                                                start = max(as_date(min(pub_dates)) - min_lag,
                                                            as_date(pub_date) - max_lag),
                                                end = min(as_date(max(pub_dates)) - max_lag,
                                                          as_date(pub_date) - min_lag),
                                                pub_date = pub_date,
                                                location = country,
                                                verbose = F) %>%
              mutate(target_weekday = weekdays(date)) %>%
              dplyr::select("target_weekday", "R_pub")
            R_est <- R_est[match(wds, R_est$target_weekday),]
          },
          error = function(e) {R_est <<- data.frame(target_weekday = wds,
                                                    R_pub = rep(NA, 7))}
        )
        
        weekday <- weekdays(as_date(pub_date))
        R_est <- R_est_empty %>% mutate(!!weekday := R_est$R_pub)
          
        if (!exists("R_est_list")){
          R_est_list <- R_est
        } else {
          R_est_list <- rbindlist(list(R_est_list, R_est))
        }
      }
      if (exists("R_est_list")){
        
        R_est_mean <- R_est_list[, lapply(.SD, function(x) exp(mean(log(x), na.rm = TRUE))),
                                 by = .(target_weekday)]
        
        plot_weekday_effects(R_est_mean,
                             ylim_l=0.75, ylim_u=1.25,
                             plot_title = paste0("Mean over latest 7 estimates per pub date (",
                                                 method, " ", country, ")"),
                             min_lag = min_lag, max_lag = max_lag,
                             filenames = paste0(method, "_", country, ".pdf"),
                             verbose = F)
      }
    } else {
      print(paste("No estimates from", method, "for", country, "."))
    }
  }
}

# Braunschweig: mean over pub weekdays, otherwise estimates not published often enough
{
  method <- "Braunschweig"
  country <- "DE"
  print(method)
  print(country)
    
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which((as_date(pub_dates) <= end_date) &
                                 (as_date(pub_dates) >= start_date))]
  
  
    
  min_lag <- pub_delays[method, country]
  max_lag <- min_lag + 6
  
  if (exists("R_est_list")) rm(R_est_list)
  if (exists("R_est")) rm(R_est)
  
  R_est_empty <- data.frame(target_weekday = wds,
                            R_pub = rep(NA, 7))
  
  for (pub_date in pub_dates){
    tryCatch(
      {
        R_est <- load_published_R_estimates(method,
                                            start = max(as_date(min(pub_dates)) - min_lag,
                                                        as_date(pub_date) - max_lag),
                                            end = min(as_date(max(pub_dates)) - max_lag,
                                                      as_date(pub_date) - min_lag),
                                            pub_date = pub_date,
                                            location = country,
                                            verbose = F) %>%
          mutate(target_weekday = weekdays(date)) %>%
          dplyr::select("target_weekday", "R_pub")
        R_est <- R_est[match(wds, R_est$target_weekday),]
      },
      error = function(e) {R_est <<- data.frame(target_weekday = wds,
                                                R_pub = rep(NA, 7))}
    )
    
    if (!exists("R_est_list")){
      R_est_list <- R_est
    } else {
      R_est_list <- rbindlist(list(R_est_list, R_est))
    }
  }
  if (exists("R_est_list")){
    
    R_est_list <- R_est_list %>% na.omit()
    
    R_est_mean <- R_est_list[, lapply(.SD, function(x) exp(mean(log(x), na.rm = TRUE))),
                             by = .(target_weekday)]
    
    plot_title <- paste0("Mean over latest 7 estimates per pub date (",
                         method, " ", country, ")")
    filenames <- paste0(method, "_", country, ".pdf")
    ylim_l <- 0.75
    ylim_u <- 1.25

    R_est <- R_est_mean %>%
      rename(mean_R = R_pub)
    
    wds <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
    
    # plot
    R_plot <- ggplot(data = R_est, aes(x = factor(target_weekday,
                                                  levels = wds),
                                       y = mean_R,
                                       group = 1)) +
      geom_line() +
      geom_hline(aes(yintercept = 1))
    
    R_plot <- R_plot +
      labs(x = "weekday (target date)", y = "Rt estimate",
           subtitle = paste(min_lag, "-", max_lag, "days previous to pub date")) +
      scale_x_discrete()
    
    R_plot <- R_plot +
      theme_minimal() +
      theme(
        plot.title = element_text(size=18),
        axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        axis.line = element_line(),
        axis.line.y.right = element_line(),
        axis.line.x.top = element_line(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent")
      ) +
      ggtitle(plot_title)
    
    R_plot <- R_plot +
      coord_cartesian(ylim = c(ylim_l, ylim_u), expand = FALSE)
    
    ggsave(R_plot, filename = paste0("Figures/estimates_realtime_target_date_influence/", filenames),  bg = "transparent",
           width = 13.1, height = 5.8)
    print(R_plot)
  }
}
