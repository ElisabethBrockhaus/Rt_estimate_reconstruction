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
                                   "ETHZ_sliding_window_deaths", "ETHZ_step_deaths",
                                   "zidatalab")]
methods

available_countries <- read.csv("Rt_estimate_reconstruction/otherFiles/available_countries.csv", row.names = 1)
pub_delays <- read.csv("Rt_estimate_reconstruction/otherFiles/pub_delays_mode.csv", row.names = 1)


#############################################
# plot coverage of 95% confidence intervals #
#############################################

start_date <- as_date("2020-11-16")
end_date <- as_date("2021-05-01")
country <- "DE"

n <- length(methods)
CI_coverage <- data.frame(matrix(rep(NA, n*24), nrow = n), row.names = methods)
colnames(CI_coverage) <- c(0:20, "num_CIs", "min_pub_date", "max_pub_date")

for (method in methods){
  print(method)
  pub_dates <- list.files(paste0(path_estimates, method),
                          full.names = F) %>% substr(1, 10)
  pub_dates <- pub_dates[which(as_date(pub_dates) <= end_date &
                                 as_date(pub_dates) >= start_date)]
  
  final_version <- "2021-07-16"
  if (method == "Braunschweig") final_version <- "2021-07-18"
  
  if (available_countries[method, country]) {
    if (exists("R_est_ts")) rm(R_est_ts)
    if (exists("R_est")) rm(R_est)
    
    min_lag <- pub_delays[method, country]
    max_lag <- min_lag + 6
    
    R_est_ts <- load_published_R_estimates(method,
                                           start = start_date - max_lag,
                                           end = as_date(final_version),
                                           pub_date = final_version,
                                           location = country,
                                           verbose = F) %>%
      dplyr::select("date", "R_pub")
    names(R_est_ts)[2] <- "R_final"
    
    for (pub_date in pub_dates){
      tryCatch(
        {
          R_est <- load_published_R_estimates(method,
                                              start = as_date(pub_date) - max_lag,
                                              end = as_date(pub_date) - min_lag,
                                              pub_date = pub_date,
                                              location = country,
                                              verbose = F) %>%
            dplyr::select("date", "lower", "upper")
          
          names(R_est) <- c("date",
                            paste0("02.5q_", pub_date),
                            paste0("97.5q_", pub_date))
        },
        error = function(e) {R_est <<- data.frame(date = seq(as_date("2019-12-28"),
                                                             as_date(pub_date),
                                                             by = "day"))}
      )
      R_est_ts <- R_est_ts %>% full_join(R_est, by="date")
    }
    if (exists("R_est_ts")){
      available_pub_dates <- colnames(R_est_ts[3:dim(R_est_ts)[2]]) %>%
        substr(7, 16) %>%
        unique()
      
      R_covered <- data.frame(date = R_est_ts$date)
      R_covered_difftime <- data.frame(estimated_after = make_difftime(day = seq(min_lag, max_lag),
                                                                       units = "day"))
      
      for (pd in available_pub_dates){
        R_covered[pd] <- ifelse(is.na(R_est_ts[, paste0("02.5q_", pd)]), NA,
                                ifelse((R_est_ts$R_final >= R_est_ts[, paste0("02.5q_", pd)]) &
                                         (R_est_ts$R_final <= R_est_ts[, paste0("97.5q_", pd)]),
                                       TRUE,
                                       FALSE))
        if (dim(na.omit(R_covered[pd]))[1] == 7) {
          R_covered_difftime[pd] <- na.omit(R_covered[pd])
        }
      }
      
      ea <- R_covered_difftime$estimated_after
      CI_coverage[method, as.character(ea)] <- rowMeans(R_covered_difftime[,-1])
      CI_coverage[method, "num_CIs"] <- dim(R_covered_difftime)[2]
      CI_coverage[method, "min_pub_date"] <- min(colnames(R_covered_difftime[,-1]))
      CI_coverage[method, "max_pub_date"] <- max(colnames(R_covered_difftime[,-1]))

    }
  } else {
    print(paste("No estimates from", method, "for", country))
  }
}

View(CI_coverage)
write.csv(CI_coverage, "Rt_estimate_reconstruction/otherFiles/CI_coverage.csv")


CI_coverage <- read_csv("Rt_estimate_reconstruction/otherFiles/CI_coverage.csv") %>%
  as.data.frame() %>%
  column_to_rownames("...1")

coverage_data <- CI_coverage[,1:21] %>%
  rownames_to_column("method") %>%
  dplyr::filter(!(method %in% methods[c(1,5,6,7,9,10)])) %>% # drop AGES and most of globalrt variations
  mutate(method = plyr::mapvalues(method,
                                  c("ETHZ_sliding_window", "globalrt_7d", "ilmenau", "RKI_7day"),
                                  c("ETH", "globalrt", "Ilmenau", "RKI"))) %>%
  arrange(method)

coverage_data <- coverage_data %>%
  gather("variable", "value", 2:dim(coverage_data)[2]) %>%
  mutate(variable = -1 * as.numeric(variable))

coverage_plot <- ggplot() +
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
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("95%-CI coverage rates") +
  labs(x = "target date - pub date", y = "coverage rate")

col_values <- get_colors(methods = unique(coverage_data$method), palette = "methods")

coverage_plot <- coverage_plot + 
  geom_line(data=coverage_data,
            aes(x = variable, y = value, color = method),
            size = .8, na.rm = T) +
  scale_color_manual(values=col_values, name="method")

print(coverage_plot)


