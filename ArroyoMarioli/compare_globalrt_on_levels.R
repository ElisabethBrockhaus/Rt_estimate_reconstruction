library(readr)

setwd("../..")
source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

path <- "Rt_estimate_reconstruction/ArroyoMarioli/estimates/"

estimates_org <- load_published_R_estimates("globalrt_7d")
names(estimates_org)[2] <- "R_calc"

file_org_rebuilt <- "bayesian_smoother_7_org.csv"
estimates_org_rebuilt <- read_csv(paste0(path, file_org_rebuilt))
names(estimates_org_rebuilt) <- c("date", "R_calc", "lower", "upper")

file_adj_data <- "bayesian_smoother_7.csv"
estimates_adj_data <- read_csv(paste0(path, file_adj_data))
names(estimates_adj_data)[1:4] <- c("date", "R_calc", "lower", "upper")

file_adj_data_gt <- "bayesian_smoother_4.csv"
estimates_adj_data_gt <- read_csv(paste0(path, file_adj_data_gt))
names(estimates_adj_data_gt)[1:4] <- c("date", "R_calc", "lower", "upper")

estimates_joined <- estimates_org[,c("date", "R_calc")] %>%
  full_join(estimates_org_rebuilt[,c("date", "R_calc")], by = "date") %>% 
  full_join(estimates_adj_data[,c("date", "R_calc")], by = "date") %>% 
  full_join(estimates_adj_data_gt[,c("date", "R_calc")], by = "date")

methods <- c("org", "org rebuilt", "adj data", "adj data & gt")
plot_for_comparison(estimates_joined, methods,
                    legend_name = "type", filenames = "_TEST.pdf")
