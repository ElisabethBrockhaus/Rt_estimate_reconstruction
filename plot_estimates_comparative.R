source("Rt_estimate_reconstruction/prepared_plots.R")

org_methods <- c("RKI", "ETH", "Ilmenau", "SDSC", "Zi",
                 "globalrt", "epiforecasts", "rtlive")
plot_for_comparison(estimates_pub, org_methods, filenames = "_real-time.pdf",
                    method = "original", variation = "parameters")


comp_methods <- c("EpiEstim", "ETH", "AGES", "Ilmenau", "epiforecasts", "globalrt", "rtlive")
plot_for_comparison(estimates_input, comp_methods, filenames = "_input-data.pdf",
                    method = "rtlive", variation = "input data")


plot_for_comparison(estimates_GTD, comp_methods, filenames = "_gtd.pdf",
                    method = method, variation = "generation time")


plot_for_comparison(estimates_delays, comp_methods, filenames = "_delays.pdf",
                    method = method, variation = "delays")


plot_for_comparison(estimates_window_Cori, comp_methods, filenames = "_windowCori.pdf",
                    method = method, variation = "generation time")


plot_for_comparison(estimates_window, comp_methods, filenames = "_windowAll.pdf",
                    method = method, variation = "generation time")


path <- "Rt_estimate_reconstruction/rtlive/summaries/"
file <- "DE_2021-07-31_all_trace_summary.csv"
R_rtlive_gt_4.7 <- read.csv(paste0(path, file)) %>%
  dplyr::select(X, mean, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "mean", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
file <- "gamma_4_de_2021-07-10_all_summary.csv"
R_rtlive_gt_4 <- read.csv(paste0(path, file)) %>%
  dplyr::select(X, mean, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "mean", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
file <- "gamma_3.4_de_2021-07-10_all_summary.csv"
R_rtlive_gt_3.4 <- read.csv(paste0(path, file)) %>%
  dplyr::select(X, mean, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "mean", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))

path <- "Rt_estimate_reconstruction/rtlive/summaries/"
file <- "gamma_5.6_de_2021-07-10_all_summary.csv"
R_rtlive_gt_5.6 <- read.csv(paste0(path, file)) %>%
  dplyr::select(X, mean, X2.5., X97.5.) %>%
  rename(c("date" = "X", "R_calc" = "mean", "lower" = "X2.5.", "upper" = "X97.5.")) %>%
  mutate(date = as_date(date))

estimates_rtlive_gt <- R_rtlive_gt_4.7 %>%
  full_join(R_rtlive_gt_4, by = "date") %>% 
  full_join(R_rtlive_gt_3.4, by = "date") %>% 
  full_join(R_rtlive_gt_5.6, by = "date")

comp_methods <- c("rtlive", "rtlive_4", "rtlive_3.4", "rtlive_5.6")
plot_for_comparison(estimates_rtlive_gt, comp_methods, filenames = "_rtlive-various-gtds.pdf",
                    method = "rtlive", variation = "generation time")
