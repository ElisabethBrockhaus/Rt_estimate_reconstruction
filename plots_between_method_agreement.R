library(readr)
library(dplyr)

setwd("..")
getwd()

source("Rt_estimate_reconstruction/prepared_plots.R")

# load estimates
path_estimates <- "Rt_estimate_reconstruction/estimates/"
{
  estimates_window <- read_csv(paste0(path_estimates, "R_cmp_Window_2021-07-10.csv"))
  estimates_gtd <- read_csv(paste0(path_estimates, "R_cmp_GTD_2021-07-10.csv"))
  estimates_input <- read_csv(paste0(path_estimates, "R_cmp_Input_2021-11-23.csv"))
  estimates_preprocess <- read_csv(paste0(path_estimates, "R_cmp_Preprocess_2021-07-10.csv"))
  estimates_SD_gtd <- read_csv(paste0(path_estimates, "R_cmp_GTD_SD_2021-07-10.csv"))
  estimates_pub <- read_csv(paste0(path_estimates, "R_pub_2021-07-10.csv"))
  estimates_adjInput <- read_csv(paste0(path_estimates, "R_adjInput_2021-07-10.csv"))
  estimates_adjInputWindow <- read_csv(paste0(path_estimates, "R_adjInputWindow_2021-07-10.csv"))
  estimates_adjInputWindowGTD <- read_csv(paste0(path_estimates, "R_adjInputWindowGTD_2021-07-10.csv"))
  estimates_adjInputWindowGTDDelay <- read_csv(paste0(path_estimates, "R_adjInputWindowGTDDelays_2021-07-10.csv"))
  estimates_adjInputWindowGTDIncRep <- read_csv(paste0(path_estimates, "R_adjInputWindowGTDIncRep_2021-07-10.csv"))
  estimates_adjInputWindowGTDIncRepType <- read_csv(paste0(path_estimates, "R_adjInputWindowGTDIncRepType_2021-07-10.csv"))
  estimates_pub_ci <- read_csv(paste0(path_estimates, "R_pub_ci_2021-07-10.csv"))
  estimates_adjInput_ci <- read_csv(paste0(path_estimates, "R_adjInput_ci_2021-07-10.csv"))
  estimates_adjInputWindow_ci <- read_csv(paste0(path_estimates, "R_adjInputWindow_ci_2021-07-10.csv"))
  estimates_adjInputWindowGTD_ci <- read_csv(paste0(path_estimates, "R_adjInputWindowGTD_ci_2021-07-10.csv"))
  estimates_adjInputWindowGTDDelay_ci <- read_csv(paste0(path_estimates, "R_adjInputWindowGTDDelays_ci_2021-07-10.csv"))
  estimates_adjInputWindowGTDIncRep_ci <- read_csv(paste0(path_estimates, "R_adjInputWindowGTDIncRep_ci_2021-07-10.csv"))
  estimates_adjInputWindowGTDIncRepType_ci <- read_csv(paste0(path_estimates, "R_adjInputWindowGTDIncRepType_ci_2021-07-10.csv"))
}

# functions for extracting ylims
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMin <- function(data) sapply(data, min, na.rm = TRUE)


##########################################
# Influence of parameters and input data #
##########################################

# window size
{
  methods_window <- c("Ilmenau", "ETH", "SDSC", "RKI")
  windows <-         c(1,         3,     4,      7)
  names(windows) <- methods_window
  window_strs <- c()
  
  for (src in names(windows)){
    window_str <- paste0(windows[src], ", ", src)
    window_strs <- c(window_strs, window_str)
  }
  # find ylim
  ylim_window_l <- min(colMin(estimates_window %>%
                                dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>%
                                dplyr::select(ends_with(as.character(window_strs)))))
  ylim_window_u <- max(colMax(estimates_window %>%
                                dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>%
                                dplyr::select(ends_with(as.character(window_strs)))))
  # plot
  plot_for_comparison(estimates_window, comp_methods = window_strs,
                      col_palette = "YlOrRd", name_consensus = "7, RKI",
                      legend_name = "window size", filenames = "_influence_window.pdf",
                      sort_numerically = TRUE,  plot_diff_matrices=T,
                      ylim_l = ylim_window_l, ylim_u = ylim_window_u)
}

# generation time distribution
{
  methods_gtd <- c("epiforecasts", "RKI",      "consensus", "rtlive",  "ETH/SDSC", "Ilmenau", "globalrt", "HZI")
  distrs <-      c("gamma",        "constant", "gamma",     "lognorm", "gamma",    "ad hoc",  "gamma",    "?")
  means <-       c( 3.6,            4.0,        4.0,         4.7,       4.8,        5.6,       7.0,        10.3)
  sds <-         c( 3.1,            0.0,        4.0,         2.9,       2.3,        4.2,       7.0,        7.6)
  gtds <- cbind("type"=distrs, "mean"=means, "sd"=sds)
  rownames(gtds) <- methods_gtd
  gtd_strs <- c()
  for (src in rownames(gtds)){
    gtd_str <- paste0(gtds[src, "mean"], "(", gtds[src, "sd"], "), ", src)
    gtd_strs <- c(gtd_strs, gtd_str)
  }
  # find ylim
  ylim_gtd_l <- min(colMin(estimates_gtd %>%
                             dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>%
                             dplyr::select(ends_with(gtd_strs))))
  ylim_gtd_u <- max(colMax(estimates_gtd %>%
                             dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>%
                             dplyr::select(ends_with(gtd_strs))))
  # plot
  plot_for_comparison(estimates_gtd, comp_methods = gtd_strs,
                      col_palette = "YlGn", name_consensus = "4(4), consensus",
                      legend_name = "GTD", filenames = "_influence_GTD.pdf",
                      sort_numerically = TRUE, plot_diff_matrices=T,
                      ylim_l = ylim_gtd_l, ylim_u = ylim_gtd_u)
  
  
  # standard deviation of the GTD 
  sds <- c(3.1, 0.001, 4.0, 2.9, 2.3, 4.2, 7.0, 7.6)
  # plot
  plot_for_comparison(estimates_SD_gtd, comp_methods = as.character(sds),
                      col_palette = "YlGn", name_consensus = 4.0,
                      legend_name = "SD of the GTD", filenames = "_influence_SD_GTD.pdf",
                      sort_numerically = TRUE, plot_diff_matrices=T,
                      ylim_l = ylim_gtd_l, ylim_u = ylim_gtd_u)
}

# input data
{
  data_sources <- c("JHU", "RKI, positive test", "RKI, symptom onset", "WHO")
  # find ylim
  ylim_input_l <- min(colMin(estimates_input %>% dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>% dplyr::select(ends_with(data_sources))))
  ylim_input_u <- max(colMax(estimates_input %>% dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>% dplyr::select(ends_with(data_sources))))
  # plot
  plot_for_comparison(estimates_input, comp_methods = data_sources,
                      col_palette = "incidence data", name_consensus = "RKI, positive test",
                      legend_name = "data source", filenames = "_influence_input_data.pdf",
                      sort_numerically = FALSE, plot_diff_matrices=T,
                      ylim_l = ylim_input_l, ylim_u = ylim_input_u)
}

# preprocessing
{
  preprocessing <- c("none", "RKI", "SDSC", "ETH")
  # find ylim
  ylim_preprocess_l <- min(colMin(estimates_preprocess %>% dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>% dplyr::select(ends_with(preprocessing))))
  ylim_preprocess_u <- max(colMax(estimates_preprocess %>% dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>% dplyr::select(ends_with(preprocessing))))
  # plot
  plot_for_comparison(estimates_preprocess, comp_methods = preprocessing,
                      col_palette = "incidence data", name_consensus = "none",
                      legend_name = "preprocessing", filenames = "_influence_preprocessing.pdf",
                      sort_numerically = FALSE, plot_diff_matrices=T,
                      ylim_l = ylim_preprocess_l, ylim_u = ylim_preprocess_u)
}

# method
{
  model <- c("RKI", "ETH", "Ilmenau", "SDSC", "globalrt", "epiforecasts", "rtlive")
  estimates_model <- estimates_adjInputWindowGTDDelay %>% dplyr::select(!c("consensus", "HZI"))
  # find ylim
  ylim_model_l <- min(colMin(estimates_model %>% dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>% dplyr::select(ends_with(model))))
  ylim_model_u <- max(colMax(estimates_model %>% dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>% dplyr::select(ends_with(model))))
  # plot
  plot_for_comparison(estimates_model, comp_methods = model,
                      col_palette = "Dark2", name_consensus = "none",
                      legend_name = "method", filenames = "_influence_model.pdf",
                      sort_numerically = FALSE, plot_diff_matrices=T,
                      ylim_l = ylim_model_l, ylim_u = ylim_model_u)
}



####################################
# Influence of the choice of model #
####################################

org_methods <- c("RKI", "ETH", "Ilmenau", "SDSC", "globalrt", "epiforecasts", "rtlive", "HZI")
org_methods_CI <- c("RKI", "ETH", "Ilmenau", "SDSC", "globalrt", "epiforecasts", "rtlive")
comp_methods <- c("consensus", "RKI", "SDSC", "ETH",
                  "Ilmenau", "epiforecasts", "globalrt", "rtlive", "HZI")
comp_CI <- c("consensus", "SDSC", "ETH",
             "Ilmenau", "epiforecasts", "globalrt", "rtlive")

# find ylim
{
  max_pub <- max(colMax(estimates_pub %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(!starts_with("date"))))
  max_adjInput <- max(colMax(estimates_adjInput %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(!starts_with("date"))))
  max_adjInputWindow <- max(colMax(estimates_adjInputWindow %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(!starts_with("date"))))
  max_adjInputWindowGTD <- max(colMax(estimates_adjInputWindowGTD %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(!starts_with("date"))))
  max_adjInputWindowGTDDelay <- max(colMax(estimates_adjInputWindowGTDDelay %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(!starts_with("date"))))
  max_pub_ci <- max(colMax(estimates_pub_ci %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInput_ci <- max(colMax(estimates_adjInput_ci %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInputWindow_ci <- max(colMax(estimates_adjInputWindow_ci %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInputWindowGTD_ci <- max(colMax(estimates_adjInputWindowGTD_ci %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInputWindowGTDDelay_ci <- max(colMax(estimates_adjInputWindowGTDDelay_ci %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  ylim_u <- max(max_pub, max_adjInput, max_adjInputWindow, max_adjInputWindowGTD, max_adjInputWindowGTDDelay,
                max_pub_ci, max_adjInput_ci, max_adjInputWindow_ci, max_adjInputWindowGTD_ci, max_adjInputWindowGTDDelay_ci)
}

source("Rt_estimate_reconstruction/prepared_plots.R")

# plot median
{
  plot_for_comparison(estimates_pub, org_methods, legend_name = "method",
                      filenames = "_real-time.pdf", ylim_u = ylim_u,
                      plot_diff_matrices=T, include_consensus=F)
  plot_for_comparison(estimates_adjInput, comp_methods, legend_name = "method",
                      filenames = "_adjInput.pdf", ylim_u = ylim_u,
                      plot_diff_matrices=T, include_consensus=F)
  plot_for_comparison(estimates_adjInputWindow, comp_methods, legend_name = "method",
                      filenames = "_adjInputWindow.pdf",
                      ylim_l = 0.5, ylim_u = 1.5,
                      plot_diff_matrices=T, include_consensus=F)
  plot_for_comparison(estimates_adjInputWindowGTD, comp_methods, legend_name = "method",
                      filenames = "_adjInputWindowGTD.pdf",
                      ylim_l = 0.5, ylim_u = 1.5,
                      plot_diff_matrices=T, include_consensus=F)
  plot_for_comparison(estimates_adjInputWindowGTDDelay, comp_methods, legend_name = "method",
                      filenames = "_adjInputWindowGTDDelays.pdf",
                      ylim_l = 0.5, ylim_u = 1.5,
                      plot_diff_matrices=T, include_consensus=F)
}




# plot CI
{
  plot_for_comparison(estimates_pub_ci, org_methods_CI, include_CI = T, legend_name = "method",
                      filenames = "_CI_real-time.pdf", ylim_u = ylim_u,
                      include_consensus=F, plot_width_diff_matrices=F)
  plot_for_comparison(estimates_adjInput_ci, comp_CI, include_CI = T,
                      filenames = "_CI_adjInput.pdf", ylim_u = ylim_u,
                      include_consensus=F, plot_width_diff_matrices=F)
  plot_for_comparison(estimates_adjInputWindow_ci, comp_CI, include_CI = T, legend_name = "method",
                      filenames = "_CI_adjInputWindow.pdf",
                      ylim_l = 0.5, ylim_u = 1.5,
                      include_consensus=F, plot_width_diff_matrices=F)
  plot_for_comparison(estimates_adjInputWindowGTD_ci, comp_CI, include_CI = T, legend_name = "method",
                      filenames = "_CI_adjInputWindowGTD.pdf",
                      ylim_l = 0.5, ylim_u = 1.5,
                      include_consensus=F, plot_width_diff_matrices=F)
  plot_for_comparison(estimates_adjInputWindowGTDDelay_ci, comp_CI, include_CI = T, legend_name = "method",
                      filenames = "_CI_adjInputWindowGTDDelays.pdf",
                      ylim_l = 0.5, ylim_u = 1.5,
                      include_consensus=F, plot_width_diff_matrices=F)
}





############
# APPENDIX #
############

# find ylim for longer time period
{
  max_pub <- max(colMax(estimates_pub %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(-date)))
  max_adjInput <- max(colMax(estimates_adjInput %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(-date)))
  max_adjInputWindow <- max(colMax(estimates_adjInputWindow %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(-date)))
  max_adjInputWindowGTD <- max(colMax(estimates_adjInputWindowGTD %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(-date)))
  max_adjInputWindowGTDDelay <- max(colMax(estimates_adjInputWindowGTDDelay %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(-date)))
  max_pub_ci <- max(colMax(estimates_pub_ci %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInput_ci <- max(colMax(estimates_adjInput_ci %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInputWindow_ci <- max(colMax(estimates_adjInputWindow_ci %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInputWindowGTD_ci <- max(colMax(estimates_adjInputWindowGTD_ci %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInputWindowGTDDelay_ci <- max(colMax(estimates_adjInputWindowGTDDelay_ci %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  ylim_u_long <- max(max_pub, max_adjInput, max_adjInputWindow, max_adjInputWindowGTD, max_adjInputWindowGTDDelay,
                max_pub_ci, max_adjInput_ci, max_adjInputWindow_ci, max_adjInputWindowGTD_ci, max_adjInputWindowGTDDelay_ci)
}
ylim_u_long <- 4

# some plots over more days
{
  plot_for_comparison(estimates_pub, org_methods,
                      #start_absdiff = "2020-03-16",
                      start_date = "2020-03-01", end_date = "2021-07-10",
                      filenames = "_real-time_long.pdf",
                      long_time_frame = T,
                      ylim_u = ylim_u_long,
                      plot_diff_matrices = F,
                      include_consensus = F)
  plot_for_comparison(estimates_adjInput, comp_methods,
                      #start_absdiff = "2020-03-16",
                      start_date = "2020-03-01", end_date = "2021-07-10",
                      filenames = "_adjInput_long.pdf",
                      long_time_frame = T,
                      ylim_u = ylim_u_long,
                      plot_diff_matrices = F,
                      include_consensus = F)
  plot_for_comparison(estimates_adjInputWindow, comp_methods,
                      #start_absdiff = "2020-03-16",
                      start_date = "2020-03-01", end_date = "2021-07-10",
                      filenames = "_adjInputWindow_long.pdf",
                      long_time_frame = T,
                      ylim_u = ylim_u_long,
                      plot_diff_matrices = F,
                      include_consensus = F)
  plot_for_comparison(estimates_adjInputWindowGTD, comp_methods,
                      #start_absdiff = "2020-03-16", 
                      start_date = "2020-03-01", end_date = "2021-07-10",
                      filenames = "_adjInputWindowGTD_long.pdf",
                      long_time_frame = T,
                      ylim_u = ylim_u_long,
                      plot_diff_matrices = F,
                      include_consensus = F)
  plot_for_comparison(estimates_adjInputWindowGTDDelay, comp_methods,
                      #start_absdiff = "2020-03-16", 
                      start_date = "2020-03-01", end_date = "2021-07-10",
                      filenames = "_adjInputWindowGTDDelays_long.pdf",
                      long_time_frame = T,
                      ylim_u = ylim_u_long,
                      plot_diff_matrices = F,
                      include_consensus = F)
}


# plot median with old delay adjustments
{
  BA_methods <- c("consensus", "RKI", "SDSC", "ETH", "Ilmenau", "epiforecasts", "globalrt", "rtlive")
  plot_for_comparison(estimates_adjInputWindowGTDIncRep, BA_methods, legend_name = "method",
                      filenames = "_adjInputWindowGTDIncRep.pdf",
                      ylim_l = 0.5, ylim_u = 1.5,
                      plot_diff_matrices=T, include_consensus=F)
  plot_for_comparison(estimates_adjInputWindowGTDIncRepType, BA_methods, legend_name = "method",
                      filenames = "_adjInputWindowGTDIncRepType.pdf",
                      ylim_l = 0.5, ylim_u = 1.5,
                      plot_diff_matrices=T, include_consensus=F)
}
