source("Rt_estimate_reconstruction/prepared_plots.R")

epiforecasts_07_10 <- load_published_R_estimates("epiforecasts")
epiforecasts_07_08 <- load_published_R_estimates("epiforecasts", pub_date = "2021-07-08")
epiforecasts_07_07 <- load_published_R_estimates("epiforecasts", pub_date = "2021-07-07")
epiforecasts_07_06 <- load_published_R_estimates("epiforecasts", pub_date = "2021-07-06")
epiforecasts_07_03 <- load_published_R_estimates("epiforecasts", pub_date = "2021-07-03")
epiforecasts_07_02 <- load_published_R_estimates("epiforecasts", pub_date = "2021-07-02")
epiforecasts_07_01 <- load_published_R_estimates("epiforecasts", pub_date = "2021-07-01")

# merge estimates and plot for comparison
estimates_epiforecasts <- epiforecasts_07_10[,c("date", "R_pub")] %>%
  full_join(epiforecasts_07_08[,c("date", "R_pub")], by = "date") %>% 
  full_join(epiforecasts_07_07[,c("date", "R_pub")], by = "date") %>% 
  full_join(epiforecasts_07_06[,c("date", "R_pub")], by = "date") %>% 
  full_join(epiforecasts_07_03[,c("date", "R_pub")], by = "date") %>% 
  full_join(epiforecasts_07_02[,c("date", "R_pub")], by = "date") %>% 
  full_join(epiforecasts_07_01[,c("date", "R_pub")], by = "date")

org_methods <- c("10", "8", "7", "6", "3", "2", "1")
plot_for_comparison(estimates_epiforecasts, org_methods,
                    legend_name = "research group", filenames = "_real-time_epiforecasts.pdf")




# find ylim
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMin <- function(data) sapply(data, min, na.rm = TRUE)

# for estimates with different window sizes
max_pub <- max(colMax(estimates_pub %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_pub"))))
min_pub <- min(colMin(estimates_pub %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_pub"))))
max_adjInput <- max(colMax(estimates_adjInput %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
min_adjInput <- min(colMin(estimates_adjInput %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
ylim_diffWindow <- c(min(min_pub, min_adjInput), max(max_pub, max_adjInput))

# for estimates with adjusted window sizes
max_adjInputWindow <- max(colMax(estimates_adjInputWindow %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
min_adjInputWindow <- min(colMin(estimates_adjInputWindow %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
max_adjInputWindowGTD <- max(colMax(estimates_adjInputWindowGTD %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
min_adjInputWindowGTD <- min(colMin(estimates_adjInputWindowGTD %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
max_adjAll <- max(colMax(estimates_adjAll %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
min_adjAll <- min(colMin(estimates_adjAll %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
ylim_adjWindow <- c(min(min_adjInputWindow, min_adjInputWindowGTD, min_adjAll), max(max_adjInputWindow, max_adjInputWindowGTD, max_adjAll))


source("Rt_estimate_reconstruction/prepared_plots.R")
plot_for_comparison(estimates_pub, org_methods, legend_name = "research group",
                    filenames = "_real-time.pdf", ylims = ylim_diffWindow, plot_diff_matrices=T)
plot_for_comparison(estimates_adjInput, comp_methods,
                    filenames = "_adjInput.pdf", ylims = ylim_diffWindow, plot_diff_matrices=T)
plot_for_comparison(estimates_adjInputWindow, comp_methods,
                    filenames = "_adjInputWindow.pdf", ylims = ylim_diffWindow, plot_diff_matrices=T)
plot_for_comparison(estimates_adjInputWindowGTD, comp_methods,
                    filenames = "_adjInputWindowGTD.pdf", ylims = ylim_diffWindow, plot_diff_matrices=T)
plot_for_comparison(estimates_adjAll, comp_methods,
                    filenames = "_adjAll.pdf", ylims = ylim_diffWindow, plot_diff_matrices=T)

plot_for_comparison(estimates_adjInputWindow, comp_methods,
                    filenames = "_adjInputWindow_smallYlim.pdf", ylims = ylim_adjWindow)
plot_for_comparison(estimates_adjInputWindowGTD, comp_methods,
                    filenames = "_adjInputWindowGTD_smallYlim.pdf", ylims = ylim_adjWindow)
plot_for_comparison(estimates_adjAll, comp_methods,
                    filenames = "_adjAll_smallYlim.pdf", ylims = ylim_adjWindow)

# same for CI plots
max_pub <- max(colMax(estimates_pub_ci %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
min_pub <- min(colMin(estimates_pub_ci %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("lower"))))
max_adjInput <- max(colMax(estimates_adjInput %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
min_adjInput <- min(colMin(estimates_adjInput %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("lower"))))
max_adjInputWindow <- max(colMax(estimates_adjInputWindow %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
min_adjInputWindow <- min(colMin(estimates_adjInputWindow %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("lower"))))
max_adjInputWindowGTD <- max(colMax(estimates_adjInputWindowGTD %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
min_adjInputWindowGTD <- min(colMin(estimates_adjInputWindowGTD %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("lower"))))
max_adjAll <- max(colMax(estimates_adjAll %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
min_adjAll <- min(colMin(estimates_adjAll %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("lower"))))
ylim_CI <- c(min(min_pub, min_adjInput, min_adjInputWindow, min_adjInputWindowGTD, min_adjAll),
             max(max_pub, max_adjInput, max_adjInputWindow, max_adjInputWindowGTD, max_adjAll))


plot_for_comparison(estimates_pub_ci, org_methods, include_CI = T, legend_name = "research group",
                    filenames = "_CI_real-time.pdf", ylims = ylim_CI)
plot_for_comparison(estimates_adjInput_CI, comp_CI, include_CI = T,
                    filenames = "_CI_adjInput.pdf", ylims = ylim_CI)
plot_for_comparison(estimates_adjInputWindow_CI, comp_CI, include_CI = T,
                    filenames = "_CI_adjInputWindow.pdf", ylims = ylim_CI)
plot_for_comparison(estimates_adjInputWindowGTD_CI, comp_CI, include_CI = T,
                    filenames = "_CI_adjInputWindowGTD.pdf", ylims = ylim_CI)
plot_for_comparison(estimates_adjAll_CI, comp_CI, include_CI = T,
                    filenames = "_CI_adjAll.pdf", ylims = ylim_CI)

