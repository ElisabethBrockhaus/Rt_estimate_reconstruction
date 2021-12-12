# TODO: source saved estimates


# find ylim
colMax <- function(data) sapply(data, max, na.rm = TRUE)
{
  max_pub <- max(colMax(estimates_pub %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_pub"))))
  max_adjInput <- max(colMax(estimates_adjInput %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
  max_adjInputWindow <- max(colMax(estimates_adjInputWindow %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
  max_adjInputWindowGTD <- max(colMax(estimates_adjInputWindowGTD %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
  max_adjAll <- max(colMax(estimates_adjAll %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
  max_pub_ci <- max(colMax(estimates_pub_ci %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInput_ci <- max(colMax(estimates_adjInput %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInputWindow_ci <- max(colMax(estimates_adjInputWindow %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInputWindowGTD_ci <- max(colMax(estimates_adjInputWindowGTD %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjAll_ci <- max(colMax(estimates_adjAll %>% dplyr::filter(date>="2021-01-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  ylim_u <- max(max_pub, max_adjInput, max_adjInputWindow, max_adjInputWindowGTD, max_adjAll,
                max_pub_ci, max_adjInput_ci, max_adjInputWindow_ci, max_adjInputWindowGTD_ci, max_adjAll_ci)
}

source("Rt_estimate_reconstruction/prepared_plots.R")

# plot median
plot_for_comparison(estimates_pub, org_methods, legend_name = "research group",
                    filenames = "_real-time.pdf", ylim_u = ylim_u)
plot_for_comparison(estimates_adjInput, comp_methods,
                    filenames = "_adjInput.pdf", ylim_u = ylim_u, plot_diff_matrices=T)
plot_for_comparison(estimates_adjInputWindow, comp_methods,
                    filenames = "_adjInputWindow.pdf", ylim_u = ylim_u, plot_diff_matrices=T)
plot_for_comparison(estimates_adjInputWindowGTD, comp_methods,
                    filenames = "_adjInputWindowGTD.pdf", ylim_u = ylim_u, plot_diff_matrices=T)
plot_for_comparison(estimates_adjAll, comp_methods,
                    filenames = "_adjAll.pdf", ylim_u = ylim_u, plot_diff_matrices=T)


# plot CI
plot_for_comparison(estimates_pub_ci, org_methods, include_CI = T, legend_name = "research group",
                    filenames = "_CI_real-time.pdf", ylim_u = ylim_u)
plot_for_comparison(estimates_adjInput_CI, comp_CI, include_CI = T,
                    filenames = "_CI_adjInput.pdf", ylim_u = ylim_u)
plot_for_comparison(estimates_adjInputWindow_CI, comp_CI, include_CI = T,
                    filenames = "_CI_adjInputWindow.pdf", ylim_u = ylim_u)
plot_for_comparison(estimates_adjInputWindowGTD_CI, comp_CI, include_CI = T,
                    filenames = "_CI_adjInputWindowGTD.pdf", ylim_u = ylim_u)
plot_for_comparison(estimates_adjAll_CI, comp_CI, include_CI = T,
                    filenames = "_CI_adjAll.pdf", ylim_u = ylim_u)

# find ylim for longer time period
{
  max_pub <- max(colMax(estimates_pub %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_pub"))))
  max_adjInput <- max(colMax(estimates_adjInput %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
  max_adjInputWindow <- max(colMax(estimates_adjInputWindow %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
  max_adjInputWindowGTD <- max(colMax(estimates_adjInputWindowGTD %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
  max_adjAll <- max(colMax(estimates_adjAll %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("R_calc"))))
  max_pub_ci <- max(colMax(estimates_pub_ci %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInput_ci <- max(colMax(estimates_adjInput %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInputWindow_ci <- max(colMax(estimates_adjInputWindow %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjInputWindowGTD_ci <- max(colMax(estimates_adjInputWindowGTD %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  max_adjAll_ci <- max(colMax(estimates_adjAll %>% dplyr::filter(date>="2020-03-01", date<"2021-06-10") %>% dplyr::select(starts_with("upper"))))
  ylim_u_long <- max(max_pub, max_adjInput, max_adjInputWindow, max_adjInputWindowGTD, max_adjAll,
                max_pub_ci, max_adjInput_ci, max_adjInputWindow_ci, max_adjInputWindowGTD_ci, max_adjAll_ci)
}
ylim_u_long <- 4

# some plots over more days
source("Rt_estimate_reconstruction/prepared_plots.R")
plot_for_comparison(estimates_adjInput, comp_methods,
                    #start_absdiff = "2020-03-16",
                    start_date = "2020-03-01", end_date = "2021-07-10",
                    filenames = "_adjInput_long.pdf", ylim_u = ylim_u_long, plot_diff_matrices=F)
plot_for_comparison(estimates_adjInputWindow, comp_methods,
                    #start_absdiff = "2020-03-16",
                    start_date = "2020-03-01", end_date = "2021-07-10",
                    filenames = "_adjInputWindow_long.pdf", ylim_u = ylim_u_long, plot_diff_matrices=F)
plot_for_comparison(estimates_adjInputWindowGTD, comp_methods,
                    #start_absdiff = "2020-03-16", 
                    start_date = "2020-03-01", end_date = "2021-07-10",
                    filenames = "_adjInputWindowGTD_long.pdf", ylim_u = ylim_u_long, plot_diff_matrices=F)
plot_for_comparison(estimates_adjAll, comp_methods,
                    #start_absdiff = "2020-03-16", 
                    start_date = "2020-03-01", end_date = "2021-07-10",
                    filenames = "_adjAll_long.pdf", ylim_u = ylim_u_long, plot_diff_matrices=F)
