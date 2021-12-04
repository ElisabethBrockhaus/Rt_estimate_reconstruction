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


plot_for_comparison(estimates_pub, org_methods, legend_name = "research group",
                    col_palette = "Set1", start_col=1,
                    filenames = "_real-time.pdf", ylims = ylim_diffWindow, plot_diff_matrices=T)
plot_for_comparison(estimates_adjInput, comp_methods,
                    col_palette = "Set1", start_col=0,
                    filenames = "_adjInput.pdf", ylims = ylim_diffWindow, plot_diff_matrices=T)
plot_for_comparison(estimates_adjInputWindow, comp_methods,
                    col_palette = "Set1", start_col=0,
                    filenames = "_adjInputWindow.pdf", ylims = ylim_diffWindow, plot_diff_matrices=T)
plot_for_comparison(estimates_adjInputWindowGTD, comp_methods,
                    col_palette = "Set1", start_col=0,
                    filenames = "_adjInputWindowGTD.pdf", ylims = ylim_diffWindow, plot_diff_matrices=T)
plot_for_comparison(estimates_adjAll, comp_methods,
                    col_palette = "Set1", start_col=0,
                    filenames = "_adjAll.pdf", ylims = ylim_diffWindow, plot_diff_matrices=T)

plot_for_comparison(estimates_adjInputWindow, comp_methods,
                    col_palette = "Set1", start_col=0,
                    filenames = "_adjInputWindow_smallYlim.pdf", ylims = ylim_adjWindow)
plot_for_comparison(estimates_adjInputWindowGTD, comp_methods,
                    col_palette = "Set1", start_col=0,
                    filenames = "_adjInputWindowGTD_smallYlim.pdf", ylims = ylim_adjWindow)
plot_for_comparison(estimates_adjAll, comp_methods,
                    col_palette = "Set1", start_col=0,
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
                    col_palette = "Set1", start_col=1,
                    filenames = "_CI_real-time.pdf", ylims = ylim_CI)
plot_for_comparison(estimates_adjInput_CI, comp_CI, include_CI = T,
                    col_palette = "Set1", start_col=0, skip_RKI=T,
                    filenames = "_CI_adjInput.pdf", ylims = ylim_CI)
plot_for_comparison(estimates_adjInputWindow_CI, comp_CI, include_CI = T,
                    col_palette = "Set1", start_col=0, skip_RKI=T,
                    filenames = "_CI_adjInputWindow.pdf", ylims = ylim_CI)
plot_for_comparison(estimates_adjInputWindowGTD_CI, comp_CI, include_CI = T,
                    col_palette = "Set1", start_col=0, skip_RKI=T,
                    filenames = "_CI_adjInputWindowGTD.pdf", ylims = ylim_CI)
plot_for_comparison(estimates_adjAll_CI, comp_CI, include_CI = T,
                    col_palette = "Set1", start_col=0, skip_RKI=T,
                    filenames = "_CI_adjAll.pdf", ylims = ylim_CI)

