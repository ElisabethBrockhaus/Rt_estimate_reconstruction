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
                    start_date = "2021-03-10", end_date = "2021-07-10",
                    legend_name = "research group", filenames = "_real-time_epiforecasts.pdf")
