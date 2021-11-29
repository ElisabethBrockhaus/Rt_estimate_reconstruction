path <- "Rt_estimate_reconstruction/epiforecasts/estimates/"

estimates_org <- load_published_R_estimates("epiforecasts")
colnames(estimates_org)[2] <- "R_calc"

estimates_adj_data <- qread(paste0(path,"R_calc_2021-07-10_final_adjInput.qs"))[,c("date", "mean", "lower_95", "upper_95")]
names(estimates_adj_data) <- c("date", "R_calc", "lower", "upper")

estimates_adj_data_gt <- qread(paste0(path,"R_calc_2021-07-10_final_adjInputWindowGTD.qs"))[,c("date", "mean", "lower_95", "upper_95")]
names(estimates_adj_data_gt) <- c("date", "R_calc", "lower", "upper")

estimates_adj_data_gt_delay <- qread(paste0(path,"R_calc_2021-07-10_final_adjAll.qs"))[,c("date", "mean", "lower_95", "upper_95")]
names(estimates_adj_data_gt_delay) <- c("date", "R_calc", "lower", "upper")

estimates_joined <- estimates_org[,c("date", "R_calc")] %>%
  full_join(estimates_adj_data[,c("date", "R_calc")], by = "date") %>% 
  full_join(estimates_adj_data_gt[,c("date", "R_calc")], by = "date") %>% 
  full_join(estimates_adj_data_gt_delay[,c("date", "R_calc")], by = "date")

methods <- c("org", "adj data", "adj data & gt", "adj data, gt & delay")
plot_for_comparison(estimates_joined, methods,
                    legend_name = "type", filenames = "_TEST.pdf")

