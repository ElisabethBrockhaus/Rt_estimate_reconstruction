setwd("../..")
# needs to be the directory with the repos "Rt_estimate_reconstruction"
getwd()

source("Rt_estimate_reconstruction/calculate_estimates.R")


# save generation time distributions for rtlive
gt_main_analysis <- get_infectivity_profile(gt_type="gamma", gt_mean = 4, gt_sd = 4, n_days = 30) # consenus for comparison
gt_lower_example <- get_infectivity_profile(gt_type="gamma", gt_mean = 3.4, gt_sd = 1.8, n_days = 30) # AGES
gt_upper_example <- get_infectivity_profile(gt_type="gamma", gt_mean = 5.6, gt_sd = 4.2, n_days = 30) # Ilmenau
gt_distributions <- cbind(0:30, gt_main_analysis, gt_lower_example, gt_upper_example)
colnames(gt_distributions) <- c("days_after_transmission", "gamma 4 (4)", "gamma 3.4 (1.8)", "gamma 5.6 (4.2)")

write.csv(gt_distributions, "Rt_estimate_reconstruction/rtlive/gt_distributions.csv", row.names = F)

plot(gt_main_analysis, type="l", ylim = c(0, 0.32))
lines(gt_lower_example, col="red")
lines(gt_upper_example, col="blue")
#lines(get_infectivity_profile(gt_type="exponential", gt_mean = 4.7, gt_sd = 2.9, n_days = 30), col = "green")
lines(c(0.00000000e+00, 7.33148180e-03, 1.03868478e-01, 1.95147643e-01,
        1.93797118e-01, 1.52828476e-01, 1.09576993e-01, 7.54229927e-02,
        5.10921821e-02, 3.44849062e-02, 2.33407210e-02, 1.58958932e-02,
        1.09121047e-02, 7.55711343e-03, 5.28165515e-03, 3.72529635e-03,
        2.65131301e-03, 1.90352772e-03, 1.37822529e-03, 1.00599866e-03,
        7.40014799e-04, 5.48403039e-04, 4.09289125e-04, 3.07532780e-04,
        2.32567920e-04, 1.76961382e-04, 1.35443284e-04, 1.04249379e-04,
        8.06714362e-05, 6.27474302e-05), col="orange")
legend(x="topright", legend=c("gamma_4", "gamma_3.4", "gamma_5.6", "lnorm_discr_4.7"),
       lty=1, col=c("black", "red", "blue", "orange"))
