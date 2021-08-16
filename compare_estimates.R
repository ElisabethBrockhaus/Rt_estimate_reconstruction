wd <- "D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code"
setwd(wd)

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

##############
# RKI
##############

# load data
RKI_data <- load_RKI_data()

# estimation
RKI_est <- estimate_RKI_R(RKI_data)
RKI_EpiEstim_est <- estimate_RKI_R_EpiEstim(RKI_data)

# plots for comparison
plot_published_vs_calculated(RKI_data, RKI_est, method_name="RKI")
plot_published_vs_calculated(RKI_data, RKI_EpiEstim_est, method_name="RKI with EpiEstim")



##############
# ETH
##############

# load data
ETH_data <- load_ETH_data()

# estimation
ETH_est <- estimate_ETH_R(ETH_data)

# plots for comparison
plot_published_vs_calculated(published=data.frame(dates=ETH_data$date,
                                                  R=ETH_data$R),
                             calculated=ETH_est,
                             method_name="ETH")



##############
# Ilmenau
##############

# load data
Ilmenau_data <- load_Ilmenau_data()

# estimation
Ilmenau_est <- estimate_Ilmenau_R(Ilmenau_data)

# plots for comparison
plot_published_vs_calculated(Ilmenau_data, Ilmenau_est, method_name="Ilmenau")



##############
# epiforecasts
##############
source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/calculate_estimates.R")

# load data
epiforecasts_data <- load_epiforecasts_data()

# estimation
EpiNow2_est <- estimate_EpiNow2_R(epiforecasts_data)
EpiNow2_est <- out$estimates$summarised
EpiNow2_est <- EpiNow2_est[EpiNow2_est$variable=="R", c("date", "mean")]
EpiNow2_est <- EpiNow2_est[EpiNow2_est$date<="2020-12-14",]

# plots for comparison
plot_published_vs_calculated(epiforecasts_data, EpiNow2_est$mean, method_name="EpiNow2")



