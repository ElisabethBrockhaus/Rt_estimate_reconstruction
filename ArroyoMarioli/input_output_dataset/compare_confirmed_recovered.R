library(readr)

setwd("D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code")

# read data
confirmed_cases <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/time_series_covid19_confirmed_global.csv")
recovered_cases <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/time_series_covid19_recovered_global.csv")

# filter for german time series
confirmed_cases <- confirmed_cases[confirmed_cases$`Country/Region` == "Germany", 6:dim(confirmed_cases)[2]]
recovered_cases <- recovered_cases[recovered_cases$`Country/Region` == "Germany", 6:dim(recovered_cases)[2]]

# transpose, bind and reformat date
confirmed_cases <- t(confirmed_cases)
colnames(confirmed_cases) <- "confirmed"
recovered_cases <- t(recovered_cases)
colnames(recovered_cases) <- "recovered"
cases <- data.frame(date = as.Date(rownames(confirmed_cases), format = "%m/%d/%y"),
                    confirmed = confirmed_cases,
                    recovered = recovered_cases)
rownames(cases) <- 1:dim(cases)[1]

plot(cases$date, cases$confirmed*0.975, type="l")
lines(cases$date-17, cases$recovered, col="blue")

plot(cases$confirmed[1:500]*0.975 - cases$recovered[18:517], type="l")

# read dataset
data <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/dataset.csv")
# filter german data
data_DE <- data[data$`Country/Region` == "Germany",]

# plot new cases and new recovered
plot(data_DE$Date, data_DE$new_cases, type="l")
lines(data_DE$Date, data_DE$new_recovered, col="red")


estimates_KF <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/estimates/estimated_R_KF.csv")
estimates_KF <- estimates_KF[estimates_KF$`Country/Region` == "Germany",]
estimates_STAN <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/estimates/estimated_R_STAN.csv")
estimates_STAN <- estimates_STAN[estimates_STAN$`Country/Region` == "Germany",]
estimates <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/estimates/estimated_R.csv")
estimates <- estimates[estimates$`Country/Region` == "Germany",]
estimates_adj <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/estimates/estimated_R_adj.csv")
estimates_adj <- estimates_adj[estimates_adj$`Country/Region` == "Germany",]

plot(estimates_KF$Date, estimates_KF$R, type="l")
lines(estimates_KF$Date, estimates_KF$ci_95_l)
lines(estimates_KF$Date, estimates_KF$ci_95_u)
lines(estimates_STAN$Date, estimates_STAN$R, col="red")
lines(estimates_STAN$Date, estimates_STAN$ci_95_l, col="red")
lines(estimates_STAN$Date, estimates_STAN$ci_95_u, col="red")
lines(estimates$Date, estimates$R, col="blue")
lines(estimates$Date, estimates$ci_95_l, col="blue")
lines(estimates$Date, estimates$ci_95_u, col="blue")
lines(estimates_adj$Date, estimates_adj$R, col="darkgreen")
lines(estimates_adj$Date, estimates_adj$ci_95_l, col="darkgreen")
lines(estimates_adj$Date, estimates_adj$ci_95_u, col="darkgreen")
