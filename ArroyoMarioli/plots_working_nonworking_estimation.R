path <- "Rt_estimate_reconstruction/ArroyoMarioli/"

# estimates
file_base <- "estimates/estimated_R_globalrt_data_from_rtlive"
long <- read_csv(paste0(path, file_base, ".csv"))
short <- read_csv(paste0(path, file_base, "_short.csv"))
dates <- seq(min(long$Date), max(long$Date), by="month")

plot(short$Date, short$R, type="l",
     col="red", xlab="date", ylab="estimated R",
     xlim=c(min(long$Date), max(long$Date)), xaxt="n")
axis(1, dates, format(dates, "%b %Y"), cex.axis=.9)
abline(h=1, col="grey", lty=2)
lines(long$Date, long$R)
lines(long$Date, long$ci_95_u, lty=2)
lines(long$Date, long$ci_95_l, lty=2)
legend(x="topright", legend=c("on long time series", "95% confidence bounds", "on shorter time series"),
       lty=c(1,2,1), col=c("black", "black", "red"))

# growth rates
file_base <- "input_output_dataset/dataset_rtlive"
long <- read_csv(paste0(path, file_base, ".csv"))
short <- read_csv(paste0(path, file_base, "_short.csv"))
dates <- seq(min(long$Date), max(long$Date), by="month")

plot(long$Date, long$gr_infected_4, type="l",
     xlab="date", ylab="growth rate", xaxt="n")
axis(1, dates, format(dates, "%b %Y"), cex.axis=.9)
abline(h=0, col="grey", lty=2)
lines(short$Date, short$gr_infected_4, col="red")
legend(x="topright", legend=c("growth rate in both datasets","additional for longer time series"),
       lty=1, col=c("red", "black"))
