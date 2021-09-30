library(readr)
setwd("..")
getwd()
ETH_delays <- read_csv("Rt_estimate_reconstruction/ETH/data/DEU/DEU_data_delays.csv")
delay_mean <- mean(ETH_delays$delay)
delay_var <- var(ETH_delays$delay)
shape <- (delay_mean^2)/delay_var
scale <- delay_var/delay_mean

hist(ETH_delays$delay, breaks=0:30, freq=F, ylim = c(0, 0.25))
lines(dgamma(1:30, shape = shape, scale = scale), lty=2, col="red")

incubation <- function(x) dgamma(x, shape = shape/4, scale = scale)
reporting <- function(y) dgamma(y, shape = shape*3/4, scale = scale)
# convolution integral
delay <- function(z) integrate(function(x,z) reporting(z-x)*incubation(x),-Inf,Inf,z)$value
delay <- Vectorize(delay)
z <- seq(0,30,1)
lines(z,delay(z),lty=2,col="blue")

incubation_mean <- (shape/4)*scale
incubation_mean
incubation_sd <- sqrt((shape/4)*scale^2)
incubation_sd
reporting_mean <- (shape*3/4)*scale
reporting_mean
reporting_sd <- sqrt((shape*3/4)*scale^2)
reporting_sd


rtlive_delays <- read_csv("rtlive/p_delay.csv")
rtlive_delay <- rtlive_delays$p_delay[7:66]
delay_mean <- rtlive_delay %*% (1:60)
delay_var <- rtlive_delay %*% ((1:60) - c(delay_mean))^2
#delay_var <- delay_var - 7
shape <- ((delay_mean)^2)/delay_var
scale <- delay_var/(delay_mean)

plot(rtlive_delays$p_delay[2:31], type="l")
#barplot(rtlive_delays$p_delay, ylim = c(0, 0.15), xlim = c(0,65), names.arg=0:65)
lines(5:30, dgamma(0:25, shape = shape, scale = scale), lty=2, col="red")

#lines(dgamma(0:65, shape = (5^2)/0.1, scale = 0.1/5), col="red")

incubation <- function(x) dgamma(x, shape = (5^2)/0.1, scale = 0.1/5)
reporting <- function(y) dgamma(y, shape = shape, scale = scale)
# convolution integral
delay <- function(z) integrate(function(x,z) reporting(z-x)*incubation(x),-Inf,Inf,z)$value
delay <- Vectorize(delay)
z <- seq(0,65,1)
lines(z,delay(z),lty=2,col="blue")

# parameters of resulting convolution
mean <- 0:65 %*% delay(z)
mean
sd <- delay(z) %*% ((0:65) - c(mean))^2
sd
