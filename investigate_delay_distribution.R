library(readr)
setwd("..")
getwd()

#######
# ETH #
#######
########### start copied ###########
# load functions
source("Rt_estimate_reconstruction/ETH/otherScripts/2_utils_getInfectionIncidence.R")
# load parameter
source("Rt_estimate_reconstruction/ETH/otherScripts/2_params_InfectionIncidencePars.R")
# load empirical delays
delays_data_path <- "Rt_estimate_reconstruction/ETH/data/all_delays.csv"
delays_onset_to_count <- read_csv(delays_data_path,
                                  col_types = cols(
                                    data_type = col_character(),
                                    onset_date = col_date(format = ""),
                                    count_date = col_date(format = ""),
                                    delay = col_number()))

# constant delay distribution
constant_delay_distributions <- list()
for (type_i in unique(names(shape_onset_to_count))) {
  m <- get_vector_constant_waiting_time_distr(
    shape_incubation,
    scale_incubation,
    shape_onset_to_count[[type_i]],
    scale_onset_to_count[[type_i]])
  
  constant_delay_distributions <- c(constant_delay_distributions, list(m))
}
names(constant_delay_distributions) <- unique(names(shape_onset_to_count))

constant_delay_symptom_to_report_distributions <- list()
for (type_i in unique(names(shape_onset_to_count))) {
  m <- get_vector_constant_waiting_time_distr(
    0,
    0,
    shape_onset_to_count[[type_i]],
    scale_onset_to_count[[type_i]])
  
  constant_delay_symptom_to_report_distributions <- c(constant_delay_symptom_to_report_distributions, list(m))
}
names(constant_delay_symptom_to_report_distributions) <- paste0('Onset to ',  unique(names(shape_onset_to_count)))

constant_delay_distributions <- c(constant_delay_distributions, constant_delay_symptom_to_report_distributions)
########### end copied ###########

View(constant_delay_distributions)

ETH_incubation <- constant_delay_distributions$Symptoms
incubation_mean <- 0:199 %*% ETH_incubation
incubation_mean/7
incubation_var <- ETH_incubation %*% ((0:199) - c(incubation_mean))^2
sqrt(incubation_var)/7
incubation_shape <- (incubation_mean^2)/incubation_var
incubation_scale <- incubation_var/incubation_mean

ETH_report_delays <- constant_delay_distributions$`Onset to Confirmed cases`
report_mean <- 0:199 %*% ETH_report_delays
report_mean/7
report_var <- ETH_report_delays %*% ((0:199) - c(report_mean))^2
sqrt(report_var)/7
report_shape <- (report_mean^2)/report_var
report_scale <- report_var/report_mean

plot(ETH_incubation[1:30], type="l")
lines(ETH_report_delays)

lines(dgamma(0:199, shape = incubation_shape, scale = incubation_scale), lty=2, col="red")
lines(dgamma(0:199, shape = report_shape, scale = report_scale), lty=2, col="red")

incubation <- function(x) dgamma(x, shape = incubation_shape, scale = incubation_scale)
reporting <- function(y) dgamma(y, shape = report_shape, scale = report_scale)
# convolution integral
delay <- function(z) integrate(function(x,z) reporting(z-x)*incubation(x),-Inf,Inf,z)$value
delay <- Vectorize(delay)
z <- 0:199
lines(z,delay(z),lty=2,col="blue")

# parameters of resulting convolution
mean <- z %*% delay(z)
mean
sd <- sqrt(delay(z) %*% (z - c(mean))^2)
sd

##########
# rtlive #
##########
rtlive_delays <- read_csv("rtlive/p_delay.csv")
rtlive_delay <- rtlive_delays$p_delay[7:66]
delay_mean <- rtlive_delay %*% (1:60)
delay_var <- rtlive_delay %*% ((1:60) - c(delay_mean))^2
delay_sd <- sqrt(delay_var)
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
sd <- sqrt(delay(z) %*% ((0:65) - c(mean))^2)
sd

################
# epiforecasts #
################
incubation_mean <- 1.62
incubation_var <- 0.418^2
reporting_mean <- 0.832
reporting_var <- 1.44^2
inc_shape <- (incubation_mean^2)/incubation_var
inc_scale <- incubation_var/incubation_mean
report_shape <- (reporting_mean^2)/reporting_var
report_scale <- reporting_var/reporting_mean

incubation <- function(x) dgamma(x, shape = inc_shape, scale = inc_scale)
reporting <- function(y) {
  dist <- dgamma(y, shape = report_shape, scale = report_scale)
  dist[dist==Inf] <- 1
  return(dist)
}
# convolution integral
delay <- function(z) integrate(function(x,z) reporting(z-x)*incubation(x),-Inf,Inf,z)$value
delay <- Vectorize(delay)
z <- 0:60
plot(z,delay(z), type="l", xlim=c(0,20))

# parameters of resulting convolution
mean <- z %*% delay(z)
mean
sd <- delay(z) %*% (z - c(mean))^2
sd



