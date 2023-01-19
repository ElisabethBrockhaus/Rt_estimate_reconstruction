library(readr)
library(EpiNow2)
library(dplyr)
library(lubridate)

setwd("../..")
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

# empirical moments of delay distribution
mean((delays_onset_to_count %>% dplyr::filter(onset_date <= as_date("2021-06-01")))$delay)
sd((delays_onset_to_count %>% dplyr::filter(onset_date <= as_date("2021-06-01")))$delay)

# constant distribution (not used for Germany)
View(constant_delay_distributions)

ETH_incubation <- constant_delay_distributions$Symptoms
incubation_mean <- 0:199 %*% ETH_incubation
incubation_mean
incubation_var <- ETH_incubation %*% ((0:199) - c(incubation_mean))^2
sqrt(incubation_var)
incubation_shape <- (incubation_mean^2)/incubation_var
incubation_scale <- incubation_var/incubation_mean

ETH_report_delays <- constant_delay_distributions$`Onset to Confirmed cases`
report_mean <- 0:199 %*% ETH_report_delays
report_mean
report_var <- ETH_report_delays %*% ((0:199) - c(report_mean))^2
sqrt(report_var)
report_shape <- (report_mean^2)/report_var
report_scale <- report_var/report_mean

plot(ETH_incubation[1:30], type="l", ylim=c(0,0.2))
lines(ETH_report_delays)

lines(dgamma(0:199, shape = incubation_shape, scale = incubation_scale), lty=2, col="red")
lines(dgamma(0:199, shape = report_shape, scale = report_scale), lty=2, col="red")

lines(dlnorm(0:199, meanlog = convert_to_logmean(incubation_mean, sqrt(incubation_var)),
             sdlog = convert_to_logsd(incubation_mean, sqrt(incubation_var))),
      lty=2, col="blue")
lines(dlnorm(0:199, meanlog = convert_to_logmean(report_mean, sqrt(report_var)),
             sdlog = convert_to_logsd(report_mean, sqrt(report_var))),
      lty=2, col="blue")

incubation <- function(x) dgamma(x, shape = incubation_shape, scale = incubation_scale)
reporting <- function(y) dgamma(y, shape = report_shape, scale = report_scale)
# convolution integral
delay <- function(z) integrate(function(x,z) reporting(z-x)*incubation(x),-Inf,Inf,z)$value
delay <- Vectorize(delay)
z <- 0:199
lines(z,delay(z),lty=2,col="blue")

# parameters of resulting convolution
z <- 0:10e4
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
shape <- ((delay_mean)^2)/delay_var
scale <- delay_var/(delay_mean)

plot(rtlive_delays$p_delay[2:31], type="l", ylim=c(0,0.14))
lines(5:30, dgamma(0:25, shape = shape, scale = scale), lty=2, col="red")
lines(0:30, dlnorm(0:30, meanlog = convert_to_logmean(5, 1),
                   sdlog = convert_to_logsd(5, 1)), lty=2, col="blue")
lines(5:30, dlnorm(0:25, meanlog = convert_to_logmean(delay_mean, delay_sd),
                   sdlog = convert_to_logsd(delay_mean, delay_sd)), lty=2, col="blue")

incubation <- function(x) dlnorm(x, meanlog = convert_to_logmean(5,1),
                                 sdlog = convert_to_logsd(5,1))
reporting <- function(y) dlnorm(y, meanlog = convert_to_logmean(delay_mean, delay_sd),
                                sdlog = convert_to_logsd(delay_mean, delay_sd))
# convolution integral
delay <- function(z) integrate(function(x,z) reporting(z-x)*incubation(x),-Inf,Inf,z)$value
delay <- Vectorize(delay)
z <- 0:30
lines(z,delay(z),lty=2,col="darkgreen")

# parameters of resulting convolution
z <- 0:10e4
mean <- z %*% delay(z)
mean
sd <- sqrt(delay(z) %*% (z - c(mean))^2)
sd

################
# epiforecasts #
################
convert_to_mean <- function(logmean, logsd){
  return(exp(logmean) * sqrt(exp(logsd^2)))
}
convert_to_sd <- function(logmean, logsd){
  return(sqrt(exp(logsd^2)-1) * exp(logmean) * sqrt(exp(logsd^2)))
}

incubation_logmean <- 1.62
incubation_logsd <- 0.418
incubation_mean <- convert_to_mean(incubation_logmean, incubation_logsd)
incubation_var <- convert_to_sd(incubation_logmean, incubation_logsd)^2

report_delay_logmean <- 0.832
report_delay_logsd <- 1.44
report_delay_mean <- convert_to_mean(report_delay_logmean, report_delay_logsd)
report_delay_var <- convert_to_sd(report_delay_logmean, report_delay_logsd)^2

incubation <- function(x) dlnorm(x, meanlog = incubation_logmean, sdlog = incubation_logsd)
report_delay <- function(y) dlnorm(y, meanlog = report_delay_logmean, sdlog = report_delay_logsd)
# convolution integral
delay <- function(z) integrate(function(x,z) report_delay(z-x)*incubation(x),-Inf,Inf,z)$value
delay <- Vectorize(delay)
z <- 0:40
plot(z,delay(z), type="l")

# parameters of resulting convolution
z <- 0:10e4
mean <- z %*% delay(z)
mean
sd <- sqrt(delay(z) %*% (z - c(mean))^2)
sd


############################################
# plot distributions used for epiforecasts #
############################################

o <- 0.1 # replacement for zeros that would lead to mathematical issues
oo <- 0.5
m <- 2 # mean incubation period and reporting delay in case of constant delay distributions (shifted back manually after the estimation)
mean_incubation <-   c(5.3, m,  m,  m,  m,  m,  5.5,  5,   m)
sd_incubation <-     c(3.2, oo, oo, oo, oo, oo, 2.4,  oo,  oo)
mean_report_delay <- c(5.5, m,  m,  m,  m,  m,  6.5,  7.1, m)
sd_report_delay <-   c(3.8, o,  o,  o,  o,  o,  17.1, 5.9, o)
delay_shift <-       c(0,   3,  -3, -6, 4,  4,  0,    0,   4)

params <- data.frame(incubation_mean=mean_incubation, incubation_sd=sd_incubation,
                     report_delay_mean=mean_report_delay, report_delay_sd=sd_report_delay,
                     delay_shift=delay_shift)
methods <- c("ETH", "RKI", "Ilmenau", "SDSC", "Zi", "AGES", "epiforecasts", "rtlive", "globalrt")
rownames(params) <- methods

incubation <- function(x) dlnorm(x, meanlog = convert_to_logmean(params[method, "incubation_mean"],
                                                                 params[method, "incubation_sd"]),
                                 sdlog = convert_to_logsd(params[method, "incubation_mean"],
                                                          params[method, "incubation_sd"]))
reporting <- function(y) dlnorm(y, meanlog = convert_to_logmean(params[method, "report_delay_mean"],
                                                                params[method, "report_delay_sd"]),
                                sdlog = convert_to_logsd(params[method, "report_delay_mean"],
                                                         params[method, "report_delay_sd"]))
# convolution integral
delay <- function(z) integrate(function(x, z) reporting(z-x)*incubation(x),-Inf,Inf,z)$value
delay <- Vectorize(delay)
z <- 0:1000

par(mfrow=c(3,3))
for (method in methods){
  # calculate delay distributions as sum over two log-normal distributions
  d <- delay(z)
  
  # delays used originally in source of parameters
  if (method == "ETH"){
    d_org <- constant_delay_distributions$`Confirmed cases` # read at the top of this script
  } else d_org <- delays_ETH[[method]]$`Confirmed cases`
  plot(0:30, d_org[1:31], pch=4,
       xlim=c(0,30), ylim=c(0, max(d_org, d)),
       main=method, xlab="days", ylab="p")
  
  # resulting delays used for epiforecasts estimation after manual shift
  z_shifted <- z-params[method, "delay_shift"]
  lines(z_shifted[1:31], d[1:31], col="red")
  
  # delays used for epiforecasts estimation
  lines(z[1:31], d[1:31], col="blue")
  
  # add mean and sd as legend
  mean <- z_shifted %*% d
  sd <- sqrt(d %*% (z_shifted - c(mean))^2)
  legend(x="topright", legend=c(paste("mean", round(mean, digits=1)),
                                paste("sd     ", round(sd, digits=1))))
}
par(mfrow=c(1,1))



#################################################################
# Reporting delay RKI line list (Nowcast vs rtlive aggregation) #
#################################################################
raw_data <- read_csv("Rt_estimate_reconstruction/rtlive/rtlive-global/data/RKI_COVID19_21_11_23.csv")
raw_data$Meldedatum <- as_date(raw_data$Meldedatum)
raw_data$Refdatum <- as_date(raw_data$Refdatum)
data <- raw_data %>%
  dplyr::select(AnzahlFall, Meldedatum, Refdatum) %>%
  dplyr::filter(Meldedatum <= as_date("2021-07-10")) %>%
  group_by(Meldedatum) %>%
  summarise(AnzahlFall=sum(AnzahlFall), Refdatum=mean(Refdatum))
View(data)

(data$AnzahlFall %*% as.numeric(data$Meldedatum - data$Refdatum)) / sum(data$AnzahlFall)


