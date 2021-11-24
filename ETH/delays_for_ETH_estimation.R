library(readr)

# more complex for ETH estimation
delay_RKI <- vector("list", 1)
names(delay_RKI) <- c("Confirmed cases")
delay_RKI["Confirmed cases"] <- list(c(rep(0, params["RKI", "delay"]+3),
                                       1,
                                       rep(0, 199-params["RKI", "delay"])))
delay_Ilmenau <- vector("list", 1)
names(delay_Ilmenau) <- c("Confirmed cases")
delay_Ilmenau["Confirmed cases"] <- list(c(rep(0, params["Ilmenau", "delay"]),
                                           1,
                                           rep(0, 199-params["Ilmenau", "delay"])))
delay_SDSC <- vector("list", 1)
names(delay_SDSC) <- c("Confirmed cases")
delay_SDSC["Confirmed cases"] <- list(c(rep(0, params["SDSC", "delay"]),
                                        1,
                                        rep(0, 199-params["SDSC", "delay"])))

delay_AGES <- vector("list", 1)
names(delay_AGES) <- c("Confirmed cases")
delay_AGES["Confirmed cases"] <- list(c(1, rep(0,199)))
delay_globalrt <- delay_AGES

delay_rtlive <- vector("list", 1)
names(delay_rtlive) <- c("Confirmed cases")
delay_rtlive["Confirmed cases"] <- list(c(read_csv("rtlive/p_delay.csv")$p_delay,
                                          rep(0, 200-66)))

# delays used in epiforecasts estimation
incubation_logmean <- 1.62
incubation_logsd <- 0.418
report_delay_logmean <- 0.832
report_delay_logsd <- 1.44

incubation <- function(x) dlnorm(x, meanlog = incubation_logmean, sdlog = incubation_logsd)
report_delay <- function(y) dlnorm(y, meanlog = report_delay_logmean, sdlog = report_delay_logsd)
# convolution integral
delay <- function(z) integrate(function(x,z) report_delay(z-x)*incubation(x),-Inf,Inf,z)$value
delay <- Vectorize(delay)

delay_epiforecasts <- vector("list", 1)
names(delay_epiforecasts) <- c("Confirmed cases")
delay_epiforecasts["Confirmed cases"] <- list(delay(0:199))

delays_ETH <- list(list(), delay_RKI, delay_Ilmenau, delay_SDSC, delay_AGES,
                   delay_epiforecasts, delay_rtlive, delay_globalrt)
names(delays_ETH) <- methods
