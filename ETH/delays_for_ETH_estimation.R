library(readr)

# more complex for ETH estimation
delay_RKI <- vector("list", 1)
names(delay_RKI) <- c("Confirmed cases")
delay_RKI["Confirmed cases"] <- list(c(rep(0, params["RKI", "delay"]),
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
delay_Zi <- vector("list", 1)
names(delay_Zi) <- c("Confirmed cases")
delay_Zi["Confirmed cases"] <- list(c(1, rep(0,199)))
delay_AGES <- delay_globalrt <- delay_Zi

delay_rtlive <- vector("list", 1)
names(delay_rtlive) <- c("Confirmed cases")
delay_rtlive["Confirmed cases"] <- list(c(read_csv("rtlive/p_delay.csv")$p_delay,
                                          rep(0, 200-66)))

# TODO: add correct delay distributions for epiforecasts
delays_ETH <- list(list(), delay_RKI, delay_Ilmenau, delay_SDSC, delay_Zi,
                   delay_AGES, delay_globalrt, delay_rtlive, list())
names(delays_ETH) <- methods