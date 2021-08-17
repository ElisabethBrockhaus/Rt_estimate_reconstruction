library(EpiEstim)

### estimate as in RKI2020
estimate_RKI_R <- function(incid){
  
  # calculate estimates
  estimate <- rep(NA, nrow(incid))
  for (t in 11:nrow(incid)){
    estimate[t-1] <- round(sum(incid$I[t-0:6]) / sum(incid$I[t-4:10]), digits = 2)
  }
  
  return(estimate)
}



### estimate as in RKI2020 but using EpiEstim (Bayes)
estimate_RKI_R_EpiEstim <- function(incid){
  
  # deterministic serial interval (4 days)
  serial_interval <- c(0, 0, 0, 0, 1)
  
  # start and end for estimations at each time point (window size = 7)
  start <- 5:(nrow(incid)-6)
  end <- start + 6
  
  # estimation with EpiEstim function
  r_EpiEstim <- estimate_R(incid$I,
                           method = "non_parametric_si",
                           config = make_config(list(t_start=start, t_end=end,
                                                     si_distr=serial_interval)))
  
  # assign estimates to time points properly
  estimate <- rep(NA, nrow(incid))
  estimate[10:(nrow(incid)-1)] <- round(r_EpiEstim$R$`Mean(R)`, digits = 2)
  
  return(estimate)
}



### estimate as in Huisman2021
estimate_ETH_R <- function(incid,
                           window=3,
                           delays=c(Cori = 10, WallingaTeunis = 5),
                           truncations=list(
                             left = c(Cori = 5, WallingaTeunis = 0),
                             right = c(Cori = 0, WallingaTeunis = 8))){
  
  # load functions for estimation
  source("ETH/covid-19-re-shiny-app/app/otherScripts/3_utils_doReEstimates.R")
  
  # Run EpiEstim
  estimateRaw <- doReEstimation(
    incid,
    slidingWindow = window,
    methods = c("Cori"),
    variationTypes = c("slidingWindow"),
    delays = delays,
    truncations = truncations
  )
  
  #estimate <- cleanCountryReEstimate(estimateRaw,
  #                                   method = 'bootstrap',
  #                                   rename_types = F) 
  #View(estimate)
  
  #estimate <- estimate %>%
  #  left_join(
  #    dplyr::select(popData, region, countryIso3),
  #    by = c("region")
  #  )
  estimate<-estimateRaw
  estimate <- estimate[estimate$variable=="R_mean",]
  
  estimate <- full_join(incid[,c("date", "deconvoluted")], estimate)$value
  
  return(estimate)
}



### estimate as in Hotz2020
estimate_Ilmenau_R <- function(incid){
  
  # infectivity profile
  infectivity <- c((0:3)/3, 1, (5:0)/5)
  names(infectivity) <- seq_along(infectivity)
  infectivity <- infectivity / sum(infectivity)
  
  # other parameters
  width <- 1
  report.delay <- 7
  alpha <- 0.05
  
  # calculate estimates
  estimate <- repronum(
    new.cases = incid$new.cases,
    profile = infectivity,
    window = width,
    delay = report.delay,
    conf.level = 1 - alpha,
    pad.zeros = TRUE
  )$repronum

  return(estimate)
}

### function from Ilmenau repo for Rt estimation
repronum <- function(
  new.cases, # I
  profile, # w
  window = 1, # H
  delay = 0, # Delta
  conf.level = 0.95, # 1-alpha
  pad.zeros = FALSE,
  min.denominator = 5,
  min.numerator = 5
) {
  # pad zeros if desired
  if(pad.zeros) new.cases <- c(rep(0, length(profile) - 1), new.cases)
  
  # compute convolutions over h, tau and both, respectively
  sum.h.I <- as.numeric(stats::filter(new.cases, rep(1, window),
                                      method = "convolution", sides = 1))
  sum.tau.wI <- as.numeric(stats::filter(new.cases, c(0, profile),
                                         method = "convolution", sides = 1))
  sum.htau.wI <- as.numeric(stats::filter(sum.tau.wI, rep(1, window),
                                          method = "convolution", sides = 1))
  
  # estimators
  repronum <- ifelse(sum.h.I < min.numerator, NA, sum.h.I) / ifelse(sum.htau.wI < min.denominator, NA, sum.htau.wI)
  
  # standard errors
  repronum.se <- sqrt(repronum / sum.htau.wI)
  
  # shift by delay
  repronum <- c(repronum, rep(NA, delay))[(1:length(repronum)) + delay]
  repronum.se <- c(repronum.se,
                   rep(NA, delay))[(1:length(repronum.se)) + delay]
  
  # standard normal qunatile
  q <- qnorm(1 - (1-conf.level) / 2)
  
  # return data.frame with as many rows as new.cases
  ret <- data.frame(
    repronum = repronum,
    repronum.se = repronum.se,
    ci.lower = repronum - q * repronum.se,
    ci.upper = repronum + q * repronum.se
  )
  if(pad.zeros) ret[-(1:(length(profile)-1)),] else ret
}


### estimate as in Richter2020
estimate_AGES_R <- function(incid, window = 13){
  
  # filter incid to dates where incidence data is available
  start_i <- min(which(!is.na(incid$I)))
  end_i <- max(which(!is.na(incid$I)))
  incid_wo_na <- incid[start_i:end_i,]
  
  # deterministic serial interval (4 days)
  serial_interval <- c(0, dgamma(1:16, shape = 2.88, scale = 1.55))
  
  # start and end for estimations at each time point
  start <- 2:(nrow(incid_wo_na) + 1 - window)
  end <- start - 1 + window
  
  # estimation with EpiEstim function
  r_EpiEstim <- estimate_R(incid_wo_na$I,
                           method = "non_parametric_si",
                           config = make_config(list(t_start=start, t_end=end,
                                                     si_distr=serial_interval)))
  
  # assign estimates to time points properly
  estimate <- rep(NA, nrow(incid))
  estimate[(start_i+window):end_i] <- r_EpiEstim$R$`Mean(R)`
  
  return(estimate)
}












