library(EpiEstim)

### estimate as in RKI2020
estimate_RKI_R <- function(incid, window=7, gt_mean=4, gt_sd=0){
  
  estimate <- data.frame(date=incid$date, R_calc=rep(NA, nrow(incid)))
  
  if (gt_sd==0){
    # calculate estimates as in RKI2020
    for (t in (gt_mean+window):nrow(incid)){
      estimate[t-1, "R_calc"] <- round(sum(incid$I[t-0:(window-1)])
                                       / sum(incid$I[t-gt_mean:(gt_mean+window-1)]),
                                       digits = 2)
    }
  } else {
    # use package EpiEstim to deal with serial interval distribution
    
    # start and end for estimations at each time point (window size = 7)
    start <- 2:(nrow(incid) + 1 - window)
    end <- start - 1 + window
    
    # estimation with EpiEstim function
    r_EpiEstim <- estimate_R(incid$I,
                             method = "parametric_si",
                             config = make_config(list(t_start=start, t_end=end,
                                                       mean_si=gt_mean, std_si=gt_sd)))
    len_est <- length(r_EpiEstim$R$`Mean(R)`)
    
    # assign estimates to time points properly
    estimate[window+0:(len_est-1), "R_calc"] <- round(r_EpiEstim$R$`Mean(R)`, digits = 2)
  }
  
  return(estimate)
}



### estimate as in Huisman2021
estimate_ETH_R <- function(incid, window=3){
  
  region <- unique(incid$region)
  
  parameter_list <- get_parameters_ETH(incid, region)
  
  # load functions for estimation
  source("Rt_estimate_reconstruction/ETH/otherScripts/3_utils_doReEstimates.R")
  
  # Run EpiEstim
  countryEstimatesRaw <- doAllReEstimations(
    incid,
    slidingWindow = window,
    methods = c("Cori"),
    variationTypes = c("slidingWindow"),
    all_delays =  parameter_list[["all_delays"]],
    truncations =  parameter_list[["truncations"]],
    interval_ends =  parameter_list[["interval_ends"]],
    swissRegions = parameter_list[["swissRegions"]]
  )
  
  gc()
  cat("raw estimates done \n")
  estimatePath <- "Rt_estimate_reconstruction/ETH/data/temp"
  if (!dir.exists(estimatePath)) {
    dir.create(estimatePath)
  }
  qs::qsave(countryEstimatesRaw, file = str_c(estimatePath, "/countryEstimatesRaw.qs"))
  gc()
  
  countryEstimates <- cleanCountryReEstimate(countryEstimatesRaw, method = 'bootstrap')
  
  # add extra truncation of 4 days for all Swiss cantonal estimates due to consolidation
  if (region %in% c("CHE")) {
    days_truncated <- 4
    canton_list <- c("AG", "BE", "BL", "BS", "FR", "GE", "GR", "JU", "LU", "NE", "SG", "SO", "SZ", "TG", "TI",
                     "VD", "VS", "ZG", "ZH", "SH", "AR", "GL", "NW", "OW", "UR", "AI")
    
    countryEstimates_cantons <- countryEstimates %>%
      filter(region %in% canton_list) %>%
      group_by(country, region, source, data_type, estimate_type) %>%
      filter(row_number() <= (n() - days_truncated)) %>%
      ungroup()
    
    countryEstimates_CH <- countryEstimates %>%
      filter(!(region %in% canton_list))
    
    countryEstimates <- bind_rows(countryEstimates_cantons, countryEstimates_CH)
    
  }
  
  countryDataPath <- file.path("Rt_estimate_reconstruction/ETH/data/countryData",
                               str_c(region, "-Estimates.rds"))
  saveRDS(countryEstimates, file = countryDataPath)
  
  estimate <- full_join(unique(incid[, c("date")]),
                        countryEstimates,
                        by="date")[, c("date", "median_R_mean")]
  names(estimate) <- c("date", "R_calc")
  
  return(estimate)
}


get_parameters_ETH <- function(deconvolvedCountryData, region){
  
  stringencyDataPath <- file.path("Rt_estimate_reconstruction/ETH/data/countryData",
                                  str_c(region, "-OxCGRT.rds"))
  stringencyIndex <- readRDS(stringencyDataPath)
  
  swissRegions <- deconvolvedCountryData %>%
    filter(country %in% c("Switzerland", "Liechtenstein")) %>%
    dplyr::select(region) %>%
    distinct() %>%
    .$region
  
  # in Re estimation, the interval starts on interval_end + 1
  # so the intervention start dates need to be shifted to - 1
  interval_ends_df <- stringencyIndex %>%
    filter(c(0, diff(stringencyIndex$value, 1, 1)) != 0) %>%
    mutate(interval_ends = date - 1) %>%
    dplyr::select(region, interval_ends)
  
  interval_ends <- split(interval_ends_df$interval_ends, interval_ends_df$region)
  interval_ends[["default"]] <- interval_ends[[region]]
  
  ### add additional interval 7 days before last Re estimate, for country level only
  # discarding interval ends more recent than that.
  lastIntervalEnd <- deconvolvedCountryData %>%
    filter(data_type == "infection_Confirmed cases", region == region, replicate == 0) %>%
    slice_max(date) %>%
    pull(date)
  lastIntervalStart <- lastIntervalEnd - 7
  
  if (length(interval_ends) == 0) {
    interval_ends[[region]] <- lastIntervalStart
  } else {
    interval_ends[[region]] <- c(
      interval_ends[[region]][interval_ends[[region]] < lastIntervalStart],
      lastIntervalStart)
  }
  
  if (region == "CHE") {
    additionalRegions <- setdiff(swissRegions, names(interval_ends))
    for (iregion in additionalRegions) {
      interval_ends[[iregion]] <- interval_ends[[region]]
    }
  }
  
  ### Delays applied
  all_delays <- list(
    "infection_Confirmed cases" = c(Cori = 0, WallingaTeunis = -5),
    "infection_Confirmed cases / tests" = c(Cori = 0, WallingaTeunis = -5),
    "infection_Deaths" = c(Cori = 0, WallingaTeunis = -5),
    "infection_Hospitalized patients" = c(Cori = 0, WallingaTeunis = -5),
    "Confirmed cases" = c(Cori = 10, WallingaTeunis = 5),
    "Confirmed cases / tests" = c(Cori = 10, WallingaTeunis = 5),
    "Deaths" = c(Cori = 20, WallingaTeunis = 15),
    "Hospitalized patients" = c(Cori = 8, WallingaTeunis = 3),
    "infection_Excess deaths" = c(Cori = 0, WallingaTeunis = -5),
    "Excess deaths" = c(Cori = 20, WallingaTeunis = 15))
  
  truncations <- list(
    left = c(Cori = 5, WallingaTeunis = 0),
    right = c(Cori = 0, WallingaTeunis = 8))
  
  parameter_list <- list(swissRegions, interval_ends, all_delays, truncations)
  names(parameter_list) <- c("swissRegions", "interval_ends", "all_delays", "truncations")
  
  return(parameter_list)
}



### estimate as in Hotz2020
estimate_Ilmenau_R <- function(incid, window = 1,
                               gt_type=c("org", "gamma"),
                               gt_mean=5.61, gt_sd=4.24){
  
  # infectivity profile
  if (gt_type == "org"){
    infectivity <- c((0:3)/3, 1, (5:0)/5)
    names(infectivity) <- seq_along(infectivity)
    infectivity <- infectivity / sum(infectivity)
  } else {
    infectivity <- dgamma(1:11, shape = (gt_mean^2)/gt_sd, scale = gt_sd/gt_mean)
    infectivity <- infectivity / sum(infectivity)
  }

  # other parameters
  report.delay <- 7
  alpha <- 0.05
  
  # calculate estimates
  estimate <- repronum(
    new.cases = incid$I,
    profile = infectivity,
    window = window,
    delay = report.delay,
    conf.level = 1 - alpha,
    pad.zeros = TRUE
  )$repronum
  
  estimate <- data.frame(date=incid$date, R_calc=estimate)

  return(estimate)
}


### function from Ilmenau repo for their Rt estimation
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
  
  # standard normal quantile
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
estimate_AGES_R <- function(incid, window = 13, mean_si = 4.46, std_si = 2.63){
  
  # filter incid to dates where incidence data is available
  start_i <- min(which(!is.na(incid$I)))
  end_i <- max(which(!is.na(incid$I)))
  incid_wo_na <- incid[start_i:end_i,]
  
  # deterministic serial interval
  #serial_interval <- c(0, dgamma(1:16, shape = 2.88, scale = 1.55))
  
  # start and end for estimations at each time point
  start <- 2:(nrow(incid_wo_na) + 1 - window)
  end <- start - 1 + window
  
  # estimation with EpiEstim function
  r_EpiEstim <- estimate_R(incid_wo_na$I,
                           method = "parametric_si",
                           #method = "non_parametric_si",
                           config = make_config(list(t_start=start, t_end=end,
                                                     mean_si=mean_si, std_si=std_si)))
                                                     #si_distr=serial_interval)))
  len_est <- length(r_EpiEstim$R$`Mean(R)`)
  
  # assign estimates to time points properly
  estimate <- data.frame(date=incid$date, R_calc=c(rep(NA, (start_i+window-1)),
                                                   r_EpiEstim$R$`Mean(R)`))
  return(estimate)
}












