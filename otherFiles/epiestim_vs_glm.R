# Johannes Bracher, johannes.bracher@kit.edu

library(EpiEstim)
library(MASS)
library(zoo)
library(dplyr)

setwd("../..")
getwd()

source("Rt_estimate_reconstruction/calculate_estimates.R")
source("Rt_estimate_reconstruction/prepared_plots.R")


# read in incidence data:
dat <- read.csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_07_10.csv", colClasses = c("date" = "Date"))

# function to estimate R_t using glm or glm.nb:
# Arguments:
# incid: the incidence time series
# t: the time indices for which to compute R_t
# si_distr: the serial interval distribution; as in EpiEstim starting with a value for 0 days
# window: the window size (estimates are assigned to the last day of the window)
# family: the type of distribution; either "Poisson" or "NegBin"
estimate_R_glm <- function(incid, t, si_distr, window = 7, family = "Poisson"){
  si_distr <- si_distr[-1] # remove entry for zero days
  max_si <- length(si_distr)
  
  # add zeros at the beginning of the time series
  incid <- c(rep(0, max_si), incid)
  
  # compute weighted sums of lagged incidence values:
  lagged0 <- rollapply(incid, function(x) sum(rev(si_distr)*x), width = length(si_distr))
  lagged <- c(rep(NA, length(si_distr)), lagged0[-length(lagged0)])
  
  # vectors to store results:
  R <- L <- U <- rep(NA, length(t))
  
  # initialize progress bar:
  pb = txtProgressBar(min = 0, max = length(t), initial = 0) 
  
  for(i in seq_along(t)){
    # subset to observations inside window
    subs <- t[i] + max_si - (window - 1):0
    obs_temp <- incid[subs]
    lagged_temp <- lagged[subs]
    
    # fit GLM:
    if(family == "Pois"){
      fit <- glm(obs_temp ~ -1 + lagged_temp, family = poisson(link = "identity"))
    }
    if(family == "NegBin"){
      fit <- glm.nb(obs_temp ~ -1 + lagged_temp, link = identity)
    }
    
    # extract coefficients and CIs:
    R[i] <- fit$coefficients[1]
    ci <- suppressMessages(confint(fit))
    L[i] <- ci[1]
    U[i] <- ci[2]
    
    # update progress bar:
    setTxtProgressBar(pb,i)
  }
  
  # format results:
  result <- data.frame(t = t, R = R, R_lower = L, R_upper = U)
  
  return(result)
}

# define serial interval distribution
gtd <- get_infectivity_profile("gamma", 4, 4)

## estimate using EpiEstim                     
R_EpiEstim <- estimate_R(incid = dat, 
                         method = "non_parametric_si",
                         config = make_config(list(si_distr = gtd)))$R
R_EpiEstim$date <- tail(dat$date, nrow(R_EpiEstim))

# estimate using Poisson GLM
timepoints <- 8:nrow(dat)
R_glm_pois <- estimate_R_glm(dat$I,
                             si_distr = gtd,
                             t = timepoints,
                             family = "Pois", window = 7)
R_glm_pois$date <- dat$date[timepoints]


# estimate using NegBin GLM
R_glm_NB <- estimate_R_glm(dat$I,
                           si_distr = gtd,
                           t = timepoints,
                           family = "NegBin", window = 7)
R_glm_NB$date <- dat$date[timepoints]

# aggregate estimates
estimates <- R_EpiEstim %>%
  dplyr::select("date", "Mean(R)", "Quantile.0.025(R)", "Quantile.0.975(R)") %>%
  rename(R_EpiEstim = `Mean(R)`, lower_EpiEstim = `Quantile.0.025(R)`, upper_EpiEstim = `Quantile.0.975(R)`) %>%
  # full_join(R_glm_pois %>%
  #             dplyr::select(date, R, R_lower, R_upper) %>%
  #             rename(R_GLM_Poisson = R, lower_GLM_Poisson = R_lower, upper_GLM_Poisson = R_upper),
  #           by = "date") %>%
  full_join(R_glm_NB %>%
              dplyr::select(date, R, R_lower, R_upper) %>%
              rename(R_GLM_NegBin = R, lower_GLM_NegBin = R_lower, upper_GLM_NegBin = R_upper),
            by = "date")

# find ylim
colMin <- function(data) sapply(data, min, na.rm = TRUE)
colMax <- function(data) sapply(data, max, na.rm = TRUE)
ylim_input_l <- min(colMin(estimates %>% 
                             dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>%
                             dplyr::select(starts_with("lower_"))))
ylim_input_u <- max(colMax(estimates %>%
                             dplyr::filter(date>="2021-01-01", date<="2021-06-10") %>%
                             dplyr::select(starts_with("upper_"))))

# Plot point estimates with uncertainty intervals
plot_for_comparison(estimates,
                    comp_methods = c("EpiEstim/Poisson GLM", "NegBin GLM"),
                    col_palette = "Dark2", name_consensus = "none",
                    legend_name = "estimation method", filenames = "_influence_type_of_distribution.pdf",
                    sort_numerically = FALSE, point_and_CI=T,
                    ylim_l = ylim_input_l, ylim_u = ylim_input_u)
# The point estimates agree almost perfectly. The uncertainty intervals are quite different.


