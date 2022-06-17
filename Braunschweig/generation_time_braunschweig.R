# Reconstructing the generation time distribution of the HZI model from Khailaie and Mitra et al (2020)
# https://doi.org/10.1186/s12916-020-01884-4

# Johannes Bracher, johannes.bracher@kit.edu

# General strategy:
# I sample the paths of infected individuals through the compartmental system with exponentially distributed times
# of stay in each compartment the individual reaches. In each infectious compartment, infection times are sampled
# according to a Poisson process given the time of stay.
# Then all sampled times until secondary infection from 5000 individuals are pooled to obtain the generation
# time distribution.


# set seed:
set.seed(123)

# read in sampled parameter values from repo:
param_vals <- read.csv("https://gitlab.com/simm/covid19/secir/-/raw/master/codes/settings/param_random.csv?inline=false")
colMeans(param_vals)

# set parameter values
# are all in agreement with file in repo:
(R2 <- 1/3.2)
(R3 <- 1/(5.2 - 1/R2))
(R4 <- 1/7)
(R6 <- 1/(4.25 - 3.7)) # 0.55 # corresponds to R6prime
(R9 <- 1/(1/R3 + 0.5/R4))
(R11 <- 1/3.7)
(R12 <- 1/(1/R4 - 1/R11))
alpha <- 0.22
beta <- 0.15
gamma <- 1
omega <- 1
chi <- 1
mu <- 0.085 # taken from repo
rho <- 0.13 # taken from repo

# a constant steering how any infections are generated for each individual
const <- 10

# function to sample infection times for one individual
# attention: this uses the above parameters from the global workspace rather
# than getting them as parameters.
sample_infection_times <- function(plot = FALSE){
  time_CI <- time_CR <- time_IX <- time_IR <- time_I <- time_IH <-  # times spent in different compartments
    # becomes_hospitalized <- becomes_detected <- # indicators: does the person get hospitalized / detected?
    n_infected_CI <- n_infected_I <- n_infected_IH <- n_infected_IR <- # number of infections generated during stay in different compartments
    n_infected_IX <- n_infected_CR <- 0
  
  # times at which the individual arrives in different compartments
  # remains NA if compartment is never entered
  start_CI <- start_CR <- start_IX <- start_IR <- start_I <- start_IH <-
    end_CI <- end_CR <- end_IX <- end_IR <- end_I <- end_IH <- NA
  # time points of the infections occurring during stays in the different compartments:
  infection_times_CI <- infection_times_CR <- infection_times_IX <-
    infection_times_IR <- infection_times_I <- infection_times_IH <- numeric(0)
  
  # sample time spent in compartment E:
  time_E <- rexp(1, rate = R2)
  
  # sample whether the individual becomes a carrier:
  becomes_carrier <- rbinom(1, 1, 1 - alpha)
  
  # case 1: becomes a carrier; CI
  if(becomes_carrier){
    # sample time spent in CI
    time_CI <- rexp(1, rate = R3)
    # sample number of infections caused during CI
    # note: const enters here to scale up/down transmission
    n_infected_CI <- rpois(1, const*gamma*time_CI)
    # starting time of CI is end tie of E
    start_CI <- time_E
    # compute end time of stay in CI:
    end_CI <- start_CI + time_CI
    # sample the infection times from a uniform over the already sampled time of stay
    infection_times_CI <- runif(n_infected_CI, min = start_CI, max = end_CI)
    
    # sample whether gets detected
    becomes_detected <- rbinom(1, 1, mu)
    # If detected (I):
    if(becomes_detected){ # I
      # sample time spent in I
      time_I <- rexp(1, R11)
      # sample number of infected generated during I:
      n_infected_I <- rpois(1, const*omega*time_I)
      # start time of I is end time of CI
      start_I <- end_CI
      # compute end time of I
      end_I <- start_I + time_I
      # sample infection time points
      infection_times_I <- runif(n_infected_I, min = start_I, max = end_I)
      
      # sample whether gets hospitalized (moves to IH):
      becomes_hospitalized <- rbinom(1, 1, rho)
      if(becomes_hospitalized){ # IH
        # sample time spent in IH
        time_IH <- rexp(1, R6)
        # sample number of infections in IH
        n_infected_IH <- rpois(1, const*beta*time_IH)
        # start time corresponds to end of IH
        start_IH <- end_I
        # compute end time
        end_IH  <- start_IH + time_IH
        # sample infection times:
        infection_times_IH <- runif(n_infected_IH, min = start_IH, max = end_IH)
        
      }else{ # if not hospitalized: # I_R
        # sample time spent in IR
        time_IR <- rexp(1, R12)
        # sample number of infections in IR
        n_infected_IR <- rpois(1, const*beta*time_IR)
        # start time corresponds to end of I
        start_IR <- end_I
        # compute end time
        end_IR <- start_IR + time_IR
        # smple infection times:
        infection_times_IR <- runif(n_infected_IR, min = start_IR, max = end_IR)
      }
    }else{# infected undetected IX
      # sample time spent in IX
      time_IX <- rexp(1, R4)
      # sample number of infections in IX
      n_infected_IX <- rpois(1, const*chi*time_IX)
      # start time is end of CI
      start_IX <- end_CI
      # compute end time:
      end_IX <- start_IX + time_IX
      # sample infection times:
      infection_times_IX <- runif(n_infected_IX, min = start_IX, max = end_IX)
    }
  }else{ # case 2: recovering directly; CR
    # sample time spent in CR
    time_CR <- rexp(1, rate = R9)
    # numer of infections caused in CR
    n_infected_CR <- rpois(1, const*gamma*time_CR)
    # start time is end of E
    start_CR <- time_E
    # compute end time
    end_CR <- start_CR + time_CR
    # sample infection times:
    infection_times_CR <- runif(n_infected_CR, min = start_CR, max = end_CR)
  }
  
  # make a little plot if desired
  if(plot){
    par(las = 1)
    plot(NULL, xlim = c(0, 50), ylim = c(0, 8), axes = FALSE, xlab = "time", ylab = "")
    axis(1)
    mtext(2, text = c("E", "CR", "CI", "IX", "I", "IR", "IH"), at = 1:7)
    lines(c(0, time_E), rep(1, 2), col = "black")
    
    lines(c(start_CR, end_CR), rep(2, 2), col = "black")
    points(infection_times_CR, rep(2, n_infected_CR), col = "red", cex = 0.5)
    
    lines(c(start_CI, end_CI), rep(3, 2), col = "black")
    points(infection_times_CI, rep(3, n_infected_CI), col = "red", cex = 0.5)
    
    lines(c(start_IX, end_IX), rep(4, 2), col = "black")
    points(infection_times_IX, rep(4, n_infected_IX), col = "red", cex = 0.5)
    
    lines(c(start_I, end_I), rep(5, 2), col = "black")
    points(infection_times_I, rep(5, n_infected_I), col = "red", cex = 0.5)
    
    lines(c(start_IR, end_IR), rep(6, 2), col = "black")
    points(infection_times_IR, rep(6, n_infected_IR), col = "red", cex = 0.5)
    
    lines(c(start_IH, end_IH), rep(7, 2), col = "black")
    points(infection_times_IH, rep(7, n_infected_IH), col = "red", cex = 0.5)
  }
  
  # pool all infection times:
  infection_times <- c(infection_times_CI, infection_times_CR, infection_times_I,
                       infection_times_IH, infection_times_IR, infection_times_IX)
  
  return(infection_times)
}

# can sample a few times to get an intuition:
sample_infection_times(plot = TRUE)
# usually only few short intervals, but sometimes really many really long ones.

# sample infection times for 1000 individuals:
n_sim <- 5000
all_infection_times <- sample_infection_times()
for(i in 1:n_sim){
  new_infection_times <- sample_infection_times()
  all_infection_times <- c(all_infection_times, new_infection_times)
}

# make a histogram and compute moments:
hist(all_infection_times, breaks = 0:100)
mean(all_infection_times)
sd(all_infection_times)
median(all_infection_times)

# a function to map reproductive numbers (R) to growth rates (r)
# this is based on formula (3.6) from Wallinga and Lipsitch (2007)
R_from_r <- function(r, a, y){
  # a are the histogram bounds
  # y are relative frequencies
  denominator <- sum(y*(
    exp(-r*a[-length(a)]) - exp(-r*a[-1])/(a[-1] - a[-length(a)])
  ))
  r/denominator
}

# define a grid of growth rate values:
values_r <- seq(-1, 1, by = 0.001)

# compute the HZI R values corresponding to the grid of r values (growth rates)
# a are the histogram bounds
a_hzi <- 0:ceiling(max(all_infection_times))
# histogram values / relative frequencies:
hi_hzi <- hist(all_infection_times, breaks = a_hzi)
y_hzi <- hi_hzi$density
# compute R values corresponding to the values in the grid of r values (growth rates)
values_R_hzi <- sapply(X = values_r,FUN =  R_from_r, y = y_hzi, a = a_hzi)
# need to fill in value for 0 manually (probably dividing by 0 somewhere)
values_R_hzi[values_r == 0] <- 1

# compute the consensus R values corresponding to the grid of r values (growth rates)
a_consensus <- a_hzi
# sample from generation time of consensus model:
gen_times_consensus <- rexp(100000, 1/4)
# histogram values / relative frequencies:
hi_consensus <- hist(gen_times_consensus, breaks = a_consensus)
y_consensus <- hi_consensus$density
# compute R values corresponding to the values in the grid of r values (growth rates)
values_R_consensus <- sapply(values_r, R_from_r, y = y_consensus, a = a_consensus)
# need to fill in value for 0 manually (probably dividing by 0 somewhere)
values_R_consensus[values_r == 0] <- 1

# compare the mappings of growth rate to R under the two specifications:
plot(values_r, values_R_hzi, type = "l", xlim = c(-1, 1), ylim = c(0, 10))
lines(values_r, values_R_consensus, col = "red")

# mapping of R values form HZI to R values in consensus setting:
plot(values_R_hzi, values_R_consensus, xlim = c(0, 3), ylim = c(0, 3), type = "l", col = "blue")
abline(h = 1, col = "lightgrey")
abline(v = 1, col = "lightgrey")
abline(0:1)


# write out results:
to_store <- data.frame(values_r, values_R_consensus, values_R_hzi)
write.csv(to_store, "mapping_hzi_to_consensus.csv", row.names = FALSE)
