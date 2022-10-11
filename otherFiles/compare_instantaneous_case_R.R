library(EpiEstim)

data(Flu2009)

mean <- 2.6
std <- 1.5

R_instantaneaous <- estimate_R(Flu2009$incidence, 
                               method="parametric_si",
                               config = make_config(list(
                                 t_start = seq(2, 26), t_end = seq(8, 32),
                                 mean_si = mean, std_si = std)))

R_case <- wallinga_teunis(Flu2009$incidence, 
                          method="parametric_si",
                          config = list(t_start = seq(2, 26), t_end = seq(8, 32),
                                        mean_si = mean, std_si = std,
                                        n_sim = 100))

plot(R_instantaneaous$R$t_start, R_instantaneaous$R$`Mean(R)`, type="l",
     xlab = "t_start", ylab = "R", main = "Instantaneous vs. Case Reproductive Number")
lines(R_case$R$t_start, R_case$R$`Mean(R)`, col="red")
lines(R_case$R$t_start + mean, R_case$R$`Mean(R)`, col="red", lty = 2)
legend(x="topright", legend=c("instantaneous (Cori)", "case (Wallinga-Teunis)", "case shifted by mean si"),
       col=c("black", "red", "red"), lty = c(1,1,2))

