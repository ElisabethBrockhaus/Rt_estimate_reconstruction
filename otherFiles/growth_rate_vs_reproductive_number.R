library(ggplot2)

setwd("../..")
getwd()

source("Rt_estimate_reconstruction/prepared_plots.R")


# define growth rates
r_t <- seq(-0.5, 0.5, 10e-5)

# SIR case
calc_R_t <- function(w) {
  R <- 1:length(r_t)
  for (i in R) {
    R[i] <- 1 + (c(w) * r_t[i])
    if (R[i] <= 0) {
      R[i] <- NA
    }
  }
  return(R)
}
R_t_3 <- calc_R_t(3)
R_t_4 <- calc_R_t(4)
R_t_5 <- calc_R_t(5)
R_t_6 <- calc_R_t(6)
R_t_7 <- calc_R_t(7)

# upper bounds
calc_max_R_t <- function(w) exp(w * r_t)
R_t_3_max <- calc_max_R_t(3)
R_t_4_max <- calc_max_R_t(4)
R_t_5_max <- calc_max_R_t(5)
R_t_6_max <- calc_max_R_t(6)
R_t_7_max <- calc_max_R_t(7)

# lower bound
R_t_min <- 1:length(r_t)
for (i in R_t_min) {
  if (r_t[i] < 0) {
    R_t_min[i] <- 0
  } else {
    R_t_min[i] <- 1
  }
}
R_t_3_min <- R_t_4_min <- R_t_5_min <- R_t_6_min <- R_t_7_min <- R_t_min

# join R_t
data <- data.frame(r = r_t, R3 = R_t_3, R4 = R_t_4, R5 = R_t_5, R6 = R_t_6, R7 = R_t_7)
names(data) <- c("r", "3", "4", "5", "6", "7")
data_max <- data.frame(r = r_t, R3 = R_t_3_max, R4 = R_t_4_max, R5 = R_t_5_max, R6 = R_t_6_max, R7 = R_t_7_max)
names(data_max) <- c("r", "3", "4", "5", "6", "7")
data_min <- data.frame(r = r_t, R = R_t_min)
names(data_min) <- c("r", "all")

# reshape data
data_reshaped  <- data %>%
  gather("variable", "value", 2:dim(data)[2]) %>%
  mutate(type = "SEIR case")
data_max_reshaped  <- data_max %>%
  gather("variable", "value", 2:dim(data_max)[2]) %>%
  mutate(type = "upper bound")
data_min_reshaped  <- data_min %>%
  gather("variable", "value", 2:dim(data_min)[2]) %>%
  mutate(variable = "all") %>%
  mutate(type = "lower bound")

data_joined <- data_reshaped %>% rbind(data_max_reshaped) %>% rbind(data_min_reshaped)

col_values <- get_colors(c("3", "4", "5", "6", "7", "all"), palette = "Dark2", name_consensus = "all")
line_values <- c("SEIR case" = "solid", "upper bound" = "dashed", "lower bound" = "dotted")

# plot
R_plot <- ggplot(data = data_joined, aes(x = r, y = value,
                                         colour = variable, linetype = type,
                                         group=interaction(variable, type))) +
  labs(x = "growth rate r", y = "reproductive number R") +
  theme_minimal() +
  theme(
    axis.text=element_text(size=16),
    axis.title=element_text(size=18),
    legend.text=element_text(size=16),
    legend.title=element_text(size=18),
    axis.line = element_line(),
    axis.line.y.right = element_line(),
    axis.line.x.top = element_line(),
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent")
  ) +
  geom_line() +
  #geom_line(data = data_max_reshaped, aes(group = variable, color = variable), linetype = "dashed") +
  #geom_line(data = data_min_reshaped, aes(group = variable, color = variable), linetype = "dotted") +
  scale_color_manual(values = col_values, name = "mean generation time") +
  guides(colour = guide_legend(override.aes = list(linetype = c(rep("solid", 5), "dotted")))) +
  scale_linetype_manual(values = line_values, name = "") +
  #guides(colour = guide_legend("title"), values = c("upper bound", "SIR case", "lower bound"),
  #       size="legend", legend.position = c(-0.25, 8),
  #       shape = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dotted")))) +
  coord_cartesian(ylim = c(-0.2, 10.2), expand = F)

print(R_plot)
ggsave(R_plot, filename = "Figures/growthrate_vs_reproductivenumber.pdf",  bg = "transparent",
       width = 13.1, height = 5.8)



###############################################################################
# for presentation slides
R <- 1.07
incidence <- data.frame(date=1:100, I=c(rep(1,10), rep(0,90)))
for (t in 2:100) {
  incidence$I[t] <- incidence$I[t-1] * R
}
plot(incidence, type="l")

R <- 1.07
incidence <- data.frame(date=1:100, I=c(rep(1,10), rep(0,90)))
for (t in 3:100) {
  incidence$I[t] <- incidence$I[t-2] * R
}
plot(incidence, type="l")

#incidence$I <- incid$I[1:100]

incidence$gr <- c(NA, (incidence$I[2:100] - incidence$I[1:99]) / incidence$I[1:99])

incidence$R_gt3 <- 1 + (incidence$gr / 3)
incidence$R_gt4 <- 1 + (incidence$gr / 4)
incidence$R_gt5 <- 1 + (incidence$gr / 5)
incidence$R_gt6 <- 1 + (incidence$gr / 6)
incidence$R_gt7 <- 1 + (incidence$gr / 7)

gr <- 0.03
incidence <- data.frame(date=1:100, I=c(1, rep(0,99)))
for (t in 2:100) {
  incidence$I[t] <- incidence$I[t-1] * (1+gr)
}
plot(incidence, type="l", xlab="day", ylab="new infections")

for (gr in seq(0.04, 0.07, by=0.01)){
  incidence <- data.frame(date=1:100, I=c(1, rep(0,99)))
  for (t in 2:100) {
    incidence$I[t] <- incidence$I[t-1] * (1+gr)
  }
  lines(incidence)
}
#text(locator(), labels = c("GT = 3", "GT = 4", "GT = 5", "GT = 6", "GT = 7"))
text(x=c(35.5, 43, 52, 66, 90), y=rep(18,5), pos=4,
     labels = c("GT = 3", "GT = 4", "GT = 5", "GT = 6", "GT = 7"))
text(x=c(5), y=c(18), pos=4,
     labels = c("R = 1.01"))

incidence_data <- RKI_incid %>%
  inner_join(rtlive_incid, by = "date") %>%
  inner_join(ETH_incid_mean, by = "date") %>%
  inner_join(ilmenau_incid, by = "date") %>%
  inner_join(sdsc_incid, by = "date") %>%
  inner_join(epiforecasts_incid, by = "date") %>%
  inner_join(globalrt_incid, by = "date")

names(incidence_data) <- c("date", "RKI (nowcast)", "rtlive", "ETH (deconvoluted)", "Ilmenau", "SDSC (smoothed)", "epiforecasts", "globalrt") 

# reshape data
incid  <- incidence_data %>%
  gather("variable", "value", 2:dim(incidence_data)[2]) %>%
  filter(date >= as.Date("2021-02-01")) %>%
  filter(date <= as.Date("2021-05-31"))

# plot
R_plot <- ggplot(data = incid, aes(x = date, y = value)) +
  labs(x = NULL, y = "new infections") +
  scale_x_date(limits = as.Date(c(min(incid$date),max(incid$date))),
               breaks = seq(as.Date(min(incid$date)), as.Date(max(incid$date)), "month"),
               date_labels = "%B %d", expand = c(0,1)) +
  theme_minimal() +
  theme(
    axis.text=element_text(size=16),
    axis.title=element_text(size=18),
    legend.text=element_text(size=16),
    legend.title=element_text(size=18),
    axis.line = element_line(),
    axis.line.y.right = element_line(),
    axis.line.x.top = element_line(),
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent")
  ) +
  geom_line(aes(group = variable, color = variable)) +
  scale_colour_viridis_d(end = 0.9, name="group")

print(R_plot)
ggsave(R_plot, filename = "incidence_data.pdf",  bg = "transparent",
       width = 13.1, height = 6.3)