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