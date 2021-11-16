library(readr)

setwd("..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
# and for SDSC method "covid-19-forecast" (https://renkulab.io/gitlab/covid-19/covid-19-forecast/-/tree/master)
getwd()

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

# load incidence data used by research groups
rtlive_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid.csv")

#RKI_incid <- load_incidence_data(method = "RKI")

#ETH_incid <- load_incidence_data(method = "ETHZ_sliding_window")[, c("date", "value")]
#ETH_incid_mean <- ETH_incid$value %>% aggregate(by = list(ETH_incid$date), FUN = mean)
#names(ETH_incid_mean) <- c("date", "I")

ilmenau_incid <- load_incidence_data(method = "ilmenau")

#sdsc_incid <- load_incidence_data(method = "sdsc")

epiforecasts_incid <- load_incidence_data(method = "epiforecasts")

globalrt_data <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/time_series_covid19_confirmed_global.csv")
dates_globalrt <- as.Date(names(globalrt_data[,6:595]), format = "%m/%d/%y")
globalrt_total <- data.frame(date = dates_globalrt, cases = t(globalrt_data[1,6:595]))
globalrt_incid <- data.frame(date = dates_globalrt,
                             I =t(globalrt_data[1,6:595]) - dplyr::lag(t(globalrt_data[1,6:595]), n = 1))

# join incidence time series
incidence_data <- rtlive_incid %>%
  #inner_join(RKI_incid, by = "date") %>%
  #inner_join(ETH_incid_mean, by = "date") %>%
  inner_join(ilmenau_incid, by = "date") %>%
  #inner_join(sdsc_incid, by = "date") %>%
  inner_join(epiforecasts_incid, by = "date") %>%
  inner_join(globalrt_incid, by = "date")
names(incidence_data) <- c("date", "rtlive", #"RKI (nowcast)", "ETH (deconvoluted)",
                           "Ilmenau", #"SDSC (smoothed)",
                           "epiforecasts", "globalrt") 

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
