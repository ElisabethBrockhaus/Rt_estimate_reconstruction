library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# set system locale time to English for correct labelling of x axes
Sys.setlocale("LC_TIME", "English")

setwd("../..")
# needs to be the directory with the repos "Rt_estimate_reconstruction"
getwd()


######################
# get incidence data #
######################

# read RKI Nowcast data for RKI estimation
RKI_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/RKI_nowcast_21_07_10.csv")

# read smoothed RKI incidence data for SDSC estimation
SDSC_incid <-  read_csv("Rt_estimate_reconstruction/incidence_data/SDSC_incid_21_07_10.csv")

# read deconvolved RKI incidence data for ETH estimation
incid_for_ETH <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_for_ETH_21_07_10.csv")

# read incidence data used by rtlive (sourced from RKI line-list data) for other estimations
incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_07_10.csv")

# compare incidence time series used by RKI vs. SDSC vs. rtlive vs. ETH
incid_ETH_repl0 <- incid_for_ETH %>%
  dplyr::filter(replicate == 0) %>%
  dplyr::select(date, value) %>%
  rename("I" = "value")
plot(RKI_incid, type="l")
lines(incid, col="blue")
lines(SDSC_incid, col="red")
lines(incid_ETH_repl0, col="green")

# adjust by mean delay considered
RKI_incid_d0 <- RKI_incid %>%
  mutate(date = date + 3)
incid_ETH_repl0_d0 <- incid_ETH_repl0 %>%
  mutate(date = date + 11)

# join preprocessed incidence time series
incidence_data <- RKI_incid_d0 %>%
  inner_join(incid, by = "date") %>%
  inner_join(SDSC_incid, by = "date") %>%
  inner_join(incid_ETH_repl0_d0, by = "date")
names(incidence_data) <- c("date", "nowcast (RKI)", "none", "smoothing (SDSC)", "deconvolution (ETH)")

# reshape data
incidence  <- incidence_data %>%
  gather("variable", "value", 2:dim(incidence_data)[2]) %>%
  filter(date >= as.Date("2021-01-01")) %>%
  filter(date <= as.Date("2021-06-10"))

# plot
R_plot <- ggplot(data = incidence, aes(x = date, y = value)) +
  labs(x = NULL, y = "new infections") +
  scale_x_date(limits = as.Date(c(min(incidence$date),max(incidence$date))),
               breaks = seq(as.Date(min(incidence$date)), as.Date(max(incidence$date)), "month"),
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
  scale_color_manual(values=brewer.pal(name="Dark2",n=6)[3:6], name="data source")

print(R_plot)
ggsave(R_plot, filename = "Figures/preprocessed_incidence_data.pdf",  bg = "transparent",
       width = 13.1, height = 6.3)
