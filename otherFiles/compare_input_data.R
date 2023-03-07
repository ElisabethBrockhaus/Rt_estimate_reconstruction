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

#source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

# load incidence data used by research groups

########################
# Ilmenau/rtlive (RKI) #
########################
#rki_data <- read_csv("Rt_estimate_reconstruction/rtlive/rtlive-global/data/rtlive_data_21_11_23.csv")
#rki_incid <- rki_data %>%
#  dplyr::filter(region == "all") %>%
#  dplyr::select(date, new_cases) %>%
#  rename(c("I" = "new_cases"))
#write_csv(rki_incid, "Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_11_23.csv")
rki_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_11_23.csv")

rki_incid_nowcast <- read_csv("https://raw.githubusercontent.com/robert-koch-institut/SARS-CoV-2-Nowcasting_und_-R-Schaetzung/main/Archiv/Nowcast_R_2021-11-23.csv") %>%
  dplyr::select(Datum, PS_COVID_Faelle) %>%
  rename(date = Datum)

######################
# epiforecasts (WHO) #
######################
#who_incid <- load_incidence_data(method = "epiforecasts")
#write_csv(who_incid, "Rt_estimate_reconstruction/incidence_data/epiforecasts_incid_21_11_23.csv")
who_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/epiforecasts_incid_21_11_23.csv")

#######################
# globalrt/SDSC (JHU) #
#######################
#jhu_incid_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
#jhu_incid <- jhu_incid_global %>%
#  filter(`Country/Region` == "Germany") %>%
#  unlist(., use.names=TRUE)
#jhu_incid <- as.numeric(jhu_incid[5:dim(jhu_incid_global)[2]])
#jhu_incid <- diff(jhu_incid, lag=1)
#jhu_incid <- data.frame("date"=as.Date(names(jhu_incid_global)[6:dim(jhu_incid_global)[2]], format = "%m/%d/%y"),
#                        "I"=jhu_incid)
#write_csv(jhu_incid, "Rt_estimate_reconstruction/incidence_data/jhu_incid_21_11_23.csv")
jhu_incid <- read_csv("Rt_estimate_reconstruction/incidence_data/jhu_incid_21_11_23.csv")


# join incidence time series
incidence_data <- rki_incid %>%
  inner_join(who_incid, by = "date") %>%
  inner_join(jhu_incid, by = "date") %>%
  inner_join(rki_incid_nowcast, by = "date") %>%
  dplyr::filter(date<=as_date("2021-07-10"))
names(incidence_data) <- c("date", "RKI, positive test", "WHO", "JHU", "RKI, symptom onset")

# reshape data
incid  <- incidence_data %>%
  gather("variable", "value", 2:dim(incidence_data)[2]) %>%
  filter(date >= as.Date("2021-01-01")) %>%
  filter(date <= as.Date("2021-06-10"))

col_values <- get_colors(c("JHU", "RKI, positive test", "WHO", "RKI, symptom onset"), palette = "incidence data", name_consensus = "RKI, positive test")

# plot
R_plot <- ggplot(data = incid, aes(x = date, y = value)) +
  labs(x = NULL, y = "new infections") +
  scale_x_date(limits = as.Date(c(min(incid$date),max(incid$date))),
               breaks = seq(as.Date(min(incid$date)), as.Date(max(incid$date)), "month"),
               date_labels = "%B %d", expand = c(0,1)) +
  theme_minimal(base_family = "serif") +
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
  scale_color_manual(values=col_values, name="data source")

print(R_plot)
ggsave(R_plot, filename = "Figures/incidence_data.pdf",  bg = "transparent",
       width = 13.1, height = 5.8)
