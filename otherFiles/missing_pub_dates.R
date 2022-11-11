library(lubridate)
library(dplyr)
library(ggplot2)

setwd("../..")
getwd()

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

path_estimates <- "reproductive_numbers/data-processed/"
methods <- c("RKI_7day", "ETHZ_sliding_window", "ilmenau", "sdsc", "globalrt_7d", "epiforecasts", "rtlive", "Braunschweig")

start_date <- as_date("2020-10-01")
end_date <- as_date("2021-10-28")
target_dates <- seq(start_date, end_date, by = "day")

df_missing_dates <- data.frame(matrix(rep(NA, length(methods)*length(target_dates)),
                                      ncol = length(methods)),
                               row.names = target_dates) %>%
  setNames(methods)

for (method in methods){
  dates <- list.files(paste0(path_estimates, method), pattern = "\\d{4}-\\d{2}-\\d{2}") %>%
    substr(1,10) %>%
    as_date
  dates <- dates[dates %in% target_dates]
  df_missing_dates[as.character(dates), method] <- F
}
df_missing_dates[is.na(df_missing_dates)] <- T

plot_dates <- df_missing_dates %>%
  rownames_to_column("date") %>%
  mutate(date = as_date(date)) %>%
  gather(method, value, -date) %>%
  filter(value==TRUE) %>%
  mutate(method = plyr::mapvalues(method,
                                  c("Braunschweig", "ETHZ_sliding_window", "globalrt_7d", "ilmenau", "RKI_7day", "sdsc"),
                                  c("HZI",          "ETH",                 "globalrt",    "Ilmenau", "RKI",      "SDSC"))) %>%
  arrange(desc(method))

methods_legend <- unique(plot_dates$method)
col_values <- get_colors(methods = methods_legend, palette = "methods")

plot_dates <- plot_dates %>%
  mutate(method = factor(method, levels = methods_legend))

plot <- ggplot(data=plot_dates, aes(date, method, col=method)) +
  theme_minimal() +
  theme(
    plot.margin = unit(c(1,9,1,2), "mm"),
    plot.title = element_text(size=18),
    axis.text=element_text(size=16),
    axis.title=element_blank(),
    legend.text=element_text(size=16),
    legend.title=element_text(size=18),
    axis.line = element_line(),
    legend.position = "None",
    panel.border = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank()
  ) +
  geom_point() +
  scale_x_date(date_labels = "%b %y",
               limits = c(start_date, end_date), expand = expansion(0.01)) +
  scale_color_manual(values=col_values, name="method")
print(plot)
ggsave(plot, filename = paste0("Figures/missing_pub_dates.pdf"),
       bg = "transparent", width = 14, height = 5)
