library(lubridate)
library(dplyr)
library(ggplot2)

setwd("../..")
getwd()

source("Rt_estimate_reconstruction/load_data.R")
source("Rt_estimate_reconstruction/prepared_plots.R")

path_estimates <- "reproductive_numbers/data-processed/"
methods <- c("RKI_7day", "ETHZ_sliding_window", "ilmenau", "sdsc", "globalrt_7d", "epiforecasts", "rtlive", "Braunschweig")

start_plot <- as_date("2020-10-01")
end_plot <- as_date("2021-09-30")

target_dates <- seq(start_plot, end_plot, by = "day")

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
  mutate(value = ifelse(value, 1, 2)) %>%
  mutate(method = plyr::mapvalues(method,
                                  c("Braunschweig", "ETHZ_sliding_window", "globalrt_7d", "ilmenau", "RKI_7day", "sdsc"),
                                  c("HZI",          "ETH",                 "globalrt",    "Ilmenau", "RKI",      "SDSC"))) %>%
  arrange(desc(method))

methods_legend <- unique(plot_dates$method)
col_values <- get_colors(methods = methods_legend, palette = "Set2")

plot_dates <- plot_dates %>%
  mutate(method = factor(method, levels = methods_legend))

# add column whether data is used in Fig. 5
plot_dates$used <- FALSE

start_default <- "2020-10-01"
start_globalrt <- "2021-02-15"
start_ilmenau <- "2020-11-16"
days_until_final <- 70

for (method in unique(plot_dates$method)){
  if (method == "globalrt"){
    start_date <- as_date(start_globalrt)
  } else if (method == "Ilmenau"){
    start_date <- as_date(start_ilmenau)
  } else {
    start_date <- as_date(start_default)
  }
  end_date <- as_date("2021-09-30") - days(days_until_final)
  if (method == "rtlive") end_date <- as_date("2021-07-31") - days(days_until_final)
  plot_dates[(plot_dates$method == method) &
               (plot_dates$date >= start_date) &
               (plot_dates$date <= end_date), "used"] <- TRUE
}

plot <- ggplot(data=plot_dates, aes(x=date, y=method, col=method, size=value)) +
  theme_minimal(base_family = "serif") +
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
  geom_point(shape=108) +
  scale_size(range = c(0,2)) +
  scale_x_date(date_labels = "%b %y",
               limits = c(start_plot, end_plot), expand = expansion(0.01)) +
  scale_color_manual(values=col_values, name="method") +
  geom_line(data=plot_dates[plot_dates$used,], aes(x=date, y=method),
            color="black", size=10, alpha=0.15)
  
print(plot)
ggsave(plot, filename = paste0("Figures/missing_pub_dates.pdf"),
       bg = "transparent", width = 14, height = 5)
