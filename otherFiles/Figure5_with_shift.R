library(ggpubr)

setwd("../..")
# needs to be the directory with the repos "Rt_estimate_reconstruction", "reproductive_numbers" 
getwd()

source("Rt_estimate_reconstruction/prepared_plots.R")

conf_level <- "95"
days_until_final <- 70
methods_CI <- c("epiforecasts", "ETHZ_sliding_window", "globalrt_7d",
                "ilmenau", "RKI_7day", "rtlive", "SDSC", "ETH_old", "ETH_new")
methods <- c("Braunschweig", "epiforecasts", "ETHZ_sliding_window", "globalrt_7d",
             "ilmenau", "RKI_7day", "rtlive", "SDSC", "ETH_old", "ETH_new")

labels <- read_csv(paste0("Rt_estimate_reconstruction/otherFiles/consistence_measures/",
                          days_until_final, "/estimate_labels.csv")) %>%
  as.data.frame() %>%
  dplyr::select(!num_est) %>%
  rename(method = ...1) %>%
  mutate(method = plyr::mapvalues(method,
                                  c("Braunschweig", "ETHZ_sliding_window", "globalrt_7d", "ilmenau", "RKI_7day"),
                                  c("HZI",          "ETH",                 "globalrt",    "Ilmenau", "RKI")))

labels[nrow(labels)+1:2,] <- labels[labels$method=="ETH",]
labels[nrow(labels)-1,"method"] <- "ETH_old"
labels[nrow(labels),"method"] <- "ETH_new"

labels <- labels %>%
  pivot_longer(!method, names_to = "variable", values_to = "label")

optimal_shift <- c(10, 10, 10, 4, 7, 7, 10, 19, 3, -3)
names(optimal_shift) <- c("ETH", "ETH_new", "ETH_old", "RKI", "Ilmenau", "SDSC",
                          "epiforecasts", "rtlive", "globalrt", "HZI")

##################
# coverage rates #
##################

CI_coverage <- read_csv(paste0("Rt_estimate_reconstruction/otherFiles/consistence_measures/",
                               days_until_final, "/",
                               conf_level, "_CI_coverage.csv")) %>%
  as.data.frame() %>%
  column_to_rownames("...1")

coverage_data <- CI_coverage[, c(as.character(0:39), "min_lag")] %>%
  rownames_to_column("method") %>%
  dplyr::filter(method %in% methods) %>%
  mutate(method = plyr::mapvalues(method,
                                  c("ETHZ_sliding_window", "globalrt_7d", "ilmenau", "RKI_7day"),
                                  c("ETH",                 "globalrt",    "Ilmenau", "RKI"))) %>%
  arrange(method)

coverage_data <- coverage_data %>%
  gather("variable", "value", 2:(dim(coverage_data)[2] - 1))

coverage_data <- merge(coverage_data, labels, by=c("method", "variable")) %>%
  mutate(variable = -1 * as.numeric(variable),
         label = ifelse(label=="estimate", 1,
                        ifelse(label=="estimate based on partial data", 2,
                               ifelse(label=="forecast", 3, NA)))) %>%
  mutate(variable = variable+optimal_shift[method]) %>%
  filter(variable %in% -20:0)

coverage_data_ETH <- coverage_data %>% filter(method %in% c("ETH_old", "ETH_new"))
coverage_data <- coverage_data %>% filter(!(method %in% c("ETH_old", "ETH_new")))

coverage_plot <- ggplot() +
  theme_minimal(base_family="serif") +
  theme(
    plot.margin = unit(c(3,14,2,3), "mm"),
    plot.title = element_text(size=18),
    axis.text=element_text(size=16),
    axis.title.y=element_text(size=18),
    axis.title.x = element_blank(),
    legend.text=element_text(size=16),
    legend.title=element_text(size=18),
    legend.position = "bottom",
    panel.border= element_rect(fill = "transparent", size = 0.5),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank()
  ) +
  ylab("coverage of consolidated est.") +
  coord_cartesian(xlim = c(-20.3, 0.3), ylim = c(-0.05, 1.05), expand = FALSE, clip = "off") +
  scale_x_continuous(labels = paste0(seq(20, 0, -5), "d back"))

methods_legend <- unique(coverage_data$method)
col_values <- get_colors(methods = methods_legend, palette = "methods")

for (l in 1:3){
  coverage_plot <-  coverage_plot +
    geom_line(data=coverage_data[coverage_data$label<=l,],
              aes(x = variable, y = value, color=method),
              linewidth = 1.5, linetype=l, na.rm = T)
}

coverage_plot <-  coverage_plot +
  geom_text(data=subset(coverage_data, (variable == -1 * min_lag + optimal_shift[method]) |
                          (variable == 0 & !(method %in% c("ETH", "Ilmenau", "RKI")))),
            aes(label = method, color = method,
                x = ifelse(method%in%c("ETH", "Ilmenau", "RKI"), variable, variable-0.8),
                y = ifelse(method=="globalrt", 1.025, 
                           ifelse(method=="epiforecasts", value+0.07, value-0.02))),
            size = 4, family="serif",
            show.legend = FALSE)

coverage_plot <- coverage_plot +
  geom_line(data=coverage_data_ETH,
            aes(x = variable, y = value, group=method),
            size = 1, color=col_values["ETH"],
            na.rm = T, show.legend=FALSE) +
  geom_text(data=subset(coverage_data_ETH, variable == -4),
            aes(label = ifelse(method == "ETH_old", "old", "new"),
                x = variable+0.6,
                y = value),
            color = col_values["ETH"],
            size = 4, family="serif",
            show.legend = FALSE)

coverage_plot <- coverage_plot +
  scale_color_manual(values=col_values, name="method")

coverage_plot <- coverage_plot +
  geom_hline(yintercept=0.95, linetype="dashed") +
  annotate("text", label = "nominal level: 95%", x = -17.5, y = 0.92, size = 5, family = "serif")

coverage_plot <- coverage_plot +
  guides(color = guide_legend(override.aes = list(size=3, label="")))

print(coverage_plot)

ggsave(coverage_plot, filename = paste0("Figures/CI/", conf_level, "_coverage_rates_shifted.pdf"),
       bg = "transparent", width = 8, height = 5.8)




#############
# CI widths #
#############

CI_width <- read_csv(paste0("Rt_estimate_reconstruction/otherFiles/consistence_measures/",
                            days_until_final, "/", conf_level, "_CI_width.csv")) %>%
  as.data.frame() %>%
  column_to_rownames("...1")

width_data <- CI_width[, c(as.character(0:39), "min_lag")] %>%
  rownames_to_column("method") %>%
  dplyr::filter(method %in% methods_CI) %>%
  mutate(method = plyr::mapvalues(method,
                                  c("ETHZ_sliding_window", "globalrt_7d", "ilmenau", "RKI_7day"),
                                  c("ETH",                 "globalrt",    "Ilmenau", "RKI"))) %>%
  arrange(method)

width_data <- width_data %>%
  gather("variable", "value", 2:(dim(width_data)[2] - 1))

width_data <- merge(width_data, labels, by=c("method", "variable")) %>%
  mutate(variable = -1 * as.numeric(variable),
         label = ifelse(label=="estimate", 1,
                        ifelse(label=="estimate based on partial data", 2,
                               ifelse(label=="forecast", 3, NA)))) %>%
  mutate(variable = variable+optimal_shift[method]) %>%
  filter(variable %in% -20:0)

width_data_ETH <- width_data %>% filter(method %in% c("ETH_old", "ETH_new"))
width_data <- width_data %>% filter(!(method %in% c("ETH_old", "ETH_new")))

width_plot <- ggplot() +
  theme_minimal(base_family="serif") +
  theme(
    plot.margin = unit(c(3,14,2,3), "mm"),
    plot.title = element_text(size=18),
    axis.text=element_text(size=16),
    axis.title.y=element_text(size=18),
    axis.title.x = element_blank(),
    legend.text=element_text(size=16),
    legend.title=element_text(size=18),
    legend.position = "bottom",
    panel.border= element_rect(fill = "transparent", size = 0.5),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank()
  ) +
  ylab("width of 95%-CI") +
  coord_cartesian(xlim = c(-20.3, 0.3), ylim = c(-0.03, 0.83), expand = FALSE, clip = "off") +
  scale_x_continuous(labels = paste0(seq(20, 0, -5), "d back")) +
  scale_y_continuous(labels=function(x) sprintf("    %.1f", x))

methods_legend <- unique(width_data$method)
col_values <- get_colors(methods = methods_legend, palette = "methods")

for (l in 1:3){
  width_plot <-  width_plot +
    geom_line(data=width_data[width_data$label<=l,],
              aes(x = variable, y = value, color=method),
              size = 1.5, linetype=l, na.rm = T)
}

width_plot <-  width_plot +
  geom_text(data=subset(width_data, (variable == -1 * min_lag + optimal_shift[method]) |
                          (variable == 0 & !(method %in% c("ETH", "Ilmenau", "RKI")))),
            aes(label = method, color = method,
                x = ifelse(method%in%c("epiforecasts", "rtlive", "SDSC"), variable-1,
                           ifelse(method=="globalrt", variable-1.5, variable)),
                y = ifelse(method %in% c("ETH", "globalrt", "RKI", "SDSC"), value-0.02, value+0.02)),
            size = 4, family="serif",
            show.legend = FALSE)

width_plot <- width_plot +
  geom_line(data=width_data_ETH,
            aes(x = variable, y = value, group=method),
            size = 1, color=col_values["ETH"],
            na.rm = T, show.legend=FALSE) +
  geom_text(data=subset(width_data_ETH, variable == -4),
            aes(label = ifelse(method == "ETH_old", "old", "new"),
                x = ifelse(method == "ETH_old", variable, variable+0.6),
                y = ifelse(method == "ETH_old", value-0.02, value)),
            color = col_values["ETH"],
            size = 4, family="serif",
            show.legend = FALSE)


width_plot <- width_plot +
  scale_color_manual(values=col_values, name="method") +
  guides(color = guide_legend(override.aes = list(size=3)))

print(width_plot)

ggsave(width_plot, filename = paste0("Figures/CI/", conf_level, "_CI_widths_shifted.pdf"),
       bg = "transparent", width = 8, height = 5.8)

  
  
  
#################
# diff to final #
#################
  
for(diff_type in c("abs_diff", "diff", "switch", "switch2")) {
  ylim <- switch (diff_type,
              "abs_diff" = c(-0.01, 0.175),
              "diff" = c(-0.018, 0.06),
              "switch" = c(-0.01, 0.4),
              "switch2" = c(-0.01, 0.4))
  
  diff_to_final <- read_csv(paste0("Rt_estimate_reconstruction/otherFiles/consistence_measures/",
                                   days_until_final, "/",
                                   diff_type, "_to_final.csv")) %>%
    as.data.frame() %>%
    column_to_rownames("...1") %>%
    dplyr::select(as.character(0:39), "num_est") %>%
    rownames_to_column("method") %>%
    full_join(read_csv(paste0("Rt_estimate_reconstruction/otherFiles/consistence_measures/",
                              days_until_final, "/95_CI_width.csv")) %>%
                as.data.frame() %>%
                column_to_rownames("...1") %>%
                dplyr::select("min_lag") %>%
                rownames_to_column("method"),
              by = "method") %>%
    column_to_rownames("method")
  diff_to_final["Braunschweig", "min_lag"] <- 3
  
  diff_data <- diff_to_final %>%
    rownames_to_column("method") %>%
    dplyr::filter(method %in% methods) %>%
    mutate(method = plyr::mapvalues(method,
                                    c("Braunschweig", "ETHZ_sliding_window", "globalrt_7d", "ilmenau", "RKI_7day"),
                                    c("HZI",          "ETH",                 "globalrt",    "Ilmenau", "RKI"))) %>%
    arrange(method)
  
  diff_data <- diff_data %>%
    mutate(legend = paste0(method, " (n=", num_est, ")")) %>%
    dplyr::select(!num_est) %>%
    gather("variable", "value", as.character(0:39))
  
  diff_data <- merge(diff_data, labels, by=c("method", "variable")) %>%
    mutate(variable = -1 * as.numeric(variable),
           label = ifelse(label=="estimate", 1,
                          ifelse(label=="estimate based on partial data", 2,
                                 ifelse(label=="forecast", 3, NA)))) %>%
    mutate(variable = variable+optimal_shift[method]) %>%
    filter(variable %in% -20:0)
  
  diff_data <- diff_data %>% filter(!(method %in% c("ETH_old", "ETH_new")))
  
  diff_plot_temp <- ggplot() +
    theme_minimal(base_family="serif") +
    theme(
      plot.margin = unit(c(3,14,2,3), "mm"),
      plot.title = element_text(size=18),
      axis.text=element_text(size=16),
      axis.title.y=element_text(size=18),
      axis.title.x = element_blank(),
      legend.text=element_text(size=16),
      legend.title=element_text(size=18),
      legend.position = "bottom",
      panel.border= element_rect(fill = "transparent", size = 0.5),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major = element_line(),
      panel.grid.minor = element_blank()
    ) +
    ylab(switch(diff_type,
                "diff" = "MSD to consolidated est.",
                "abs_diff" = "MAD to consolidated est.",
                "switch" = "proportion flip above / below 1",
                "switch2" = "proportion flip above / below [0.97, 1.03]")) +
    coord_cartesian(xlim = c(-20.3, 0.3), ylim = ylim, expand = FALSE) +
    scale_x_continuous(labels = paste0(seq(20, 0, -5), "d back"))
  
  if (diff_type == "diff") {
    diff_plot_temp <- diff_plot_temp + geom_hline(yintercept = 0)
  }
  
  # get colors and legend names
  methods_legend <- unique(diff_data$method)
  col_values <- get_colors(methods = methods_legend, palette = "methods")
  
  for (l in 1:3){
    diff_plot_temp <-  diff_plot_temp +
      geom_line(data=diff_data[diff_data$label<=l,],
                aes(x = variable, y = value, color=method),
                size = 1.5, linetype=l, na.rm = T)
  }
  
  if(diff_type=="abs_diff"){
    diff_plot_temp <-  diff_plot_temp +
      geom_text(data=subset(diff_data, (variable == -1 * min_lag + optimal_shift[method]) |
                              (variable == 0 & !(method %in% c("ETH", "Ilmenau", "RKI")))),
                aes(label = method, color = method,
                    x = ifelse(method %in% c("epiforecasts", "SDSC", "globalrt"), variable-0.8, 
                               ifelse(method=="Ilmenau", variable-1.4, 
                                      ifelse(method=="RKI", variable, 
                                             ifelse(method=="rtlive", -10, variable+0.7)))),
                    y = ifelse(method == "Ilmenau", 0.17, 
                               ifelse(method%in%c("epiforecasts", "RKI"), value-0.01, 
                                      ifelse(method=="rtlive", 0.046, value+0.001)))),
                size = 4, family="serif",
                show.legend = FALSE)
  } else {
    diff_plot_temp <-  diff_plot_temp +
      geom_text(data=subset(diff_data, (variable == -1 * min_lag + optimal_shift[method]) |
                              (variable == 0 & !(method %in% c("ETH", "Ilmenau", "RKI")))),
                aes(label = method, color = method,
                    x = ifelse(method %in% c("globalrt", "SDSC", "rtlive", "epiforecasts"), variable-0.8, 
                               ifelse(method=="Ilmenau", variable-1.8, variable+0.3)),
                    y = ifelse(method == "Ilmenau", 0.055,
                               ifelse(method =="rtlive", value+0.003,
                                      ifelse(method=="SDSC", value-0.003, 
                                             ifelse(method=="HZI", value, value+0.001))))),
                size = 4, family="serif",
                show.legend = FALSE)
  }
  
  diff_plot_temp <- diff_plot_temp + 
    scale_color_manual(values=col_values, name="method", labels=unique(diff_data$legend)) +
    guides(color = guide_legend(override.aes = list(size=3)))
  
  if (dim(subset(diff_data, value > ylim[2]))[1] > 0) {
    diff_plot_temp <- diff_plot_temp + 
      geom_text(data=subset(diff_data, value > ylim[2]),
                aes(label = round(value, 2),
                    x = variable,
                    y = ylim[2] - 0.06 * (ylim[2] + abs(ylim[1])),
                    color = method),
                angle=90,
                size = 5, family="serif",
                show.legend = FALSE)
  }
  if (dim(subset(diff_data, value < ylim[1]))[1] > 0) {
    diff_plot_temp <- diff_plot_temp + 
      geom_text(data=subset(diff_data, value < ylim[1]),
                aes(label = round(value, 2),
                    x = variable,
                    y = ylim[1] + 0.06 * (ylim[2] + abs(ylim[1])),
                    color = method),
                angle=90,
                size = 5, family="serif",
                show.legend = FALSE)   
  }
  
  print(diff_plot_temp)
  
  ggsave(diff_plot_temp, filename = paste0("Figures/CI/", diff_type, "_to_final.pdf"),
         bg = "transparent", width = 8, height = 5.8)
  
  assign(paste0(diff_type, "_plot"), diff_plot_temp)
}


################
# put together #
################

consistence_plot <- ggarrange(coverage_plot, width_plot, abs_diff_plot, diff_plot,
                              switch_plot, switch2_plot, ncol=2, nrow=3,
                              labels = list("A", "B", "C", "D", "E", "F"),
                              font.label = list(size = 18, face = "bold", family = "serif"),
                              common.legend = T, legend="bottom",
                              legend.grob = get_legend(abs_diff_plot))
print(consistence_plot)
ggsave(consistence_plot, filename = "Figures/CI/consistence_plots_with_shift.pdf",
       bg = "transparent", width = 16, height = 13.6)


