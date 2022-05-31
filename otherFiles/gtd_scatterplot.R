library(ggplot2)
library(dplyr)

setwd("../..")
getwd()

# parameter combinations used in papers
methods <- c("ETH(CH)", "RKI",      "Ilmenau", "SDSC",  "AT",    "epiforecasts", "rtlive",  "globalrt",    "HZI", "FR",    "IT",    "NO", "DK", "SE", "PT",  "BE", "NL",    "NL (Omicron)", "consensus")
gt_dist <- c("gamma",   "constant", "ad hoc",  "gamma", "gamma", "gamma",        "lognorm", "exponential", "?",   "gamma", "gamma", "?",  "?",  "?",  "?",   "?",  "gamma", "gamma",        "gamma")
mean_gt <- c(4.8,        4,          5.6,       4.8,     3.4,     3.6,            4.7,       7,             10.5,  7,       6.7,     7.5,  4.7,  4.8,  3.96,  4.7,  4,       3.5,            4)
sd_gt <-   c(2.3,        0,          4.2,       2.3,     1.8,     3.1,            2.9,       7,             8,     4.5,     4.9,     2.9,  2.9,  2.3,  4.74,  2.9,  2,       1.75,           4)

gtds <- data.frame(gtd=gt_dist, gt_mean=mean_gt, gt_sd=sd_gt)
rownames(gtds) <- methods

gtd_scatter <- gtds %>%
  tibble::rownames_to_column(var = "group") %>%
  select(group, gt_mean, gt_sd) %>%
  group_by(gt_mean, gt_sd) %>% 
  mutate(groups = paste0(group, collapse = "/")) %>%
  select(!group) %>%
  distinct() %>%
  mutate(grey = (nchar(groups)==2)) %>%
  tibble::column_to_rownames(var = "groups")
gtd_scatter["NL (Omicron)", "grey"] <- T

ages_first <- c(4.46, 2.63)
names(ages_first) <- c("gt_mean", "gt_sd")
ages_data <- data.frame(gt_mean=ages_first["gt_mean"], gt_sd=ages_first["gt_sd"])

# plot assumed mean and sd of generation times
gtd_scatterplot <- ggplot(data = gtd_scatter, aes(x = gt_mean, y = gt_sd, color = grey)) +
  labs(x = "mean", y = "standard deviation") +
  theme_minimal() +
  theme(
    axis.text=element_text(size=16),
    axis.title=element_text(size=18),
    axis.line = element_line(),
    axis.line.y.right = element_line(),
    axis.line.x.top = element_line(),
    legend.position = "none",
    panel.grid.minor.x = element_line(colour = "transparent"),
    panel.background = element_rect(fill = "transparent")
  ) +
  geom_point(size = 2) +
  geom_text(label = rownames(gtd_scatter),
            position = position_nudge(y = ifelse(rownames(gtd_scatter) %in% c("AT", "rtlive/DK/BE"), +0.33, -0.28)),
            size = 4) + 
  scale_color_manual(values=c("#000000", "#909090")) +
  
  geom_point(data = data.frame(gt_mean=4, gt_sd=4), colour = "red", size = 3) + 
  
  geom_point(data = ages_data, colour = "#909090", shape=21, size = 2) +
  geom_text(mapping = aes(x=ages_data$gt_mean, y=ages_data$gt_sd),
            label="(AT)", color="#909090",
            nudge_y = -0.25, check_overlap = TRUE, size = 2.5) +
  geom_segment(aes(x=ages_first["gt_mean"], y=ages_first["gt_sd"],
                   xend=gtds["AT","gt_mean"], yend=gtds["AT","gt_sd"]),
               arrow=arrow(length=unit(2.5, "mm")), color="#909090", size=0.01, linetype=2) +
  
  coord_cartesian(xlim = c(3.25, 10.5))

print(gtd_scatterplot)

ggsave(gtd_scatterplot, filename = "Figures/gtd_scatterplot.pdf",  bg = "transparent",
       width = 8, height = 5)
