library(ggplot2)
library(dplyr)

setwd("../..")
getwd()

# parameter combinations used in papers
methods <- c("ETH(CH)", "RKI",      "Ilmenau", "SDSC",  "AT",    "epiforecasts", "rtlive",  "globalrt",    "HZI", "FR",    "IT",    "NO", "DK", "SE", "PT",  "BE", "NL",    "CZ",      "SCT",   "SI",    "consensus")
gt_dist <- c("gamma",   "constant", "ad hoc",  "gamma", "gamma", "gamma",        "lognorm", "exponential", "?",   "gamma", "gamma", "?",  "?",  "?",  "?",   "?",  "gamma", "uniform", "gamma", "gamma", "gamma")
mean_gt <- c(4.8,        4,          5.6,       4.8,     3.4,     3.6,            4.7,       7,             10.3,  7,       6.7,     7.5,  4.7,  4.8,  3.96,  4.7,  3.5,     5.5,       6.5,     6.5,     4)
sd_gt <-   c(2.3,        0,          4.2,       2.3,     1.8,     3.1,            2.9,       7,             7.6,   4.5,     4.9,     2.9,  2.9,  2.3,  4.74,  2.9,  1.75,    1.12,      4.11,    4.11,    4)

gtds <- data.frame(gtd=gt_dist, gt_mean=mean_gt, gt_sd=sd_gt)
rownames(gtds) <- methods

gtd_scatter <- gtds %>%
  tibble::rownames_to_column(var = "group") %>%
  dplyr::select(group, gt_mean, gt_sd) %>%
  group_by(gt_mean, gt_sd) %>% 
  mutate(groups = paste0(group, collapse = "/")) %>%
  dplyr::select(!group) %>%
  distinct() %>%
  mutate(grey = (nchar(groups)==2 | groups=="SCT/SI")) %>%
  tibble::column_to_rownames(var = "groups")

AT_first <- data.frame(gt_mean=4.46, gt_sd=2.63)
NL_first <- data.frame(gt_mean=4, gt_sd=2)

# plot assumed mean and sd of generation times
gtd_scatterplot <- ggplot(data = gtd_scatter, aes(x = gt_mean, y = gt_sd, color = grey)) +
  labs(x = "mean", y = "standard deviation") +
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_line(colour = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    axis.text=element_text(size=16),
    axis.title=element_text(size=18),
    axis.line = element_line(),
    axis.line.y.right = element_line(),
    axis.line.x.top = element_line(),
    legend.position = "none"
  ) +
  geom_point(size = 2) +
  geom_text(label = rownames(gtd_scatter),
            position = position_nudge(y = ifelse(rownames(gtd_scatter) %in% c("AT", "RKI", "rtlive/DK/BE"), +0.33, -0.28)),
            size = 4) + 
  scale_x_continuous(breaks = seq(0,10,2)) +
  scale_color_manual(values=c("#000000", "#909090")) +
  
  geom_point(data = data.frame(gt_mean=4, gt_sd=4), colour = "red", size = 3) + 
  
  geom_point(data = AT_first, colour = "#909090", shape=21, size = 2) +
  geom_text(mapping = aes(x=AT_first$gt_mean, y=AT_first$gt_sd),
            label="(AT)", color="#909090",
            nudge_y = -0.24, check_overlap = TRUE, size = 2.5) +
  geom_segment(aes(x=AT_first$gt_mean, y=AT_first$gt_sd,
                   xend=gtds["AT","gt_mean"], yend=gtds["AT","gt_sd"]),
               arrow=arrow(length=unit(2.5, "mm")), color="#909090", size=0.01, linetype=5) +
  
  geom_point(data = NL_first, colour = "#909090", shape=21, size = 2) +
  geom_text(mapping = aes(x=NL_first$gt_mean, y=NL_first$gt_sd),
            label="(NL)", color="#909090",
            nudge_y = -0.24, check_overlap = TRUE, size = 2.5) +
  geom_segment(aes(x=NL_first$gt_mean, y=NL_first$gt_sd,
                   xend=gtds["NL","gt_mean"], yend=gtds["NL","gt_sd"]),
               arrow=arrow(length=unit(2.5, "mm")), color="#909090", size=0.01, linetype=5) +
  
  coord_cartesian(xlim=c(0,10.6), ylim=c(-0.2,7.9), expand=F)

print(gtd_scatterplot)

ggsave(gtd_scatterplot, filename = "Figures/gtd_scatterplot.png",  bg = "transparent",
       width = 10, height = 5)
