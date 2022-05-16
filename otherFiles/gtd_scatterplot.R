library(ggplot2)

setwd("../..")
getwd()

# parameter combinations used in papers
gt_dist <- c("gamma", "constant", "ad hoc", "gamma", "gamma", "gamma", "lognorm", "exponential", "?",  "gamma")
mean_gt <- c(4.8,      4,          5.6,      4.8,    3.4,     3.6,     4.7,       7,             10.5, 4)
sd_gt <-   c(2.3,      0,          4.2,      2.3,    1.8,     3.1,     2.9,       7,             8,    4)

gtds <- data.frame(gtd=gt_dist, gt_mean=mean_gt, gt_sd=sd_gt)
methods <- c("ETH", "RKI", "Ilmenau", "SDSC", "AGES", "epiforecasts", "rtlive", "globalrt", "HZI", "consensus")
rownames(gtds) <- methods

ages_first <- c(4.46, 2.63)
names(ages_first) <- c("gt_mean", "gt_sd")
ages_data <- data.frame(gt_mean=ages_first["gt_mean"], gt_sd=ages_first["gt_sd"])

# plot assumed mean and sd of generation times
gtd_scatterplot <- ggplot(data = gtds[-which(rownames(gtds)=="rtlive"), ], aes(x = gt_mean, y = gt_sd)) +
  labs(x = "mean", y = "standard deviation") +
  theme_minimal() +
  theme(
    axis.text=element_text(size=16),
    axis.title=element_text(size=18),
    axis.line = element_line(),
    axis.line.y.right = element_line(),
    axis.line.x.top = element_line(),
    legend.position = "bottom",
    panel.grid.minor.x = element_line(colour = "transparent"),
    panel.background = element_rect(fill = "transparent")
  ) +
  geom_point(size = 2) +
  geom_text(label=c("ETH/SDSC", "RKI", "Ilmenau", "SDSC", "AGES",
                    "epiforecasts", "globalrt", "HZI", "consensus"),
            nudge_y = -0.28, check_overlap = TRUE, size = 5) + 
  
  geom_point(data = gtds["rtlive", c("gt_mean", "gt_sd")], size = 2) +
  geom_text(mapping = aes(x=gtds["rtlive", "gt_mean"], y=gtds["rtlive", "gt_sd"]), label="rtlive",
            nudge_y = +0.33, check_overlap = TRUE, size = 5) +
  
  geom_point(data = data.frame(gt_mean=4, gt_sd=4), colour = "red", size = 3) + 
  
  geom_point(data = ages_data, colour = "black", shape=21, size = 2) +
  geom_text(mapping = aes(x=ages_data$gt_mean, y=ages_data$gt_sd), label="(AGES)",
            nudge_y = -0.25, check_overlap = TRUE, size = 2.5) +
  geom_segment(aes(x=ages_first["gt_mean"], y=ages_first["gt_sd"],
                   xend=gtds["AGES","gt_mean"], yend=gtds["AGES","gt_sd"]),
               arrow=arrow(length=unit(2.5, "mm")), color="black", size=0.01, linetype=2)

print(gtd_scatterplot)

ggsave(gtd_scatterplot, filename = "Figures/gtd_scatterplot.pdf",  bg = "transparent",
       width = 8, height = 5)
