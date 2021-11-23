library(ggplot2)

# parameter combinations used in papers
gt_dist <- c("gamma", "constant", "ad hoc", "gamma", "gamma", "gamma", "lognorm", "exponential")
mean_gt <- c(4.8,      4,          5.6,      4.8,    3.4,     3.6,     4.7,       7)
sd_gt <-   c(2.3,      0,          4.2,      2.3,    1.8,     3.1,     2.9,       7)

gtds <- data.frame(gtd=gt_dist, gt_mean=mean_gt, gt_sd=sd_gt)
methods <- c("ETH", "RKI", "Ilmenau", "SDSC", "AGES", "epiforecasts", "rtlive", "globalrt")
rownames(gtds) <- methods

# plot assumed mean and sd of generation times
gtd_scatterplot <- ggplot(data = gtds, aes(x = gt_mean, y = gt_sd)) +
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
  geom_text(label=c("ETH/SDSC", "RKI", "Ilmenau", "SDSC", "AGES", "epiforecasts", "rtlive", "globalrt"),
            nudge_y = 0.3, check_overlap = TRUE, size = 5) + 
  geom_point(data = data.frame(gt_mean=4, gt_sd=4), colour = "red", size = 3)
print(gtd_scatterplot)
ggsave(gtd_scatterplot, filename = "D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code/Figures/gtd_scatterplot.pdf",  bg = "transparent",
       width = 8, height = 5)
