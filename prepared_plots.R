library(reshape2)
library(ggplot2)
library(cowplot)
library(viridis)
library(pheatmap)
library(tidyverse)

# set system locale time to English for correct labelling of x axes
Sys.setlocale("LC_TIME", "English")

# function for saving pheatmap
save_pheatmap_pdf <- function(x, filename, width=9, height=4) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

plot_published_vs_calculated <- function(published, calculated, method_name, diff_bounds=c(-0.1, 0.1)){
  
  # combine data
  ts <- full_join(published, calculated, by="date")[, c("date", "R_pub", "R_calc")]
  
  # calcuate difference between estimates
  diff <- data.frame(date=ts$date, difference=ts$R_pub-ts$R_calc)
  diff <- na.omit(diff)
  
  # reshape data
  ts <- melt(ts, id.var='date')
  ts <- na.omit(ts)
  
  # build plot with estimates
  est_plot <- ggplot(data=ts, aes(x=date, y=value, color=variable)) +
    geom_hline(aes(yintercept = 1)) +
    geom_line() +
    xlim(c(min(ts$date), max(ts$date))) + 
    scale_colour_viridis_d(end=0.7, name="Estimates", labels=c("published", "calculated")) +
    labs(title=method_name, x = "date", y = "Rt estimate") +
    theme(legend.position = "top")
  
  # build plot showing differences
  diff_plot <- ggplot(data=diff, aes(x=date, y=difference)) +
    geom_line() +
    xlim(c(min(ts$date), max(ts$date))) + 
    ylim(diff_bounds) +
    labs(x = "date", y = "published - calculated")
  
  # align plots
  plot_grid(est_plot, diff_plot, ncol = 1, align = "v", rel_heights = c(2, 1))
}



plot_multiple_estimates <- function(estimates, legend_name,
                                    include_CI=F, comparing_parameters=F) {
  
  # reshape data
  R_est  <- estimates %>%
    gather("variable", "value", 2:dim(estimates)[2]) %>%
    separate(variable,
             into = c("type", "model"),
             sep = "[.]",
             extra = "merge",
             remove = TRUE) %>%
    spread(type, value)

  if (comparing_parameters) {
    R_est <- R_est %>%
      mutate_at(vars("model"), as.numeric) %>%
      arrange(date, model) %>%
      mutate_at(vars("model"), as.factor)
  }
  
  # plot
  R_plot <- ggplot(data = R_est, aes(x = date, y = R)) +
    geom_hline(aes(yintercept = 1)) +
    labs(x = NULL, y = "Rt estimate") +
    scale_x_date(limits = as.Date(c(min(estimates$date),max(estimates$date)-1)),
                 date_labels = "%b %d", expand = c(0,1)) +
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
      )
  
  if (include_CI){
    R_plot <-  R_plot +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = .3) +
      scale_fill_viridis_d(end = 0.9, name=legend_name)
  } else {
    R_plot <-  R_plot +
      geom_line(aes(group = model, color=model)) +
      scale_colour_viridis_d(end = 0.9, name=legend_name)
  }
  
  return(R_plot)
}





########################################
# first plot plus confidence intervals #
########################################

plot_published_vs_calculated_95CI <- function(published, calculated, method_name, diff_bounds=c(-0.1, 0.1)){
  
  if ("R_pub" %in% names(published)){
    names(published)[2] <- "0.5"
  }
  
  # combine data
  ts <- full_join(published, calculated, by="date", suffix = c(".pub", ".calc"))
  
  # calcuate difference between estimates
  diff <- data.frame(date=ts$date, difference=ts$`0.5.pub`-ts$`0.5.calc`)
  diff <- na.omit(diff)
  
  # reshape data
  ts_median <- melt(ts[, c("date", "0.5.pub", "0.5.calc")], id.var='date')
  ts_median <- na.omit(ts_median)
  
  ts_quantiles <- ts[, c("date", "0.025.pub", "0.975.pub", "0.025.calc", "0.975.calc")]
  names(ts_quantiles) <- c("date", "lower.pub", "upper.pub", "lower.calc", "upper.calc")
  
  # build plot with estimates
  est_plot <- ggplot() +
    geom_ribbon(data=ts_quantiles ,aes(x=date, ymax=upper.pub, ymin=lower.pub), fill=viridis(1), alpha=.4) +
    geom_ribbon(data=ts_quantiles ,aes(x=date, ymax=upper.calc, ymin=lower.calc), fill=viridis(1, begin=0.7), alpha=.4) +
    geom_hline(aes(yintercept = 1)) +
    geom_line(data=ts_median, aes(x=date, y=value, color=variable)) +
    xlim(c(min(ts_median$date), max(ts_median$date))) + 
    scale_colour_viridis_d(end=0.7, name="Estimates", labels=c("published", "calculated")) +
    labs(title=method_name, x = "date", y = "Rt estimate") +
    theme(legend.position = "top")
  
  # build plot showing differences
  diff_plot <- ggplot(data=diff, aes(x=date, y=difference)) +
    geom_line() +
    xlim(c(min(ts_median$date), max(ts_median$date))) + 
    ylim(diff_bounds) +
    labs(x = "date", y = "published - calculated")
  
  # align plots
  plot_grid(est_plot, diff_plot, ncol = 1, align = "v", rel_heights = c(5, 1))
}



##############################
# plots for final comparison #
##############################

plot_for_comparison <- function(estimates, comp_methods,
                                legend_name="model", filenames = "latest_plot.png",
                                include_CI=F, method, variation, comparing_parameters=F){
  if (length(comp_methods)*3 == dim(estimates)[2]-1){
    names_ci <- rep(NA, length(comp_methods)*3)
    for (i in 1:length(comp_methods)){
      names_ci[(i-1)*3 + 1] <- paste0("R.", comp_methods[i])
      names_ci[(i-1)*3 + 2] <- paste0("lower.", comp_methods[i])
      names_ci[(i-1)*3 + 3] <- paste0("upper.", comp_methods[i])
    }
    names(estimates) <- c("date", names_ci)
  } else {
    names_R <- rep(NA, length(comp_methods))
    for (i in 1:length(comp_methods)){
      names_R[i] <- paste0("R.", comp_methods[i])
    }
    names(estimates) <- c("date", names_R)
  }
  
  estimates <- estimates[estimates$date < "2021-06-15",]
  estimates <- estimates[rowSums(is.na(estimates)) == 0,]
  latest_estimates <- estimates[estimates$date >= "2021-01-13",]
  
  R_plot <- plot_multiple_estimates(estimates, legend_name, include_CI = include_CI, comparing_parameters = comparing_parameters)
  R_plot_latest <- plot_multiple_estimates(latest_estimates, legend_name, include_CI = include_CI, comparing_parameters = comparing_parameters)
  
  #ggsave(R_plot, filename = "estimates_plot.png",  bg = "transparent")
  #print(R_plot)
  ggsave(R_plot_latest, filename = paste0("Figures/estimates", filenames),  bg = "transparent",
         width = 13.1, height = 6.3)
  print(R_plot_latest) 
  
  if(!comparing_parameters) {
    n <- length(comp_methods)
    matr <- matrix(rep(rep(0,n), n), ncol=n)
    corr <- matrix(rep(rep(0,n), n), ncol=n)
    colnames(matr) <- rownames(matr) <- colnames(corr) <- rownames(corr) <- comp_methods
    
    for (method1 in comp_methods) {
      for (method2 in comp_methods){
        diff <- estimates[,paste0("R.", method1)] - estimates[,paste0("R.", method2)]
        matr[method1, method2] <- mean(as.vector(abs(diff)))
        corr[method1, method2] <- cor(estimates[,paste0("R.", method1)], estimates[,paste0("R.", method2)])
      }
    }
    
    mean_abs_diff <- pheatmap(matr, color = viridis(100), breaks = seq(0,0.3,0.3/100),
                              border_color = NA, display_numbers = TRUE,
                              fontsize = 18, fontsize_number=22, number_color = "white",
                              angle_col = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE)#,
                              #main = paste("Mean absolute differences between estimates over",
                              #             max(estimates$date) - min(estimates$date),
                              #             "days using", method, variation))
    save_pheatmap_pdf(mean_abs_diff, paste0("Figures/mean_abs_diff", filenames))
    
    
    correlations <- pheatmap(corr, color = viridis(100, direction = -1), breaks = seq(0.25,1,0.75/100),
                             border_color = NA, display_numbers = TRUE,
                             fontsize = 18, fontsize_number=22, number_color = "white",
                             angle_col = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE)#,
                             #main = paste("Correlations between estimates over",
                             #             max(estimates$date) - min(estimates$date),
                             #             "days using", method, variation))
    save_pheatmap_pdf(correlations, paste0("correlations", filenames))
  }
}