library(reshape2)
library(ggplot2)
library(cowplot)
library(viridis)
library(pheatmap)


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



plot_multiple_estimates <- function(estimates, methods, include_CI=F) {
  # TODO: estimates einzeln übergeben und abhängig von Anzahl trotzdem richtig joinen.
  estimates <- as.data.frame(estimates)
  
  # pull first mean/median of estimates
  R_est <- estimates[, 1:2]
  names(R_est)[2] <- paste0("R.", methods[1])
  
  # add median/mean of further estimates in correct order for coloring
  for (i in 2:length(methods)) {
    name <- paste0("R.", methods[i])

    if (dim(estimates)[2] == (1 + 3*length(methods))){
      #view(estimates)
      R_est[[name]] <- estimates[, (2 + (i-1)*3)]
    } else if (dim(estimates)[2] == (1 + length(methods))){
      R_est[[name]] <- estimates[, 1 + i]
    } else {
      print("Make sure to pass 1 or 3 columns per method (R, (lb_ci, ub_ci))")
    } 
  }
  
  #View(R_est)
  
  # reshape data
  R_est <- melt(R_est, id.var='date')
  R_est <- na.omit(R_est)
  
  R_plot <- ggplot()
  
  if(include_CI){
    
    if (dim(estimates)[2] != (1 + 3*length(methods))){
      print("Make sure to pass 3 columns per method (R, lb_ci, ub_ci)")
      
    } else{
      # select CI's of estimates
      R_CI <- estimates %>%
        dplyr::select(starts_with(c("date", "lower", "upper", "0.025", "0.975")))
      
      View(R_CI)
      
      for (i in 1:length(methods)) {
        lower <- paste0("lower.", methods[i])
        upper <-paste0("upper.", methods[i])
        names(R_CI)[i*2] <- lower
        names(R_CI)[i*2+1] <- upper
        
        color <- viridis(1, begin = (i - 1) * 0.8 / (length(methods) - 1))
        
        # add shades to plot
        R_plot <- R_plot +
          geom_ribbon(data=R_CI ,aes(x=date, ymax=.data[[upper]], ymin=.data[[lower]]), fill=color, alpha=.3)
      }
      
      View(R_CI)
    }
  }
  
  R_plot <- R_plot +
    geom_hline(aes(yintercept = 1)) +
    geom_line(data=R_est, aes(x=date, y=value, color=variable)) +
    scale_colour_viridis_d(name="method", labels=methods) +
    labs(x = "date", y = "Rt estimate") +
    theme(legend.position = "top")
  print(R_plot)
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

plot_for_comparison <- function(estimates, comp_methods, method, variation){
  names(estimates) <- c("date", comp_methods)
  
  estimates <- estimates[rowSums(is.na(estimates)) == 0,]
  latest_estimates <- estimates[estimates$date >= "2021-04-01",]
  
  plot_multiple_estimates(estimates, methods = comp_methods)
  plot_multiple_estimates(latest_estimates, methods = comp_methods)
  
  n <- dim(estimates)[2] - 1
  matr <- matrix(rep(rep(0,n), n), ncol=n)
  corr <- matrix(rep(rep(0,n), n), ncol=n)
  colnames(matr) <- rownames(matr) <- colnames(corr) <- rownames(corr) <- comp_methods
  
  par(mfrow=c(6,6))
  for (method1 in comp_methods) {
    for (method2 in comp_methods){
      diff <- estimates[,method1] - estimates[,method2]
      if (plot_all_differences){
        plot(diff, type="l", ylab = "diff", ylim=c(-0.5, 0.5),
             main = paste(method1, "vs.", method2))
      }
      matr[method1, method2] <- mean(abs(diff))
      corr[method1, method2] <- cor(estimates[,method1], estimates[,method2])
    }
  }
  par(mfrow=c(1,1))
  
  pheatmap(matr, color = viridis(100), breaks = seq(0,0.4,0.4/100),
           border_color = NA, display_numbers = TRUE,
           fontsize = 12, fontsize_number=20, number_color = "white",
           angle_col = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
           main = paste("Mean absolute differences between estimates over",
                        max(estimates$date) - min(estimates$date),
                        "days using", method, variation))
  
  pheatmap(corr, color = viridis(100, direction = -1), breaks = seq(0.5,1,0.5/100),
           border_color = NA, display_numbers = TRUE,
           fontsize = 12, fontsize_number=20, number_color = "white",
           angle_col = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE,
           main = paste("Correlations between estimates over",
                        max(estimates$date) - min(estimates$date),
                        "days using", method, variation))
}