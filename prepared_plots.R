library(reshape2)
library(ggplot2)
library(cowplot)
library(viridis)
library(scales)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(lubridate)


# set system locale time to English for correct labelling of x axes
Sys.setlocale("LC_TIME", "English")

get_colors <- function(methods, palette, name_consensus = "consensus"){
  all_methods <- c("AGES", "HZI", "consensus", "epiforecasts", "ETH", "globalrt", "Ilmenau", "RKI", "RKI_4day", "rtlive", "SDSC")
  num_est <- length(methods)
  
  # constant colors for different models
  if (all(methods %in% all_methods)){
    cols <- c("#E41A1C", # red
              "#f0bc00", # yellow
              "#000000", # black
              "#ff8b00", # orange
              "#00d333", # light green
              "#00922c", # dark green
              "#0397d8", # light blue
              "#1502a0", # medium blue
              "#230082", # dark blue 
              "#bb2881", # berry
              "#7c007c") # violett
    names(cols) <- all_methods
    cols <- cols[methods]
  
  # if all EpiEstim use color palette
  } else {
    if (palette == "YlGn"){
      # avoid problem of too short color palette
      start_col <- 1
      cols <- brewer.pal(name=palette,n=num_est+start_col)[start_col+(1:num_est)]
      names(cols) <- methods
      # make consensus model color black and shift rest of colors, such that lightest one is not needed
      cols[1:(which(names(cols)==name_consensus)-1)] <- cols[2:which(names(cols)==name_consensus)]
      cols[names(cols) == name_consensus] <- "#000000"

    } else if (palette == "Spectral"){
      cols_weekdays <- brewer.pal(name="Spectral", n=10)[c(1:4, 8:10)]
      names(cols_weekdays) <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
      
      if (length(intersect(methods, names(cols_weekdays))) == 0) {
        cols <- cols_weekdays[weekdays(as_date(methods))]
        names(cols) <- methods
      } else if (methods[num_est] == name_consensus) {
        cols <- c(cols_weekdays[methods[-num_est]], "#000000")
        names(cols)[num_est] <- name_consensus
        cols <- na.omit(cols[c(names(cols_weekdays), name_consensus)])
      } else {
        cols <- cols_weekdays
      }
      
      # make consensus model color black and shift rest of colors, such that lightest one is not needed
      cols[names(cols) == name_consensus] <- "#000000"
      
    } else {
      start_col <- 2
      cols <- brewer.pal(name=palette,n=num_est+start_col)[start_col+(1:num_est)]
      names(cols) <- methods
      # make consensus model color black
      cols[names(cols) == name_consensus] <- "#000000"
    }
  }
  
  #show_col(cols, ncol = num_est, labels = F)
  return(cols)
}

# function for saving pheatmap
save_pheatmap_pdf <- function(x, filename, width=9, height=5) {
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



plot_multiple_estimates <- function(estimates, legend_name, plot_title="",
                                    col_palette="Set1", name_consensus="consensus",
                                    include_CI=F, sort_numerically=F, long_time_frame = F) {
  
  # reshape data
  R_est  <- estimates %>%
    gather("variable", "value", 2:dim(estimates)[2]) %>%
    separate(variable,
             into = c("type", "model"),
             sep = "[.]",
             extra = "merge",
             remove = TRUE) %>%
    spread(type, value)

  if (sort_numerically) {
    R_est <- R_est %>%
      mutate_at(vars("model"), as.numeric) %>%
      arrange(date, model) %>%
      mutate_at(vars("model"), as.factor)
  }
  
  # plot
  R_plot <- ggplot(data = R_est, aes(x = date, y = R)) +
    geom_hline(aes(yintercept = 1)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size=18),
      axis.text=element_text(size=16),
      axis.title=element_text(size=18),
      legend.text=element_text(size=16),
      legend.title=element_text(size=18),
      axis.line = element_line(),
      axis.line.y.right = element_line(),
      axis.line.x.top = element_line(),
      legend.position = "bottom",
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major = element_line(),
      panel.grid.minor = element_blank()
    ) +
    ggtitle(plot_title)
  
  if (class(R_est$date) == "difftime") {
    R_plot <- R_plot +
      labs(x = "pub date - target date", y = "Rt estimate") +
      scale_x_continuous()
  }
  
  if (long_time_frame) {
    R_plot <- R_plot +
      geom_rect(data=NULL,aes(xmin=as_date("2020-03-01"), xmax=as_date("2020-03-31"),
                              ymin=-Inf, ymax=Inf), fill="grey", alpha=0.01) +
      geom_rect(data=NULL,aes(xmin=as_date("2021-06-11"), xmax=as_date("2021-07-09"),
                              ymin=-Inf, ymax=Inf), fill="grey", alpha=0.01)
  }
  
  col_values <- get_colors(methods = unique(R_est$model), col_palette, name_consensus = name_consensus)
  
  if (include_CI){
    R_plot <-  R_plot +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = .25) +
      scale_fill_manual(values=col_values, name=legend_name)
  } else {
    R_plot <-  R_plot +
      #geom_line(aes(group = model, color=model)) +
      geom_line(data=R_est[R_est$model!=name_consensus & !is.na(R_est$R),],
                aes(x = date, y = R, color=model), size = .5, na.rm = T) +
      geom_line(data=R_est[R_est$model==name_consensus,],
                aes(x = date, y = R, color=model), size = .8) +
      scale_color_manual(values=col_values, name=legend_name)
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

plot_for_comparison <- function(estimates, comp_methods, start_absdiff = "2020-04-01",
                                start_date = "2021-01-01", end_date = "2021-06-10",
                                legend_name="model", plot_title="",
                                col_palette="Set1",
                                include_consensus = T, name_consensus="consensus",
                                filenames = "_latest_plot.png",
                                include_CI=F, plot_diff_matrices=F, sort_numerically=F,
                                long_time_frame = F,
                                ylim_l=0, ylim_u=2.05,
                                verbose = T){
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
  
  
  if (!include_consensus) {
    estimates <- dplyr::select(estimates, !contains(name_consensus))
  }

  estimates_plot <- estimates %>%
    dplyr::filter(date >= start_date, date <= end_date)
  
  if (verbose){
    print("Date range for plotting:")
    print(min(estimates_plot$date))
    print(max(estimates_plot$date))
  }

  R_plot <- plot_multiple_estimates(estimates_plot, legend_name, plot_title = plot_title,
                                    col_palette = col_palette, name_consensus = name_consensus,
                                    include_CI = include_CI,
                                    sort_numerically = sort_numerically,
                                    long_time_frame = long_time_frame) +
    coord_cartesian(ylim = c(ylim_l-0.05, ylim_u+0.05), expand = FALSE)
  
  ggsave(R_plot, filename = paste0("Figures/estimates", filenames),  bg = "transparent",
         width = 13.1, height = 5.8)
  print(R_plot)
  
  if(plot_diff_matrices) {
    estimates_absdiff <- estimates %>%
      dplyr::filter(date >= start_absdiff, date <= end_date)
    estimates_absdiff <- estimates_absdiff[rowSums(is.na(estimates_absdiff)) == 0,]
    print("Date range for mean absolute difference calculation:")
    print(min(estimates_absdiff$date))
    print(max(estimates_absdiff$date))
    print("Number of days in MAD:")
    print(as.numeric(max(estimates_absdiff$date) - min(estimates_absdiff$date)))
    
    # adjust method names
    if (include_consensus){
      methods_ <- comp_methods
    } else {
      methods_ <- comp_methods[!comp_methods %in% c(name_consensus)]
    }
    for (i in c(" ", ",", "\\(", "\\)")) {methods_ <- gsub(i, ".", methods_)}
    
    n <- length(methods_)
    matr <- matrix(rep(rep(0,n), n), ncol=n)
    corr <- matrix(rep(rep(0,n), n), ncol=n)
    colnames(matr) <- rownames(matr) <- colnames(corr) <- rownames(corr) <- methods_
    
    for (method1 in methods_) {
      for (method2 in methods_){
        estimates_absdiff <- data.frame(estimates_absdiff)
        diff <- estimates_absdiff[,paste0("R.", method1)] - estimates_absdiff[,paste0("R.", method2)]
        matr[method1, method2] <- mean(as.vector(abs(diff)))
        #corr[method1, method2] <- cor(estimates_absdiff[,paste0("R.", method1)],
        #                              estimates_absdiff[,paste0("R.", method2)])
      }
    }
    
    df <- data.frame(matr)
    row.names(df) <- colnames(df) <- methods_

    if (sort_numerically){
      df <- df[order(as.numeric(row.names(df))), order(as.numeric(colnames(df)))]
    } else{
      df <- df[order(row.names(df)), order(colnames(df))]
    }
    
    mean_abs_diff <- pheatmap(df, color = viridis(100), breaks = seq(0,0.2,0.2/100),
                              border_color = NA, cellwidth = unit(2, "cm"), cellheight = unit(1, "cm"),
                              display_numbers = TRUE,
                              fontsize = 18, fontsize_number=22, number_color = "white",
                              angle_col = 315, cluster_rows = FALSE, cluster_cols = FALSE,
                              legend = FALSE)
    save_pheatmap_pdf(mean_abs_diff, paste0("Figures/mean_abs_difference", filenames))
    
    #correlations <- pheatmap(corr, color = viridis(100, direction = -1), breaks = seq(0.25,1,0.75/100),
    #                         border_color = NA, display_numbers = TRUE,
    #                         fontsize = 18, fontsize_number=22, number_color = "white",
    #                         angle_col = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE)
    #save_pheatmap_pdf(correlations, paste0("correlations", filenames))
  }
}



plot_real_time_estimates <- function(estimates,
                                     start_date = "2021-04-01", end_date = "2021-05-01",
                                     legend_name="weekday (pub date)", plot_title="",
                                     name_consensus="2021-07-16",
                                     filenames = "_latest_plot.pdf",
                                     ylim_l=0.5, ylim_u=1.5) {
  
  estimates <- estimates %>%
    dplyr::filter(date >= start_date, date <= end_date)
  
  # reshape data
  R_est  <- estimates %>%
    gather("variable", "value", 2:dim(estimates)[2]) %>%
    separate(variable,
             into = c("type", "model"),
             sep = "[.]",
             extra = "merge",
             remove = TRUE) %>%
    spread(type, value)
  
  R_est$weekday <- weekdays(as_date(R_est$model))
  
  # plot
  R_plot <- ggplot(data = R_est, aes(x = date, y = R), color = weekday) +
    geom_hline(aes(yintercept = 1)) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size=18),
      axis.text=element_text(size=16),
      axis.title=element_text(size=18),
      legend.text=element_text(size=16),
      legend.title=element_text(size=18),
      axis.line = element_line(),
      axis.line.y.right = element_line(),
      axis.line.x.top = element_line(),
      legend.position = "bottom",
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major = element_line(),
      panel.grid.minor = element_blank()
    ) +
    ggtitle(plot_title) +
    
    labs(x = NULL, y = "Rt estimate") +
    scale_x_date(limits = as.Date(c(min(estimates$date),max(estimates$date)-1)),
                 date_labels = "%b %d", expand = c(0,1),
                 breaks = date_breaks("1 week"))
  
  col_values <- get_colors(methods = c(unique(R_est$weekday), name_consensus),
                           "Spectral", name_consensus = name_consensus)
  
  R_plot <-  R_plot +
    geom_line(data=R_est[R_est$model!=name_consensus & !is.na(R_est$R),],
              aes(x = date, y = R, group=model, color=weekday), size = .5, na.rm = T) +
    geom_line(data=R_est[R_est$model==name_consensus,],
              aes(x = date, y = R, color=model), size = .8) +
    scale_color_manual(values=col_values, name=legend_name)
  
  R_plot <- R_plot +
    coord_cartesian(ylim = c(ylim_l-0.05, ylim_u+0.05), expand = FALSE)
  
  ggsave(R_plot, filename = paste0("Figures/estimates", filenames),  bg = "transparent",
         width = 13.1, height = 5.8)
  print(R_plot)
}



plot_weekday_effects <- function(estimates,
                                 legend_name="estimated on",
                                 plot_title="Mean estimates in week previous to pub date",
                                 min_lag = 1,
                                 max_lag = 7,
                                 col_palette="Spectral",
                                 filenames = "latest_plot.png",
                                 ylim_l=0.5, ylim_u=1.5,
                                 verbose = T){
  # reshape data
  R_est <- estimates %>%
    gather("pub_weekday", "mean_R", 2:dim(estimates)[2])
  
  wds <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
  
  # plot
  R_plot <- ggplot(data = R_est, aes(x = factor(target_weekday,
                                                levels = wds),
                                     y = mean_R,
                                     group = pub_weekday,
                                     color = pub_weekday)) +
    geom_line() +
    geom_hline(aes(yintercept = 1))
  
  R_plot <- R_plot +
    labs(x = "weekday (target date)", y = "Rt estimate",
         subtitle = paste(min_lag, "-", max_lag, "days previous to pub date")) +
    scale_x_discrete()
  
  R_plot <- R_plot +
    theme_minimal() +
    theme(
      plot.title = element_text(size=18),
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
    ggtitle(plot_title)
  
  col_values <- get_colors(methods = unique(R_est$pub_weekday), col_palette)
  
  R_plot <- R_plot +
    scale_color_manual(values=col_values, name=legend_name)

  R_plot <- R_plot +
    coord_cartesian(ylim = c(ylim_l, ylim_u), expand = FALSE)
  
  ggsave(R_plot, filename = paste0("Figures/estimates_realtime_target_date_influence/", filenames),  bg = "transparent",
         width = 13.1, height = 5.8)
  print(R_plot)
  
  
  # plot with inverted legend and x-axis
  R_plot <- ggplot(data = R_est, aes(x = factor(pub_weekday,
                                                levels = wds),
                                     y = mean_R,
                                     group = target_weekday,
                                     color = target_weekday)) +
    geom_line() +
    geom_hline(aes(yintercept = 1))
  
  R_plot <- R_plot +
    labs(x = "weekday (pub date)", y = "Rt estimate",
         subtitle = paste(min_lag, "-", max_lag, "days previous to publishing")) +
    scale_x_discrete()
  
  R_plot <- R_plot +
    theme_minimal() +
    theme(
      plot.title = element_text(size=18),
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
    ggtitle(plot_title)
  
  col_values <- get_colors(methods = unique(R_est$target_weekday), col_palette)
  
  R_plot <- R_plot +
    scale_color_manual(values=col_values, name="weekday (target date)")
  
  R_plot <- R_plot +
    coord_cartesian(ylim = c(ylim_l, ylim_u), expand = FALSE)
  
  ggsave(R_plot, filename = paste0("Figures/estimates_realtime_pub_date_influence/", filenames),  bg = "transparent",
         width = 13.1, height = 5.8)
  print(R_plot)
}



plot_CI_coverage_rates <- function(conf_level = "95"){
  methods <- c("Braunschweig", "epiforecasts", "ETHZ_sliding_window", "globalrt_7d",
               "ilmenau", "RKI_7day", "rtlive", "SDSC")
  
  CI_coverage <- read_csv(paste0("Rt_estimate_reconstruction/otherFiles/",
                                 conf_level, "_CI_coverage.csv")) %>%
    as.data.frame() %>%
    column_to_rownames("...1")
  
  coverage_data <- CI_coverage[, c(as.character(0:20), "min_lag")] %>%
    rownames_to_column("method") %>%
    dplyr::filter(method %in% methods) %>%
    mutate(method = plyr::mapvalues(method,
                                    c("Braunschweig", "ETHZ_sliding_window", "globalrt_7d", "ilmenau", "RKI_7day"),
                                    c("HZI",          "ETH",                 "globalrt",    "Ilmenau", "RKI"))) %>%
    arrange(method)
  
  coverage_data <- coverage_data %>%
    gather("variable", "value", 2:(dim(coverage_data)[2] - 1)) %>%
    mutate(variable = -1 * as.numeric(variable))
  
  coverage_plot <- ggplot() +
    theme_minimal() +
    theme(
      plot.title = element_text(size=18),
      axis.text=element_text(size=16),
      axis.title=element_text(size=18),
      legend.text=element_text(size=16),
      legend.title=element_text(size=18),
      axis.line = element_line(),
      axis.line.y.right = element_line(),
      axis.line.x.top = element_line(),
      legend.position = "bottom",
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major = element_line(),
      panel.grid.minor = element_blank()
    ) +
    #ggtitle(paste0(conf_level, "%-CI coverage rates")) +
    labs(x = "target date - pub date", y = "coverage rate") +
    coord_cartesian(xlim = c(-20, 0), ylim = c(-0.03, 1.03), expand = FALSE)
  
  methods_legend <- unique(coverage_data$method)
  col_values <- get_colors(methods = methods_legend, palette = "methods")
  line_types <- rep(1, length(methods_legend))
  names(line_types) <- methods_legend
  if (("rtlive" %in% methods_legend) &
      ("globalrt" %in% methods_legend) &
      (conf_level == "95")) line_types["rtlive"] <- 2
  
  coverage_plot <- coverage_plot + 
    geom_line(data=coverage_data,
              aes(x = variable,
                  y = value,
                  color = method,
                  linetype = method),
              size = .8, na.rm = T) +
    scale_color_manual(values=col_values, name="method") +
    scale_linetype_manual(values=line_types)
    
  coverage_plot <- coverage_plot +
    geom_hline(yintercept=0.95, linetype="dashed")
  
  ggsave(coverage_plot, filename = paste0("Figures/CI/", conf_level, "_coverage_rates.pdf"),
         bg = "transparent", width = 8, height = 5.8)
  print(coverage_plot)
}


plot_CI_widths <- function(conf_level = "95"){
  methods <- c("Braunschweig", "epiforecasts", "ETHZ_sliding_window", "globalrt_7d",
               "ilmenau", "RKI_7day", "rtlive", "SDSC")
  
  CI_width <- read_csv(paste0("Rt_estimate_reconstruction/otherFiles/", conf_level, "_CI_width.csv")) %>%
    as.data.frame() %>%
    column_to_rownames("...1")
  
  width_data <- CI_width[, c(as.character(0:20), "min_lag")] %>%
    rownames_to_column("method") %>%
    dplyr::filter(method %in% methods) %>%
    mutate(method = plyr::mapvalues(method,
                                    c("Braunschweig", "ETHZ_sliding_window", "globalrt_7d", "ilmenau", "RKI_7day"),
                                    c("HZI",          "ETH",                 "globalrt",    "Ilmenau", "RKI"))) %>%
    arrange(method)
  
  width_data <- width_data %>%
    gather("variable", "value", 2:(dim(width_data)[2] - 1)) %>%
    mutate(variable = -1 * as.numeric(variable))
  
  width_plot <- ggplot() +
    theme_minimal() +
    theme(
      plot.title = element_text(size=18),
      axis.text=element_text(size=16),
      axis.title=element_text(size=18),
      legend.text=element_text(size=16),
      legend.title=element_text(size=18),
      axis.line = element_line(),
      axis.line.y.right = element_line(),
      axis.line.x.top = element_line(),
      legend.position = "bottom",
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major = element_line(),
      panel.grid.minor = element_blank()
    ) +
    #ggtitle(paste0("Mean width of ", conf_level, "%-CI")) +
    labs(x = "target date - pub date", y = "width") +
    coord_cartesian(xlim = c(-20, 0), ylim = c(0.015, 0.9), expand = FALSE)
  
  methods_legend <- unique(width_data$method)
  col_values <- get_colors(methods = methods_legend, palette = "methods")
  
  width_plot <- width_plot + 
    geom_line(data=width_data,
              aes(x = variable,
                  y = value,
                  color = method),
              size = .8, na.rm = T) +
    #scale_y_log10(breaks = trans_breaks(identity, identity)) +
    scale_color_manual(values=col_values, name="method")

  ggsave(width_plot, filename = paste0("Figures/CI/", conf_level, "_CI_widths.pdf"),
         bg = "transparent", width = 8, height = 5.8)
  print(width_plot)
}
