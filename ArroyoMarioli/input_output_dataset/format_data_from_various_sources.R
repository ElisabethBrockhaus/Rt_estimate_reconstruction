library(readr)

setwd("D:/EllasDaten/Uni/Wirtschaftsingenieurwesen/6Semester/Bachelorarbeit/Code")
source("Rt_estimate_reconstruction/load_data.R")

# load JHU data (has correct format)
target_format <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/time_series_covid19_confirmed_global.csv")

locations <- c("DE", "AT", "CH")
all <- c("Germany", "Austria", "Switzerland")
names(locations) <- all


# function to load and reformat input data from different sources
load_data_for_globalrt <- function(method, countries=c("Germany")){
  
  # load data for given countries
  if (method == "ETH"){
    # get non-deconvolved data
    for (country in countries){
      country_data <- load_incidence_data(method = "ETHZ_sliding_window",
                                          location = locations[country],
                                          deconvolved = FALSE)
      names(country_data)[2] <- country
      joined_data <- if (!exists("joined_data")) country_data else full_join(joined_data, country_data, by = "date")
    }
  } else if ((method == "rtlive") & (countries == c("Germany"))) {
    country_data <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid.csv")
    names(country_data) <- c("date", "Germany")
    joined_data <- if (!exists("joined_data")) country_data else full_join(joined_data, country_data, by = "date")
  } else {
    for (country in countries){
      country_data <- load_incidence_data(method = method, location = locations[country])
      names(country_data)[2] <- country
      joined_data <- if (!exists("joined_data")) country_data else full_join(joined_data, country_data, by = "date")
    }
  }
  
  joined_data <- joined_data[joined_data$date < "2021-10-01",]
  
  # save data in format for epiforecast estimation
  write_csv(joined_data, paste0("Rt_estimate_reconstruction/epiforecasts/input_data/incidence_data_", method, ".csv"))
  
  # make date column the index
  data <- column_to_rownames(joined_data, var = "date")
  
  # bring date into same format as in JHU data
  rownames(data) <- format(as_date(rownames(data)), "%m/%d/%y")
  rownames(data) <- gsub("(0(?=[0-9]+))", "", rownames(data), perl=TRUE)
  
  # new cases to total cases
  data <- apply(data, MARGIN = 2, FUN = cumsum)

  # transpose data
  data <- t(data)
  
  # add columns at the beginning
  data <- cbind.data.frame(target_format[target_format$`Country/Region` %in% countries ,1:5], data)
  names(data)[1] <- NA
  rownames(data) <- NULL
  
  return(data)
}

# load, format and save data from RKI (used for main analysis)
#data_rki <- load_data_for_globalrt(method = "RKI", countries=c("Germany"))
write_csv(load_data_for_globalrt(method = "RKI", countries=c("Germany")),
          "Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/time_series_covid19_confirmed_global_RKI.csv")

write_csv(load_data_for_globalrt(method = "rtlive", countries=c("Germany")),
          "Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/time_series_covid19_confirmed_global_rtlive.csv")


#####################
# different sources #
#####################
# load data from ETHZ_sliding_window, ilmenau, AGES, SDSC
#data_eth <- load_data_for_globalrt(method = "ETH", countries=c("Germany")) # TODO: make usable for Austria
#data_ilmenau <- load_data_for_globalrt(method = "ilmenau", countries=c("Germany"))
#data_ages <- load_data_for_globalrt(method = "AGES", countries=c("Austria"))
#data_sdsc <- load_data_for_globalrt(method = "sdsc", countries=all)
#data_epiforecasts <- load_data_for_globalrt(method = "epiforecasts", countries=all)

# compare input data
#plot(t(target_format[1, 46:595]), type="l", col="darkgreen")
#lines(t(data_rki[1, 6:555]))
#lines(t(data_eth[1,27:576]), col="darkblue")
#lines(t(data_sdsc[1,46:595]), col="blue")
#lines(t(data_epiforecasts[1,65:614]), col="lightblue")
#legend(x="topleft", legend=c("globalrt", "RKI", "ETH", "SDSC", "epiforecasts"), lty=1, col=c("darkgreen", "black", "darkblue", "blue", "lightblue"))

# save data for Germany as .csv
#for (method in c("ETH", "ilmenau", "sdsc", "epiforecasts")) {
#  data <- load_data_for_globalrt(method = method)
#  write_csv(data, paste0("Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/time_series_covid19_confirmed_global_", method, ".csv"))
#  print(paste0("done for ", method))
#}


#####################
# compare estimates #
#####################
#estimated_R_org <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/estimates/estimated_R.csv")
#estimated_R_sdsc <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/estimates/estimated_R_sdsc.csv")

#estimate_org_DEU <- estimated_R_org[estimated_R_org$`Country/Region` == "Germany",]
#estimate_sdsc_DEU <- estimated_R_sdsc[estimated_R_sdsc$`Country/Region` == "Germany",]

#plot(estimate_org_DEU$Date, estimate_org_DEU$R, type="l")
#lines(estimate_sdsc_DEU$Date, estimate_sdsc_DEU$R, col="red")
