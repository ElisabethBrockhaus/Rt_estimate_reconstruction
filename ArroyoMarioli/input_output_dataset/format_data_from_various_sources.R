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
    country_data <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_07_10.csv")
    names(country_data) <- c("date", "Germany_rtlive")
    country_data <- country_data %>%
      mutate(Germany_rtlive_rolling = zoo::rollapplyr(Germany_rtlive, width = 7, FUN = sum, partial = TRUE))
    #country_data$Germany_rtlive_rolling[1:6] <- NA
    joined_data <- if (!exists("joined_data")) country_data else full_join(joined_data, country_data, by = "date")

  } else if ((method == "rtlive_preprocessed") & (countries == c("Germany"))) {
    # TODO: delete this case (contains mistakes)
    country_data <- read_csv("Rt_estimate_reconstruction/incidence_data/rtlive_incid_21_07_10.csv")
    # calculate moving sum (7) over incidence for globalrt estimation with preprocessing substituting window size
    country_data$I <- as.numeric(stats::filter(country_data, rep(7, 7), method = "convolution", sides = 1)[,2])
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

  # remove NA
  data <- na.omit(data)
  
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

rtlive_data <- load_data_for_globalrt(method = "rtlive", countries=c("Germany"))
rtlive_data[1,2] <- "Germany_rtlive"
rtlive_data[2,2] <- "Germany_rtlive_rolling"

write_csv(rtlive_data, "Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/time_series_covid19_confirmed_global_final.csv")


# load, format and save data from RKI (used for main analysis)
#data_rki <- load_data_for_globalrt(method = "RKI", countries=c("Germany"))
#write_csv(load_data_for_globalrt(method = "RKI", countries=c("Germany")),
#          "Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/time_series_covid19_confirmed_global_RKI.csv")

#write_csv(load_data_for_globalrt(method = "rtlive", countries=c("Germany")),
#          "Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/time_series_covid19_confirmed_global_rtlive.csv")

#write_csv(load_data_for_globalrt(method = "rtlive_preprocessed", countries=c("Germany")),
#          "Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/time_series_covid19_confirmed_global_rtlive_preprocessed.csv")

