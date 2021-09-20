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
  if (method == "ETHZ_sliding_window"){
    # ask for non-deconvolved data
    for (country in countries){
      country_data <- load_incidence_data(method = method, location = locations[country], deconvolved = FALSE)
      names(country_data)[2] <- country
      joined_data <- if (!exists("joined_data")) country_data else full_join(joined_data, country_data, by = "date")
    }
  } else {
    for (country in countries){
      country_data <- load_incidence_data(method = method, location = locations[country])
      names(country_data)[2] <- country
      joined_data <- if (!exists("joined_data")) country_data else full_join(joined_data, country_data, by = "date")
    }
  }
  
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
  rownames(data) <- NULL
  
  # save data as .csv
  #write_csv(data, paste0("Rt_estimate_reconstruction/ArroyoMarioli/input_output_dataset/time_series_covid19_confirmed_global_", method, ".csv"))
  return(data)
}

# load data from RKI, ETHZ_sliding_window, ilmenau, AGES, SDSC
data_rki <- load_data_for_globalrt(method = "RKI", countries=c("Germany"))
data_eth <- load_data_for_globalrt(method = "ETHZ_sliding_window", countries=c("Germany")) # TODO: make usable for Austria
data_ilmenau <- load_data_for_globalrt(method = "ilmenau", countries=c("Germany"))
data_ages <- load_data_for_globalrt(method = "AGES", countries=c("Austria"))
data_scsc <- load_data_for_globalrt(method = "sdsc", countries=all)



#####################
# compare estimates #
#####################
estimated_R_org <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/estimates/estimated_R.csv")
estimated_R_sdsc <- read_csv("Rt_estimate_reconstruction/ArroyoMarioli/estimates/estimated_R_sdsc.csv")

estimate_org_DEU <- estimated_R_org[estimated_R_org$`Country/Region` == "Germany",]
estimate_sdsc_DEU <- estimated_R_sdsc[estimated_R_sdsc$`Country/Region` == "Germany",]

plot(estimate_org_DEU$Date, estimate_org_DEU$R, type="l")
lines(estimate_sdsc_DEU$Date, estimate_sdsc_DEU$R, col="red")
