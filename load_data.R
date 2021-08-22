library(readr)
library(dplyr)
library(qs)
library(lubridate)


download_OWID_data <- function(){
  
  # path to OWID incidence data
  repo_incid <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv"
  
  # load data
  data_raw <- read_csv(repo_incid)
  
  # extract relevant countries
  owid_austria <- data_raw[data_raw$location=="Austria",]
  owid_germany <- data_raw[data_raw$location=="Germany",]
  owid_switzerland <- data_raw[data_raw$location=="Switzerland",]
  
  # write data as csv
  write.csv(owid_austria, "Rt_estimate_reconstruction/incidence_data/OWID_Austria.csv")
  write.csv(owid_germany, "Rt_estimate_reconstruction/incidence_data/OWID_Germany.csv")
  write.csv(owid_switzerland, "Rt_estimate_reconstruction/incidence_data/OWID_Switzerland.csv")
}


load_published_R_estimates <- function(source, incid, location="DE"){
  # define path where Rt estimates are located
  path <- paste0("reproductive_numbers/data-processed/", source, "/")
  
  # find most recent estimates
  file <- max(list.files(path, pattern = paste0("\\d{4}-\\d{2}-\\d{2}-", source, ".csv")))
  
  # load estimates
  R_est <- read_csv(paste0(path, file), col_types = list(date = col_date()))
  R_est <- R_est[(R_est$type=="point") & (R_est$location==location), ][c("date", "value")]
  names(R_est) <- c("date", "R_pub")
  
  # merge incidence data and published R estimates
  data <- full_join(incid, R_est, by="date")
  
  return(data)
}


load_RKI_data <- function(include_R = TRUE){
  
  # path to RKI incidence data
  repo_incid <- "https://raw.githubusercontent.com/robert-koch-institut/SARS-CoV-2-Nowcasting_und_-R-Schaetzung/main/Nowcast_R_aktuell.csv"
  
  # load incidence data
  data <- read_csv(repo_incid)
  names(data) <- c("Datum", "NeuErkr", "lb_NeuErkr", "ub_NeuErkr",
                   "NeuErkr_ma4", "lb_NeuErkr_ma4", "ub_NeuErkr_ma4",
                   "R_7Tage", "lb_R_7Tage", "ub_R_7Tage")
  data <- data.frame(date=data$Datum, I=data$NeuErkr)
  
  # load Rt estimates
  if (include_R){
    data <- load_published_R_estimates("RKI_7day", data)
  }
  
  return(data)
}


load_ETH_data <- function(){
  
  # load preprocessed/deconvoluted data (got those from ETH repo)
  countryData <- qread("Rt_estimate_reconstruction/incidence_data/countryData.qs")
  
  # select only time series of confirmed cases (ignore deaths)
  caseData <- countryData[countryData$data_type=="Confirmed cases",]
  
  # assume all cases arose from local infection events
  caseData$local_infection <- rep(T, nrow(caseData))
  
  # use deconvoluted ts instead of reported values
  caseData$value <- caseData$deconvoluted
  
  # only keep data for which the deconvoluted ts is known
  caseData <- caseData[!is.na(caseData$value),]
  
  data <- load_published_R_estimates("ETHZ_sliding_window", caseData)
  
  # return data frame with case data and published R estimates
  return(data)
}


load_Ilmenau_data <- function(include_R = TRUE){
  
  # path to preprocessed RKI incidence data
  path_incid <- "Rt_estimate_reconstruction/incidence_data/data_ger_tot.qs"
  
  # load incidence data
  data <- qread(path_incid)
  data <- data.frame(date=data$date, I=data$new.cases)
  data$date <- as_date(data$date)
  
  # load Rt estimates
  if(include_R){
    data <- load_published_R_estimates("ilmenau", data)
  }
  
  return(data)
}


load_AGES_data <- function(include_R = TRUE){
  
  # path to EMS incidence data
  link_incid <- "https://covid19-dashboard.ages.at/data/CovidFaelle_Timeline.csv"

  # load data
  data_raw <- read.csv(link_incid, sep = ";", dec = ",")

  data <- aggregate(data_raw$AnzahlFaelle, by=list(data_raw$Time), FUN=sum)
  names(data) <- c("date", "I")
  data <- mutate(data, date = as.Date(date, format = "%d.%m.%Y"))
  
  # load Rt estimates
  if (include_R){
    data <- load_published_R_estimates("AGES", data, location="AT")
  }

  data <- data[order(data$date),]
  rownames(data) <- 1:nrow(data)
  return(data)
}

