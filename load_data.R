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


load_published_R_estimates <- function(source, incid){
  # define path where Rt estimates are located
  path <- paste0("reproductive_numbers/data-processed/", source, "/")
  
  # find most recent estimates
  file <- max(list.files(path, pattern = paste0("\\d{4}-\\d{2}-\\d{2}-", source, ".csv")))
  
  # load estimates
  R_est <- read_csv(paste0(path, file), col_types = list(date = col_date()))
  R_est <- R_est[R_est$type=="point",][c("date", "value")]
  names(R_est) <- c("date", "R_pub")
  
  # merge incidence data and published R estimates
  data <- full_join(incid, R_est, by="date")
  
  return(data)
}


load_RKI_data <- function(include_R = TRUE){
  
  # path to RKI incidence data
  repo_incid <- "https://raw.githubusercontent.com/robert-koch-institut/SARS-CoV-2-Nowcasting_und_-R-Schaetzung/main/Nowcast_R_aktuell.csv"
  
  # load incidence data
  incid <- read_csv(repo_incid)
  names(incid) <- c("Datum", "NeuErkr", "lb_NeuErkr", "ub_NeuErkr",
                    "NeuErkr_ma4", "lb_NeuErkr_ma4", "ub_NeuErkr_ma4",
                    "R_7Tage", "lb_R_7Tage", "ub_R_7Tage")
  incid <- data.frame(date=incid$Datum, I=incid$NeuErkr)
  
  # load Rt estimates
  if (include_R){
    data <- load_published_R_estimates("RKI_7day", incid)
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
  incid <- qread(path_incid)
  incid <- data.frame(date=incid$date, I=incid$new.cases)
  incid$date <- as_date(incid$date)
  
  # load Rt estimates
  if(include_R){
    data <- load_published_R_estimates("ilmenau", incid)
  }
  
  return(data)
}


load_AGES_data <- function(){
  
  # path to EMS incidence data
  link_incid <- "https://covid19-dashboard.ages.at/data/CovidFaelle_Timeline.csv"

  # load data
  data_raw <- read.csv(link_incid, sep = ";", dec = ",")

  incidence <- aggregate(data_raw$AnzahlFaelle, by=list(data_raw$Time), FUN=sum)
  names(incidence) <- c("dates", "I")
  incidence <- mutate(incidence, dates = as.Date(dates, format = "%d.%m.%Y"))

  # define path where Rt estimates are located
  path <- "reproductive_numbers/data-processed/AGES/"
  
  # find most recent estimates
  file <- max(list.files(path, pattern = "\\d{4}-\\d{2}-\\d{2}-AGES"))
  
  # load Rt estimates
  AGES_13day <- read_csv(paste0(path, file))
  AGES_13day <- AGES_13day[AGES_13day$type=="point",][c("date", "value")]
  names(AGES_13day) <- c("dates", "R")
  
  # return data frame with incidence and published R estimates
  data <- full_join(data.frame(incidence), AGES_13day, by="dates")
  data <- data[order(data$dates),]
  rownames(data) <- 1:nrow(data)
  return(data)
}


load_epiforecasts_data <- function(){
  
  # define path where Rt estimates are located
  path <- "reproductive_numbers/data-processed/epiforecasts/"
  
  # find most recent estimates
  file <- max(list.files(path, pattern = "\\d{4}-\\d{2}-\\d{2}-epiforecasts.csv"))
  
  # load Rt estimates
  epiforecasts_1day <- read_csv(paste0(path, file))
  epiforecasts_1day <- epiforecasts_1day[epiforecasts_1day$type=="point"&
                                           epiforecasts_1day$location=="DE",][c("date", "value")]
  names(epiforecasts_1day) <- c("dates", "R")
  
  # return data frame with incidence and published R estimates
  return(epiforecasts_1day)
}