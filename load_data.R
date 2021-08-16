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


load_RKI_data <- function(){
  
  # path to RKI incidence data
  repo_incid <- "https://raw.githubusercontent.com/robert-koch-institut/SARS-CoV-2-Nowcasting_und_-R-Schaetzung/main/Nowcast_R_aktuell.csv"
  
  # load incidence data
  incidence_raw <- read_csv(repo_incid)
  names(incidence_raw) <- c("Datum", "NeuErkr", "lb_NeuErkr", "ub_NeuErkr",
                            "NeuErkr_ma4", "lb_NeuErkr_ma4", "ub_NeuErkr_ma4",
                            "R_7Tage", "lb_R_7Tage", "ub_R_7Tage")

  # load Rt estimates published by the RKI
  path <- "reproductive_numbers/data-processed/RKI_7day/2021-06-04-RKI_7day.csv"
  RKI_7day <- read_csv(path)
  RKI_7day <- RKI_7day[RKI_7day$type=="point",][c("date", "value")]
  #path_rki4 <- "reproductive_numbers/data-processed/RKI_4day/2021-06-04-RKI_4day.csv"
  #RKI_4day <- read_csv(path_rki4)
  #RKI_4day <- RKI_4day[RKI_4day$type=="point",][c("date", "value")]
  
  # return data frame with incidence and published R estimates
  return(data.frame(dates=incidence_raw$Datum,
                    I=incidence_raw$NeuErkr,
                    R=c(rep(NA, min(RKI_7day$date)-min(incidence_raw$Datum)),
                        RKI_7day$value,
                        rep(NA, max(incidence_raw$Datum)-max(RKI_7day$date)))))
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
  
  # load Rt estimates published by the ETH
  path_R <- "reproductive_numbers/data-processed/ETHZ_sliding_window/2021-07-13-ETHZ_sliding_window.csv"
  ETH_3day <- read_csv(path_R)
  ETH_3day <- ETH_3day[ETH_3day$type=="point",][c("date", "value")]
  names(ETH_3day) <- c("date", "R")
  
  # return data frame with case data and published R estimates
  return(full_join(data.frame(caseData), ETH_3day, by="date"))
}


load_Ilmenau_data <- function(){
  
  # path to preprocessed RKI incidence data
  path_incid <- "Rt_estimate_reconstruction/incidence_data/data_ger_tot.qs"
  
  # load incidence data
  incidence <- qread(path_incid)
  incidence <- data.frame(dates=incidence$date, I=incidence$tot.cases, new.cases=incidence$new.cases)
  incidence$dates <- as_date(incidence$dates)
  
  # load Rt estimates published by TU Ilmenau
  path <- "reproductive_numbers/data-processed/ilmenau/2021-07-12-ilmenau.csv"
  Ilmenau_7day <- read_csv(path)
  Ilmenau_7day <- Ilmenau_7day[Ilmenau_7day$type=="point",][c("date", "value")]
  names(Ilmenau_7day) <- c("dates", "R")
  
  # return data frame with incidence and published R estimates
  return(full_join(data.frame(incidence), Ilmenau_7day, by="dates"))
}


load_epiforecasts_data <- function(){
  
  # path to preprocessed RKI incidence data
  #path_incid <- "Rt_estimate_reconstruction/incidence_data/RKI_COVID19.csv"
  
  # load incidence data
  #incidence <- read_csv(path_incid)
  #incidence <- data.frame(dates=incidence$date, I=incidence$tot.cases, new.cases=incidence$new.cases)
  #incidence$dates <- as_date(incidence$dates)
  
  # load Rt estimates from epiforecasts
  path <- "reproductive_numbers/data-processed/epiforecasts/2020-12-15-epiforecasts.csv"
  epiforecasts_1day <- read_csv(path)
  epiforecasts_1day <- epiforecasts_1day[epiforecasts_1day$type=="point"&
                                           epiforecasts_1day$location=="DE",][c("date", "value")]
  names(epiforecasts_1day) <- c("dates", "R")
  
  # return data frame with incidence and published R estimates
  return(epiforecasts_1day)
}