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

  # define path where Rt estimates are located
  path <- "reproductive_numbers/data-processed/RKI_7day/"
  
  # find most recent estimates
  file <- max(list.files(path, pattern = "\\d{4}-\\d{2}-\\d{2}-RKI_7day.csv"))
  
  # load Rt estimates
  RKI_7day <- read_csv(paste0(path, file))
  RKI_7day <- RKI_7day[RKI_7day$type=="point",][c("date", "value")]
  #path_rki4 <- "reproductive_numbers/data-processed/RKI_4day/"
  #file_rki4 <- max(list.files(path_rki4, pattern = "\\d{4}-\\d{2}-\\d{2}-RKI_4day.csv"))
  #RKI_4day <- read_csv(paste0(path_rki4, file_rki4))
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
  
  
  # define path where Rt estimates are located
  path <- "reproductive_numbers/data-processed/ETHZ_sliding_window/"
  
  # find most recent estimates
  file <- max(list.files(path, pattern = "\\d{4}-\\d{2}-\\d{2}-ETHZ_sliding_window.csv"))
  
  # load Rt estimates
  ETH_3day <- read_csv(paste0(path, file))
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
  
  # define path where Rt estimates are located
  path <- "reproductive_numbers/data-processed/ilmenau/"
  
  # find most recent estimates
  file <- max(list.files(path, pattern = "\\d{4}-\\d{2}-\\d{2}-ilmenau"))
  
  # load Rt estimates
  Ilmenau_7day <- read_csv(paste0(path, file))
  Ilmenau_7day <- Ilmenau_7day[Ilmenau_7day$type=="point",][c("date", "value")]
  names(Ilmenau_7day) <- c("dates", "R")
  
  # return data frame with incidence and published R estimates
  return(full_join(data.frame(incidence), Ilmenau_7day, by="dates"))
}


load_AGES_data <- function(){
  
  # path to EMS incidence data
  link_incid <- "https://covid19-dashboard.ages.at/data/CovidFaelle_Timeline.csv"

  # load data
  data_raw <- read.csv(link_incid, sep = ";")

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