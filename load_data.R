library(readr)
library(stringr)
library(dplyr)
library(tidyverse)
library(fitdistrplus)
library(qs)
library(lubridate)


# load R estimates from git repo "reproductive_numbers"
load_published_R_estimates <- function(source, start=as.Date("2019-12-28"), end=Sys.Date(), location="DE"){
  
  # define path where Rt estimates are located
  path <- paste0("reproductive_numbers/data-processed/", source, "/")
  
  tryCatch({
    # find most recent estimates
    file <- max(list.files(path, pattern = paste0("\\d{4}-\\d{2}-\\d{2}-", source, ".csv")))
    
    # load estimates
    R_est <- read_csv(paste0(path, file), col_types = list(date = col_date()))
    
  }, error=function(e) {
    print("Unmatched source, choose from:")
    print(list.dirs("reproductive_numbers/data-processed/", full.names = F))
  })
  
  R_est <- R_est[(R_est$type=="point") & (R_est$location==location), ][c("date", "value")]
  names(R_est) <- c("date", "R_pub")
  
  # return available R estimates for time between start and end
  dates <- data.frame(date=seq(start, end, by ="day"))
  data <- merge(R_est, dates, by.x = 'date', by.y = 'date', all.x = FALSE, all.y = TRUE)
  return(data)
}


# call function to load incidence data depending on method/source/paper
load_incidence_data <- function(method, location="DE"){
  
  # depending on method call functions
  if (method == "RKI") {
    if (location == "DE") {
      data <- load_RKI_data()
    } else {
      print("RKI incidence data only available for Germany, pass location = 'DE'.")
    }
    
  } else if (method == "ETHZ_sliding_window") {
    # pass country/region value pair
    if (location == "DE") {
      data <- load_ETH_deconvolvedCountryData()
    } else if (location == "AT") {
      data <- load_ETH_deconvolvedCountryData(country="Austria", region="AUT")
    } else if (location == "CH") {
      data <- load_ETH_deconvolvedCountryData(country="Switzerland", region="CHE")
    } else {
      print("Location not included in analysis, choose from [DE, AT, CH].")
    }
    
  } else if (method == "ilmenau") {
    if (location == "DE") {
      data <- load_Ilmenau_data()
    } else {
      print("Ilmenau incidence data not yet processed for locations other than Germany, pass location = 'DE'.")
    }
    
  } else if (method == "AGES") {
    if (location == "AT") {
      data <- load_AGES_data()
    } else {
      print("AGES incidence data only available for Austria, pass location = 'AT'.")
    }
    
  } else {
    print("Method unknown, choose from [RKI, ETHZ_sliding_window, ilmenau, AGES].")
  }
  
  return(data)
}



#########
# incidence data from different sources (incl. preprocessing from different papers)
#########

load_RKI_data <- function(){
  
  # path to RKI incidence data
  repo_incid <- "https://raw.githubusercontent.com/robert-koch-institut/SARS-CoV-2-Nowcasting_und_-R-Schaetzung/main/Nowcast_R_aktuell.csv"
  
  # load incidence data
  data <- read_csv(repo_incid)
  names(data) <- c("Datum", "NeuErkr", "lb_NeuErkr", "ub_NeuErkr",
                   "NeuErkr_ma4", "lb_NeuErkr_ma4", "ub_NeuErkr_ma4",
                   "R_7Tage", "lb_R_7Tage", "ub_R_7Tage")
  data <- data.frame(date=data$Datum, I=data$NeuErkr)
  
  return(data)
}


load_ETH_deconvolvedCountryData <- function(country="Germany", region="DEU",
                                            force_deconvolution = FALSE){
  
  countryDataPath <- file.path("Rt_estimate_reconstruction/ETH/data/countryData",
                               str_c(region, "-DeconvolutedData.rds"))

  # Deconvolution if necessary or wanted
  if (force_deconvolution | !file.exists(countryDataPath)){
    ETH_deconvolution(country=country, region=region)
  }
  
  # load deconvoluted data
  deconvolvedCountryData <- readRDS(countryDataPath)
  
  return(deconvolvedCountryData)
}


load_Ilmenau_data <- function(){
  
  # path to preprocessed RKI incidence data
  path_incid <- "Rt_estimate_reconstruction/incidence_data/data_ger_tot.qs"
  
  # load incidence data
  data <- qread(path_incid)
  data <- data.frame(date=data$date, I=data$new.cases)
  data$date <- as_date(data$date)
  
  return(data)
}


load_AGES_data <- function(){
  
  # path to EMS incidence data
  link_incid <- "https://covid19-dashboard.ages.at/data/CovidFaelle_Timeline.csv"

  # load data
  data_raw <- read.csv(link_incid, sep = ";", dec = ",")
  
  # sum over days
  data <- aggregate(data_raw$AnzahlFaelle, by=list(data_raw$Time), FUN=sum)
  names(data) <- c("date", "I")
  
  # convert date column
  data <- mutate(data, date = as.Date(date, format = "%d.%m.%Y"))
  
  # sort values by date
  data <- data[order(data$date),]
  rownames(data) <- 1:nrow(data)
  
  return(data)
}



#########
# additional data
#########

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



#########
# functions for preprocessing
#########

ETH_deconvolution <- function(country="Germany", region="DEU"){
  
  source("Rt_estimate_reconstruction/ETH/otherScripts/1_utils_getRawData.R")
  source("Rt_estimate_reconstruction/ETH/otherScripts/utils.R")
  
  basePath <- "Rt_estimate_reconstruction/ETH/data/countryData"
  if (!dir.exists(basePath)) {
    dir.create(basePath)
  }
  
  # fetch stringency data
  if (region == "CHE") {
    stringencyData <- read_csv(
      "https://raw.githubusercontent.com/KOF-ch/economic-monitoring/master/data/ch.kof.stringency.csv",
      col_types = cols(
        time = col_date(format = ""),
        geo = col_character(),
        variable = col_character(),
        value = col_double()
      )) %>%
      filter(
        variable == "stringency") %>%
      dplyr::transmute(
        date = time,
        countryIso3 = "CHE",
        region = recode(toupper(geo), "CH" = "CHE"),
        source = "KOF",
        StringencyIndex = value
      )
  } else {
    stringencyData <- getDataOxfordStringency(countries = region, tReload = 300) %>%
      mutate(source = "BSG Covidtracker")
  }
  
  stringencyDataPath <- file.path(basePath, str_c(region, "-OxCGRT.rds"))
  
  if (file.exists(stringencyDataPath)) {
    stringencyDataOld <- readRDS(stringencyDataPath)
    # if new data is null, keep old data (can happen because of error in reading new data)
    if (is.null(stringencyData)) {
      stringencyData <- stringencyDataOld
    }
    stringencyUnchanged <- all.equal(stringencyData, stringencyDataOld)
  } else {
    stringencyUnchanged <- FALSE
  }
  
  if (!isTRUE(stringencyUnchanged)) {
    saveRDS(stringencyData, file = stringencyDataPath)
  }
  
  stringencyIndex <- stringencyData %>%
    dplyr::transmute(
      date,
      countryIso3,
      region,
      data_type = "Stringency Index",
      source = source,
      value = StringencyIndex
    ) %>%
    filter(!is.na(value)) %>%
    arrange(countryIso3, region, date)
  
  # Fetch Country Data
  countryData <- getCountryData(
    region,
    #tempFile = "Rt_estimate_reconstruction/ETH/data/temp/ECDCdata.csv",
    HMDtemp = "Rt_estimate_reconstruction/ETH/data/temp/HMDdata.csv",
    tReload = 300)
  
  countryData <- countryData %>%
    bind_rows(
      mutate(stringencyIndex,
             date_type = if_else(
               region %in% c("CHE", "DEU", "HKG"), "report_plotting", "report"))
    )
  
  # save updated data
  countryDataPath <- file.path(basePath, str_c(region, "-Data.rds"))
  saveRDS(countryData, file = countryDataPath)
  
  # get Infection Incidence
  # load functions
  source("Rt_estimate_reconstruction/ETH/otherScripts/2_utils_getInfectionIncidence.R")
  # load parameter
  source("Rt_estimate_reconstruction/ETH/otherScripts/2_params_InfectionIncidencePars.R")
  # load empirical delays
  delays_data_path <- "Rt_estimate_reconstruction/ETH/data/all_delays.csv"
  delays_onset_to_count <- read_csv(delays_data_path,
                                    col_types = cols(
                                      data_type = col_character(),
                                      onset_date = col_date(format = ""),
                                      count_date = col_date(format = ""),
                                      delay = col_number()))
  
  # constant delay distribution
  constant_delay_distributions <- list()
  for (type_i in unique(names(shape_onset_to_count))) {
    m <- get_vector_constant_waiting_time_distr(
      shape_incubation,
      scale_incubation,
      shape_onset_to_count[[type_i]],
      scale_onset_to_count[[type_i]])
    
    constant_delay_distributions <- c(constant_delay_distributions, list(m))
  }
  names(constant_delay_distributions) <- unique(names(shape_onset_to_count))
  
  constant_delay_symptom_to_report_distributions <- list()
  for (type_i in unique(names(shape_onset_to_count))) {
    m <- get_vector_constant_waiting_time_distr(
      0,
      0,
      shape_onset_to_count[[type_i]],
      scale_onset_to_count[[type_i]])
    
    constant_delay_symptom_to_report_distributions <- c(constant_delay_symptom_to_report_distributions, list(m))
  }
  names(constant_delay_symptom_to_report_distributions) <- paste0('Onset to ',  unique(names(shape_onset_to_count)))
  
  constant_delay_distributions <- c(constant_delay_distributions, constant_delay_symptom_to_report_distributions)
  
  # filter out regions with too few cases for estimation
  countryData <- countryData %>%
    filterRegions(thresholdConfirmedCases = 500)
  # remove Oxford Stringenxy Index for Re calculation
  countryData <- countryData %>%
    filter(data_type != "Stringency Index")
  
  # filter out data_types with 0 total cases
  data_type0 <- countryData %>%
    group_by(data_type) %>%
    summarize(total = sum(value), .groups = "drop") %>%
    filter(total == 0) %>%
    .$data_type
  
  countryData <- filter(countryData, !(data_type %in% data_type0))
  
  # data filtering: ignore deaths
  countryData <- countryData %>%
    filter(data_type != "Deaths")
  
  countryData <- countryData %>%
    mutate(
      data_type = fct_drop(data_type)
    )
  
  right_truncation <- list()
  if (region %in% c("CHE", "LIE", "DEU", "HKG")) {
    right_truncation[["Confirmed cases"]] <- 0
    right_truncation[["Confirmed cases / tests"]] <- 0
    right_truncation[["Hospitalized patients"]] <- 0
    #right_truncation[["Deaths"]] <- 0
  } else {
    right_truncation["Confirmed cases"] <- 3
    right_truncation["Confirmed cases / tests"] <- 3
    right_truncation["Hospitalized patients"] <- 3
    #right_truncation["Deaths"] <- 3
  }
  
  right_truncate <- function(df, data_type, right_truncation) {
    dplyr::filter(df, date <= (max(date) - right_truncation[[unique(data_type)]]))
  }
  
  countryData$country <- country
  
  countryData <- countryData %>%
    group_by(country, region, source, data_type) %>%
    right_truncate(data_type, right_truncation) %>%
    dplyr::select(-countryIso3) %>%
    ungroup()
  
  # Deconvolution
  deconvolvedData <- list()
  
  deconvolvedData[[1]] <- get_all_infection_incidence(
    countryData,
    constant_delay_distributions = constant_delay_distributions,
    onset_to_count_empirical_delays = delays_onset_to_count,
    data_types = c("Confirmed cases",
                   "Hospitalized patients"),
    #"Deaths"),
    n_bootstrap = 100,
    verbose = FALSE)
  
  if (region %in% c("CHE")) {
    countryDataTests <- countryData %>%
      filter(data_type == "Confirmed cases / tests")
    
    deconvolvedData[[2]] <- get_all_infection_incidence(
      countryDataTests,
      constant_delay_distributions = constant_delay_distributions,
      onset_to_count_empirical_delays = delays_onset_to_count,
      data_types = c("Confirmed cases / tests"),
      n_bootstrap = 100,
      verbose = FALSE)
  }
  
  deconvolvedCountryData <- bind_rows(deconvolvedData)
  countryDataPath <- file.path(basePath, str_c(region, "-DeconvolutedData.rds"))
  if (dim(deconvolvedCountryData)[1] == 0) {
    print("no data remaining")
  } else {
    saveRDS(deconvolvedCountryData, file = countryDataPath)
  }
}
