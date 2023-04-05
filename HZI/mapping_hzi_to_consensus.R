setwd("../..")
getwd()

source("Rt_estimate_reconstruction/load_data.R")

mapping <- read.csv("Rt_estimate_reconstruction/HZI/mapping_hzi_to_consensus.csv")
hzi_R_pub <- load_published_R_estimates("Braunschweig", pub_date = "2021-07-11")

get_mapped_R <- function(R_HZI){
  if (is.na(R_HZI)){
    return(NA)
  } else if (R_HZI == 0) {
    return(0)
  } else {
    return(mapping[which.min(abs(R_HZI - mapping$values_R_hzi)), "values_R_consensus"])
  }
}

estimate <- hzi_R_pub %>%
  mutate(R_calc = sapply(R_pub, get_mapped_R),
         lower = sapply(lower, get_mapped_R),
         upper = sapply(upper, get_mapped_R)) %>%
  na.omit()

write.csv(estimate[,c("date", "R_calc", "lower", "upper")],
          "Rt_estimate_reconstruction/HZI/R_adjInputWindowGTD_2021-07-11.csv",
          row.names = F)
