if (interactive()) {
  library(here)
  library(tidyverse)
} else {
  suppressPackageStartupMessages({
    library(here)
    library(tidyverse)
  })
}

outDir <- here::here("ETH", "data")

#linelist_codes <- c("CHE", "DEU", "HKG")
linelist_codes <- c("DEU")

# EB: source format_DEU_data.R only
#for(code in linelist_codes) {
#  source(here::here("ETH", "otherScripts", paste0("format_",code,"_data.R")))
#}
source("Rt_estimate_reconstruction/ETH/otherScripts/format_DEU_data.R")

delay_data <- c()
for (code in linelist_codes) {
  country_delay_data <- read_csv(here::here("Rt_estimate_reconstruction/ETH", "data", code, paste0(code, "_data_delays.csv")))
  delay_data <- c(delay_data, list(country_delay_data))
}

all_delays <- bind_rows(delay_data)

outDir <- here::here("Rt_estimate_reconstruction/ETH", "data")
write_csv(all_delays, path = file.path(outDir, "all_delays.csv"))
