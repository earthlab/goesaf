
# Load and process ERA5 climate data for fire events ----------------------

library(tidyverse)
library(lubridate)
library(pbapply)

era5_dir <- file.path('data', 'era5', 'raw')
dir.create(era5_dir, showWarnings = TRUE, recursive = TRUE)
tar_dest <- file.path(era5_dir, 'wxdatahourly.tar.gz')
download.file('https://earthlab-mjoseph.s3-us-west-2.amazonaws.com/wxdatahourly.tar.gz', 
              destfile = tar_dest)
untar(tar_dest, exdir = era5_dir)

column_headers <- c('year', 'day_of_year', 'hour', 'temp_c', 'precipitation_m', 
                    'solar_radiation_w_per_m2', 'wind_speed_m_per_s', 
                    'vpd_kPa', 'ten_hr_dead_fuel_moisture_pct')
era5_csv_paths <- list.files(path = era5_dir, pattern = '.csv$', 
                             full.names = TRUE)
dir.create(file.path('data', 'era5', 'clean'), showWarnings = FALSE)

# load one and get it into workable shape
get_era5 <- function(csv_path) {
  out_file <- file.path('data', 'era5', 'clean', basename(csv_path))
  if (!file.exists(out_file)) {
    csv_path %>%
      read_csv(col_names = column_headers, col_types = cols()) %>%
      mutate(date = as.Date(paste0(year, '-01-01')) + day_of_year - 1, 
             datetime_utc = ymd_hms(paste0(as.character(date), " ", hour, ":00:00"))) %>%
      select(-date) %>%
      write_csv(path = out_file)
  }
  file.exists(out_file)
}

pblapply(era5_csv_paths, get_era5)
