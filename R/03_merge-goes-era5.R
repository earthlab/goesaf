library(tidyverse)
library(sf)
library(plotly)
library(viridis)
library(lubridate)
library(lutz)
library(pbapply)
library(suncalc)
library(parallel)
library(gganimate)
library(lme4)
library(effects)

af_files <- list.files(path = 'out', pattern = 'af_.*\\.shp$', 
                       full.names = TRUE)
valid_af <- af_files[file.size(af_files) > 0]
af_sizes <- file.size(valid_af)

us_natl_atlas_epsg <- 2163

fired <- st_read("data/fired/events_w_attributes.gpkg")

fire_events <- fired %>%
  filter(ignition_date > as.Date('2014-01-07')) %>%
  mutate(duration = last_date - ignition_date) %>%
  st_transform(crs = us_natl_atlas_epsg)
  
event_centroids <- fire_events %>%
  st_transform(crs = 4326) %>%
  st_centroid() %>%
  st_coordinates %>%
  as_tibble %>%
  rename(lon = X, lat = Y)

fire_events <- fire_events %>%
  bind_cols(event_centroids)

# write a csv file of fire events
fire_events %>%
  as.data.frame %>%
  select(-geom) %>%
  as_tibble %>%
  mutate(last_date_plus_7_days = last_date + 7) %>%
  select(id, ignition_date, last_date_plus_7_days, lon, lat) %>%
  write_csv('out/fire-events-for-era5.csv')

parse_detections <- function(file) {
  detections <- st_read(file, quiet = TRUE) %>%
    mutate(detection_time = substring(dscrptn, regexpr("[0-9]{2}:[0-9]{2} UTC",
                                                       dscrptn)), 
           detection_time = gsub(pattern = 'UTC.*', 'UTC', detection_time), 
           detection_datetime_utc = ymd_hm(paste(timstmp, detection_time)), 
           id = gsub(pattern = ".*_", replacement = "", file), 
           id = gsub("\\.shp", "", id), 
           id = parse_number(id)) %>%
    left_join(as.data.frame(fire_events) %>%
                as_tibble %>%
                dplyr::select(-geom)) %>%
    mutate(local_tz = tz_lookup_coords(lat, lon, method = 'accurate'), 
           detection_local = local_time(detection_datetime_utc, local_tz, 
                                        units = 'hours'), 
           detection_date = as.Date(timstmp))
  
  sunlight_times <- getSunlightTimes(date = unique(detections$detection_date), 
                                     lat = unique(detections$lat), 
                                     lon = unique(detections$lon), 
                                     keep = c('sunrise', 'sunset'), 
                                     tz = unique(detections$local_tz)) %>%
    rename(detection_date = date)
  
  sun_data <- getSunlightPosition(date = unique(detections$detection_datetime_utc), 
                                  lat = unique(detections$lat), 
                                  lon = unique(detections$lon)) %>%
    rename(detection_datetime_utc = date)
  
  out <- detections %>%
    left_join(sunlight_times) %>%
    left_join(sun_data) %>%
    rowwise %>%
    mutate(is_day = ifelse(altitude > 0, 'daytime', 'nighttime')) %>%
    ungroup %>%
    st_sf %>%
    st_cast('MULTIPOLYGON') %>%
    mutate(sensor = case_when(
      grepl("GOES", dscrptn) ~ "GOES", 
      grepl("NOAA-", dscrptn) ~ "AVHRR", 
      grepl("MODIS", dscrptn) ~ "MODIS", 
      TRUE ~ "unknown"
      ), 
      datetime_utc = round_date(detection_datetime_utc, unit = 'hour'))
  out
}

test <- parse_detections(valid_af[which.max(af_sizes)])

# sanity check - did we calculate solar altitude correctly?
# if so, it should peak in the middle of the day
test %>%
  ggplot(aes(x = detection_local, y = altitude, color = is_day)) + 
  geom_point()

cl <- makeCluster(parallel::detectCores())
clusterEvalQ(cl = cl, library(tidyverse))
clusterEvalQ(cl = cl, library(sf))
clusterEvalQ(cl = cl, library(lubridate))
clusterEvalQ(cl = cl, library(suncalc))
clusterEvalQ(cl = cl, library(lutz))
clusterExport(cl = cl, varlist = c('fire_events'))
detections <- pblapply(valid_af, parse_detections, cl = cl)
stopCluster(cl)




# Merge ERA5 data to GOES active fire data --------------------------------

joined <- sf::st_as_sf(data.table::rbindlist(detections))


mod_d <- joined %>%
  as.data.frame() %>%
  select(-geometry) %>%
  as_tibble %>%
  filter(sensor == 'GOES') %>%
  group_by(id, datetime_utc) %>%
  summarize(n = n()) %>%
  ungroup %>%
  group_by(id) %>%
  # add one day to the final hour to represent zeros
  tidyr::complete(datetime_utc = seq(min(datetime_utc), max(datetime_utc) + 60*60*24, 
                              by = 'hour'), 
           fill = list(n = 0)) %>%
  ungroup %>%
  arrange(id, datetime_utc)

ids <- unique(mod_d$id)

era5_files <- file.path('data', 'era5', 'clean', 
                        paste0('wx', ids, '.csv'))
# <1% of (40) fire ids have no era5 data - remove those
ids_to_exclude <- ids[!file.exists(era5_files)]
ids_to_keep <- ids[file.exists(era5_files)]
mod_d <- filter(mod_d, id %in% ids_to_keep)

era5_data <- file.path('data', 'era5', 'clean', 
                       paste0('wx', ids_to_keep, '.csv')) %>%
  pblapply(read_csv, col_types = cols()) %>%
  bind_rows(.id = 'file') %>%
  mutate(file = parse_number(file), 
         id = ids_to_keep[file]) %>%
  select(-file)

af_data <- mod_d %>%
  right_join(era5_data) %>%
  mutate(n = ifelse(is.na(n), 0, n))

write_csv(af_data, 'out/goes-af-era5.csv')
