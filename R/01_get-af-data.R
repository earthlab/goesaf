library(sf)
library(mapview)
library(tidyverse)
library(RCurl)
library(XML)
library(parallel)
library(pbapply)
library(assertthat)
library(data.table)

baseurls <- c("https://fsapps.nwcg.gov/afm/data_goes/kml/conus_hist/")

get_anim_urls <- function(url) {
  result <- getURL(url, verbose = TRUE, ftp.use.epsv = TRUE, dirlistonly = TRUE)
  anim_files <- getHTMLLinks(result, 
                             xpQuery = "//a/@href[contains(., 'animation')]")
  anim_urls <- paste0(url, anim_files)
  anim_urls
}

anim_urls <- lapply(baseurls, get_anim_urls) %>%
  unlist

anim_src <- anim_urls %>%
  str_split('/') %>%
  lapply(function(x) x[[5]]) %>% # fifth element defines source (MODIS vs. GOES)
  unlist %>%
  ifelse(. == 'data', 'modis', .) %>%
  gsub('data_', '', .)

anim_files <- paste(anim_src, basename(anim_urls), sep = '_')

dest_files <- file.path('data', anim_files)
dir.create('data', showWarnings = FALSE)

for (i in seq_along(anim_urls)) {
  if (!file.exists(dest_files[i])) {
    download.file(anim_urls[i], dest_files[i])
  }
}

us_natl_atlas_epsg <- 2163

dir.create('data/fired')
system("aws s3 cp s3://earthlab-natem/modis-burned-area/delineated_events/modis_event_polygons_cus.gpkg data/fired/events_w_attributes.gpkg")
fired <- st_read("data/fired/events_w_attributes.gpkg")


fire_events <- fired %>%
  filter(ignition_date > as.Date('2014-01-07')) %>%
  mutate(duration = last_date - ignition_date) %>%
  st_transform(crs = us_natl_atlas_epsg)

# for each fire event, extract relevant active fire detections
load_af <- function(x) {
  outfile <- file.path('out', paste0('af_', x$id, '.shp'))
  
  if (!file.exists(outfile)) {
    dates <- seq(x$ignition_date, x$last_date + 7, by = 1)
    af_files <- lapply(format(dates, '%Y%m%d'), function(x) {
      list.files(path = 'data', pattern = x, full.names = TRUE)
    })
    af_files <- unlist(af_files)
    corrupt_files <- c(
      'data/goes_conus_time_animation_20140520.kmz', 
      'data/goes_conus_time_animation_20140618.kmz', 
      'data/goes_conus_time_animation_20160702.kmz',
      "data/goes_conus_time_animation_20140310.kmz"
    )

    af_files <- af_files[!(af_files %in% corrupt_files)]
    af_files <- af_files[af_files %in% list.files(path = 'data', 
                                                        full.names = TRUE)]
    res <- list()
    for (f in seq_along(af_files)) {
      
      layers <- grep(pattern = 'Fire Detection Footprints \\(Last', 
           x = st_layers(af_files[f])$name, value = TRUE)
      
      # read the different layers and save any detections in an sf object
      tmp_list <- vector(mode = 'list', length = length(layers))
      names(tmp_list) <- layers
      for (layer_idx in seq_along(layers)) {
        tmp_list[[layer_idx]] <- st_read(af_files[f], layer = layers[layer_idx])
        
        if (nrow(tmp_list[[layer_idx]]) > 0) {
          tmp <- tmp_list[[layer_idx]] %>%
            st_zm %>%
            st_transform(st_crs(x)$epsg) %>%
            filter(st_intersects(., x, sparse=FALSE))
          tmp_list[[layer_idx]] <- tmp
        }
      }
      
      nrows <- lapply(tmp_list, nrow) %>%
        unlist
      if (any(nrows > 0)) {
        res[[f]] <- sf::st_as_sf(data.table::rbindlist(tmp_list[nrows > 0], 
                                                       fill = TRUE))
      } else {
        res[[f]] <- tmp_list[[1]]
      }
    }
    af_rows <- lapply(res, nrow) %>%
      unlist
    if (any(af_rows > 0)) {
      sf::st_as_sf(data.table::rbindlist(res[af_rows > 0])) %>%
        sf::st_write(outfile)
    } else {
      # create empty file
      file.create(outfile)
    }
  }
}

pb <- txtProgressBar(max = nrow(fire_events), style = 3)
split_events <- vector(mode = 'list', length = nrow(fire_events))
for (i in 1:nrow(fire_events)) {
  split_events[[i]] <- fire_events[i, ]
  setTxtProgressBar(pb, i)
}
close(pb)

dir.create("out", showWarnings = FALSE)
pboptions(use_lb = TRUE)
cl <- makeCluster(parallel::detectCores())
clusterEvalQ(cl = cl, library(tidyverse))
clusterEvalQ(cl = cl, library(sf))
pblapply(split_events, load_af, cl = cl)
stopCluster(cl)
