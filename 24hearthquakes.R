library(devtools)
library(httr)
library(dplyr)
devtools::install_github("hadley/purrr")
library(purrr)

previous.day.earthquakes <- function() {
  
  # Use HTTR to query the last 24 hours' earthquakes from USGS.
  
  now <- as.POSIXct(Sys.time(), "UTC")
  since <- now - 60 * 60 * 24
  since.str <- format(since, "%y-%m-%dT%H:%M:%S")
  resp <- GET("http://earthquake.usgs.gov/fdsnws/event/1/query", query=list(
    format="geojson", 
    starttime=since.str
  ))
  features <- httr::content(resp, "parsed")$features
  
  # Extract latitude, longitude and magnitude for each earthquake.
  
  mag <- features %>% purrr::map(~ .$properties$mag) %>% unlist
  long <- features %>% purrr::map(~ .$geom$coordinates[[1]]) %>% unlist
  lat <- features %>% purrr::map(~ .$geom$coordinates[[2]]) %>% unlist
  
  # Use the htmlwidgets globejs function to create and return an
  # interactive globe.
  
  earth <- "/Users/ggospodinov/Downloads/Land_shallow_topo_2048.jpg"
  globejs(img=earth, bodycolor="#555555", emissive="#444444",
          lightcolor="#555555", bg="#ffffff", lat=lat, long=long,
          color="#FF3333",
          value=mag * 50)
  
}