
<!-- README.md is generated from README.Rmd. Please edit that file -->

# peer

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

The goal of peer is to provide easy access to useful
[{rgee}](https://r-spatial.github.io/rgee/) functions - like composites
etc.

## Installation

You can install the development version of peerm from
[GitHub](https://github.com/):

first clone the library - SSH probably best option. Then open the
project and run the following in R:

``` r
install.packages("remotes")
remotes::install_github(repo="Permian-Global-Research/peer",
                        auth_token = "[MY-TOKEN]")
```

## Example

Quick example creating a False Colour Composite from Landsat 8 data.

``` r
library(rgee)
ee_Initialize(quiet=TRUE)
library(peer)
library(sf)

# get Landsat 8 collection from date range and cloud cover filters.
comp2021 <- l8_collect(aoi=kuamut,
                       start.date = '2021-01-01',
                       end.date= '2021-12-31',
                       min.cloud=50)|>
  l8_mask_clouds()|>   #function mask clouds
  l8_terrain_correct(aoi.sf = kuamut) # function to run terrain correction

x_img <- comp2021$
  map(l8_add_evi)$ # function to add EVI 
  map(l8_add_ndvi)$ # function to add NDVI
  median()

# x_img <- l8_mask_clouds(x)$ # landsat 8 cloud and shadow mask.
#   median()


vizParams <- list(
  bands = c('B3', 'B4', 'B5'),
  min= 0,
  max= 1000,
  gamma= c(1.1, 0.95, 0.7)
)

aoi_cent <- st_centroid(kuamut) %>%
  st_coordinates()

Map$setCenter(lon = aoi_cent[1], lat = aoi_cent[2], zoom = 10)
m<-Map$addLayer(x_img, vizParams, 'FC-composite-2021')
print(m)
```

![FCC-example](man/FCC-21.png)

``` r

s1.test <- s1_collect(kuamut, start.date = '2016-01-01',
                      end.date= '2016-12-31', orbit.pass = 'DESCENDING' )#$

x.rtc <- s1_process(s1.test)

x.no.rtc <- s1_process(s1.test, rtc = FALSE)

MapAll(x.no.rtc) +
  MapAll(x.rtc)


#--- OR USE METHODS FROM ADUGNA ET AL. 2021. ONLY DEFAULTS WORKING FOR NOW WITH:

s1.proc <- s1_adugna_process(s1.test, kuamut)
```
