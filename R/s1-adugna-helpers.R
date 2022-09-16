#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# """
# Version: v1$2
# Date: 2021-02-11
# Authors: Mullissa A$, Vollrath A$, Braun, C$, Slagter B$, Balling J$, Gou Y$, Gorelick N$,  Reiche J$
# """

# import ee

# ---------------------------------------------------------------------------//
# Linear to db scale
# ---------------------------------------------------------------------------//

lin_to_db <- function(image){
  # """
  #   Convert backscatter from linear to dB$
  #
  #   Parameters
  #   ----------
  #   image : ee$Image
  #       Image to convert
  #
  #   Returns
  #   -------
  #   ee$Image
  #       output image
  #
  #   """
  bandNames = image$bandNames()$remove('angle')
  db = ee$Image$constant(10)$
    multiply(image$select(bandNames)$
               log10())$rename(bandNames)
  image$addBands(db, NULL, TRUE)
}



db_to_lin <- function(image){
  # """
  #   Convert backscatter from dB to linear$
  #
  #   Parameters
  #   ----------
  #   image : ee$Image
  #       Image to convert
  #
  #   Returns
  #   -------
  #   ee$Image
  #       output image
  #
  #   """
  bandNames = image$bandNames()$remove('angle')
  lin = ee$Image$constant(10)$
    pow(image$select(bandNames)$
          divide(10))$rename(bandNames)
  image$addBands(lin, NULL, TRUE)
}


lin_to_db2 <- function(image){
  # """
  #   Convert backscatter from linear to dB by removing the ratio band$
  #
  #   Parameters
  #   ----------
  #   image : ee$Image
  #       Image to convert
  #
  #   Returns
  #   -------
  #   ee$Image
  #       Converted image
  #
  #   """
  db = ee$Image$constant(10)$
    multiply(image$select(c('VV', 'VH'))$
               log10())$rename(c('VV', 'VH'))
  image$addBands(db, NULL, TRUE)
}


# ---------------------------------------------------------------------------//
# Add ratio bands
# ---------------------------------------------------------------------------//

add_ratio_lin <- function(image){
  # """
  #   Adding ratio band for visualization
  #
  #   Parameters
  #   ----------
  #   image : ee$Image
  #       Image to use for creating band ratio
  #
  #   Returns
  #   -------
  #   ee$Image
  #       Image containing the ratio band
  #
  #   """
  ratio = image$addBands(image$select('VV')$
                           divide(image$select('VH'))$
                           rename('VVVH_ratio'))

  ratio$set('system:time_start', image$get('system:time_start'))

}
