#' img_band_names
#'
#' get ee_image band names as a raster
#'
#' @param img
#'
#' @return character vector of band names.
#' @export
#'
#' @examples
img_band_names <- function(img){
  img$bandNames()$getInfo()
}


#' helper to return ee_featurecollection regardless of class
#'
#' @param x either sf or ee.featurecollection.FeatureCollection
#'
#' @return  ee.featurecollection.FeatureCollection
return_ee_geom <- function(x){
  if (inherits(x, c("sf"))){ #, "sfc_POLYGON" sfc polygon behaves weird after being converted to ee...
    return(rgee::sf_as_ee(x))
  } else if(inherits(x, "ee.featurecollection.FeatureCollection")){
    return(x)
  }
  else if(inherits(x, "ee.geometry.Geometry")){
    return(x)
  }
  else {
    stop(paste0(class(x)[1], " is not a supported type use sf of ee feature collection"))
  }
}



ee_geom_bounds <- function(x){
  x$union()$geometry()$bounds()
}


#' sf_ext_as_ee
#'
#' convert an sf object into a ee Feature Collection describing the bounding box
#'  of the sf object.
#'
#' @param x sf of sfc object.
#'
#' @return
#' @export
#'
#' @examples
sf_ext_as_ee <- function(x){
  .crs <- sf::st_crs(x)
  x|>
    sf::st_bbox()|> sf::st_as_sfc(crs=.crs) |>
    sf::st_as_sf()|>
    rgee::sf_as_ee()
}



#' clip an image add properties.
#'
#' @param img
#' @param aoi
#'
#' @return clipped img
subset_bounds <- function(img, aoi) {
  # // Crop by table extension
  img$clip(aoi)$
    copyProperties(img,c('system:time_start','system:time_end'))
}
