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
  if (inherits(x, "sf")){
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
