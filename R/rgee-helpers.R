#' Get nth item from an ee collection
#'
#' @param collect
#' @param n
#'
#' @return The entire unclipped image at the nth position of the collection.
#' @export
#'
#' @examples
#' s1.test <- s1_collect(kuamut, '2019-01-01', '2019-12-31' )
#' s1_nth(s1.test, 3)
#'
get_nth <- function(collect, n){
  # p <- collect$first()$propertyNames()
  i <- collect$aggregate_array('system:id')$get(n)
  ee$Image(i$getInfo())
}



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
  .crs <- sf::st_crs(x)$wkt
  x <- x|>
    sf::st_bbox()|>
    sf::st_as_sfc(crs=.crs) |>
    sf::st_as_sf()|>
    rgee::sf_as_ee()

  x$
    geometry()$
    bounds()
}



#' clip an image add properties.
#'
#' @param img
#' @param aoi
#'
#' @return clipped img
subset_bounds <- function(img, aoi) {
  # // Crop by table extension
  x <- img$clip(aoi)
  x$copyProperties(img,c('system:time_start','system:time_end'))
  x
}

#' clip an image add properties.
#'
#' @param img
#' @param aoi
#'
#' @return clipped img
sub_bounds_func <- function(aoi) {
  # // Crop by table extension
  function(x){
    x <- x$clip(aoi)
    x$copyProperties(x,c('system:time_start','system:time_end'))
    x
  }

}

#' buffer_points
#'
#' @param pnts
#' @param radius
#' @param rectangle
#'
#' @return
#' @export
#'
#' @examples
buffer_points <- function (pnts, radius, rectangle=FALSE) {

  buff_point <- function(pt){
    if (rectangle){
      pt$buffer(radius)$bounds()
    } else {
      pt$buffer(radius)
    }
  }

  pnts$map(buff_point)
}


#' point_extract
#'
#' @param img ee image
#' @param ee_fc an ee feature collection (points)
#'
#' @return a dataframe of values for all properties
#' @export
#'
#' @examples
point_extract <- function(img, ee_fc) {
  i <- img$sampleRegions(collection = ee_fc, # feature collection here
                         scale = 10)$ # Cell size of raster
    getInfo()$
    features

  .names <- names(i[[1]]$properties)

  vals <- .names|>
    lapply(function(x) {
      sapply(i, function(y){
        y$properties[x]
      })
    }) |>
    lapply(as.numeric)

  names(vals)<- .names
  return(as.data.frame(vals))
}



#' ee_ahist
#'
#' gives approximate histograms for all bands in an ee image based on random points
#'
#' @param img an earth engine image to sample from
#' @param ee_fc an earth engine feature class of points to use. Useful when
#' comparing between images to maintain consistent sampling. Whn NULL, the points
#' are randomly generated.
#' @param n The number of random points to sample if ee_fc is NULL. Lower this value
#' when the this error occurs: "User memory limit exceeded."
#'
#' @return a ggplot object with the histogram(s)
#' @export
#'
#' @examples
ee_ahist <- function(img, ee_fc=NULL, n=1000){

  if (is.null(ee_fc)){
    ee_fc <- rando_points(img, n)
  }

  s <- point_extract(img, ee_fc)
  # browser()
  nb.cols <- length(unique(names(s)))
  mycolors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(nb.cols)


  s |>
    tidyr::pivot_longer(tidyselect::everything())|>
    ggplot2::ggplot()+
    ggplot2::aes(x=value, fill=name)+
    ggplot2::geom_histogram(ggplot2::aes(y=..density..), bins=30)+
    ggplot2::geom_density(fill=NA) +
    ggplot2::facet_wrap(~name, scale='free') +
    ggplot2::scale_fill_manual(values = mycolors, name='Band Name') +
    ggplot2::guides(fill="none") +
    ggplot2::theme_light()
}

#' rando_points
#'
#' creaets random points across the extent of an ee image.
#'
#' @param img and ee image
#' @param n default is 1000. the number of points to generate. Note Earth Engine
#' has some limits on how many values can be returned.
#'
#' @return
#' @export
#'
#' @examples
rando_points <- function(img, n=1000){
  ee$FeatureCollection$
    randomPoints(img$geometry(), n)
}


anti_selection <- function(img, bands){
  all_bands <- img$bandNames()
  img$select(all_bands$removeAll(bands))
}




# ee_reproj <- function(img){
#   prj <- ee$Projection('EPSG:4326')
#   img$reproject(proj, NULL, 30)
# }
# check_col_empty <- function(){

# }
