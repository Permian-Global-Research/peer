

#' Create Analysis ready Sentinel 1 image or collection
#'
#' Using tools from Adugna, et al 2021 (https://www.mdpi.com/2072-4292/13/10/1954)
#'
#' @param x an ee image or collection
#' @param reduce Only for image colections - should the collection be reduced to
#' an image.
#' @param mask Not used
#' @param ... not used
#' @param aoi.sf sf object to set extent for clipping.
#' @param rtc Logical - should Terrain Correction/Flattening be applied?
#'
#' @return
#' @export
#'
#' @examples
s1_adugna_process <- function(x, aoi.sf, rtc=TRUE, reduce=TRUE, mask=FALSE, ...) {
  UseMethod("s1_adugna_process")
}

# #' @rdname s1_adugna_process
# #'
# #' @export
# s1_adugna_process.ee.image.Image <- function(x, aoi.sf){
#
#   aoi <- sf_ext_as_ee(aoi.sf)
#
#   .nam <- c("VV", "VH", "Ratio")
#   # if (isTRUE(rtc)){
#   x <- s1_adugna_wrapper_rtc(x)
#   .nam <- paste0(.nam, "-RTC")
#   # } else {
#   #   x <- s1_adugna_wrapper(x)
#   # }
#   vv <- x$select('VV')$clip(aoi)$rename(.nam[1])
#   vh <- x$select('VH')$clip(aoi)$rename(.nam[2])
#   rat <- x$select('VVVHRatio')$clip(aoi)$rename(.nam[3])
#   x <- vv$addBands(vh)$addBands(rat)
#
#   x
# }

#' @rdname s1_adugna_process
#'
#' @export
s1_adugna_process.ee.imagecollection.ImageCollection <- function(x, aoi.sf, reduce=TRUE){

  aoi <- sf_ext_as_ee(aoi.sf)

  .nam <- c("VV", "VH", "Ratio")

  # if (isTRUE(rtc)){
    x <- s1_adugna_wrapper(x)
    .nam <- paste0(.nam, "-RTC")
  # } else {
    # x <- x$map(s1_wrapper)
  # }

  # x <- x$map(subset_bounds, aoi)
  if (reduce){
    vv.mean <- x$select('VV')$
      mean()$
      clip(aoi)$
      rename(.nam[1])
    vh.mean <- x$select('VH')$
      mean()$
      clip(aoi)$
      rename(.nam[2])
    rat.mean <- x$select('VVVH_ratio')$
      mean()$
      clip(aoi)$
      rename(.nam[3])

    x <- vv.mean$addBands(vh.mean)$addBands(rat.mean)
  }
  return(x)
}


# #' s1 RTC processing wrapper
# #'
# #' @param x ee image
# #'
# #' @return ee image
# # s1_adugna_wrapper_rtc <- function(x){
#
#   # x <- terrain_flatten(x)|>
#   s1_adugna_wrapper(x)|
# }

#' s1 processing wrapper
#'
#' @param x ee image
#'
#' @return ee image
s1_adugna_wrapper <- function(x){

  x <- x$map(f_mask_edges)

  x <- ee$ImageCollection(MonoTemporal_Filter(x,
                                         KERNEL_SIZE = 9,
                                         SPECKLE_FILTER = 'GAMMA MAP'))
    slope_correction(x,
      TERRAIN_FLATTENING_MODEL = 'VOLUME',
      DEM = ee$Image('USGS/SRTMGL1_003'),
      TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER = 0)$
    map(add_ratio_lin)$
    map(lin_to_db2)


}
