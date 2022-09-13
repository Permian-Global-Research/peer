#' get projection of spatial object
#'
#' cloud mask an l8 image or collection
#' @title l8_mask_clouds: cloud mask an l8 image or collection
#' @param x An ee Image or ImageCollection
#'  @param ... Not currently used.
#' @export
l8_mask_clouds <- function(x, ...) {
  UseMethod("l8_mask_clouds")
}

#' @rdname l8_mask_clouds
#'
#' @export
l8_mask_clouds.ee.image.Image <- function(x, ...){
  mask_clouds_l8_method(x)
}

#' @rdname l8_mask_clouds
#'
#' @export
l8_mask_clouds.ee.imagecollection.ImageCollection <- function(x, ...){
  x$map(mask_clouds_l8_method)
}
