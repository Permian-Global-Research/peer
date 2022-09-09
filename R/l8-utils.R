

#' getQABits
#'
#' @param img
#' @param start
#' @param end
#' @param newName
#'
#' @return
getQABits <- function(img, start, end, newName) {
  # // Compute the bits we need to extract.
  pattern <- sapply(c(start:end), function(x) 2 ^x) |>
    sum()
  # // Return a single band image of the extracted QA bits, giving the band
  # // a new name.
  img$select(0)$
  rename(newName)$bandNames()
  # img$select(0, newName)$
  img$
   bitwiseAnd(pattern)$
  rightShift(start)
}

# // A function to mask out cloudy pixels.
remove_cloud_n_shadows = function(img, base_start, base_end) {
  # // Select the QA band.
  QA <- img$select('pixel_qa')
  # // Get the internal_cloud_algorithm_flag bit.
  getQABits(QA, base_start,base_end, 'cloud_shadows')$eq(0)
  # // Return an image masking out cloudy areas.
}



#' mask_clouds_l8_method
#'
#' @param img an earth engine ee image
#'
#' @return a cloud masked ee image
mask_clouds_l8_method <- function(img){
    cl.sh <- remove_cloud_n_shadows(img, 3, 3)
    cl <- remove_cloud_n_shadows(img, 5, 5)
    img = img$updateMask(cl.sh)
    img = img$updateMask(cl)
    # img <- img$addBands(c(cl.sh, cl))
    img
}

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



#' l8_collect
#'
#' @param aoi
#' @param start.date
#' @param end.date
#' @param min.cloud
#'
#' @return
#' @export
#'
#' @examples
l8_collect <- function(aoi, start.date, end.date, min.cloud=90){

  aoi_ee <- sf_ext_as_ee(aoi)

  subset_bounds <- function(img) {
    # // Crop by table extension
    img$clip(aoi_ee)$
      copyProperties(img,c('system:time_start','system:time_end'))
  }

  ee$ImageCollection('LANDSAT/LC08/C01/T1_SR')$
    filterBounds(aoi_ee)$
    filter(ee$Filter$lt('CLOUD_COVER', min.cloud))$
    map(subset_bounds)$
    filterDate(start.date, end.date)
}


#' l8_add_ndvi
#'
#' @param img
#'
#' @return
#' @export
#'
#' @examples
l8_add_ndvi <- function(img) {
  ndvi <- img$normalizedDifference(c('B5', 'B4'))$rename('NDVI');
  img$addBands(ndvi)
}


#' l8_add_evi
#'
#' @param img
#'
#' @return
#' @export
#'
#' @examples
l8_add_evi <- function(img){
  evi <- img$expression(
    expression = '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
    opt_map =  list(
      'NIR' = img$select('B5'),
      'RED' = img$select('B4'),
      'BLUE' = img$select('B2')
    )
  )$rename('EVI')

  img$addBands(evi)

}


