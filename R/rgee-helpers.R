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



