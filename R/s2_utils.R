
#' s2_collect
#'
#' @param aoi
#' @param start.date
#' @param end.date
#' @param min.cloud
#' @param dataset
#'
#' @return
#' @export
#'
#' @examples
s2_collect <- function(aoi, start.date, end.date, min.cloud=90,
                       dataset =c("S2_SR", "S2_SR_HARMONIZED")){

  aoi_ee <- sf_ext_as_ee(aoi)


  collec <- ee$ImageCollection(paste0('COPERNICUS/', dataset[1]))$
    filterBounds(aoi_ee)$
    filter(ee$Filter$lt('CLOUDY_PIXEL_PERCENTAGE', min.cloud))$
    filterDate(start.date, end.date)

  if (length(collec$getInfo()$features)>0){
    return(collec)
  } else {
    stop(paste0("The collection is empty - perhaps try a different ",
                "min.cloud value or change the start.date and ",
                "end.date values"))
  }
}
