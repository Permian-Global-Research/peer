#' dem_as_image
#'
#' create a DEM image - choose from different sources. Currently only GLO30
#' and SRTM.
#'
#' @param aoi
#' @param dem.src
#'
#' @return
#' @export
#'
#' @examples
dem_as_image <- function(aoi, dem.src=c("glo30", "srtm"), mask=FALSE){ #
  aoi_ee <- return_ee_geom(aoi)

  if (isFALSE(mask)){
    aoi_ee <- ee_geom_bounds(aoi_ee)
  }

  if (dem.src[1]=="glo30"){
    dem.collect<- "projects/sat-io/open-datasets/GLO-30"
    dem <- rgee::ee$ImageCollection(dem.collect)$
      # filterBounds(collection$geometry())$
      mosaic()$
      setDefaultProjection('EPSG:3857',NULL,30)$
      clip(aoi_ee)

  } else if (dem.src[1]=="srtm") {
    dem <- ee$Image("USGS/SRTMGL1_003")$
      clip(aoi_ee)

  } else {
    stop('Only "glo30" or "srtm" Digital Elevation Models are supprted.')
  }
  return(dem)
}
