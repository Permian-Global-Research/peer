#' Terrain flatenning for s1 image
#'
#' Applies terrain flattening alg by Andreas Vollrath (ESA) to an ee image.
#'
#' @param img an ee image
#'
#' @return a terrain flattened image with additional "layover" mask band.
#' @details Implementation by Andreas Vollrath (ESA), inspired by Johannes Reiche (Wageningen)
#' @export
#'
#' @examples
terrain_flatten <- function(img) {

  img <- img$select(c("VV", "VH", "angle"))

  imgGeom <- img$geometry()
  srtm <- ee$Image('USGS/SRTMGL1_003')$clip(imgGeom) # 30m srtm
  sigma0Pow <- ee$Image$constant(10)$pow(img$divide(10.0))

  # Article ( numbers relate to chapters)
  # 2$1$1 Radar geometry
  theta_i <- img$select('angle')
  phi_i <- ee$Terrain$aspect(theta_i)$
    reduceRegion(ee$Reducer$mean(), theta_i$get('system:footprint'), 1000)$
    get('aspect')

  # 2$1$2 Terrain geometry
  alpha_s <- ee$Terrain$slope(srtm)$select('slope')
  phi_s <- ee$Terrain$aspect(srtm)$select('aspect')

  # 2$1$3 Model geometry
  # reduce to 3 angle
  phi_r <- ee$Image$constant(phi_i)$subtract(phi_s)

  # convert all to radians
  phi_rRad <- phi_r$multiply(pi / 180)
  alpha_sRad <- alpha_s$multiply(pi / 180)
  theta_iRad <- theta_i$multiply(pi / 180)
  ninetyRad <- ee$Image$constant(90)$multiply(pi / 180)

  # slope steepness in range (eq$ 2)
  alpha_r <- (alpha_sRad$tan()$multiply(phi_rRad$cos()))$atan()

  # slope steepness in azimuth (eq 3)
  alpha_az <- (alpha_sRad$tan()$multiply(phi_rRad$sin()))$atan()

  # local incidence angle (eq$ 4)
  theta_lia <- (alpha_az$cos()$multiply((theta_iRad$subtract(alpha_r))$cos()))$acos()
  theta_liaDeg <- theta_lia$multiply(180 / pi)
  # 2$2
  # Gamma_nought_flat
  gamma0 <- sigma0Pow$divide(theta_iRad$cos())
  gamma0dB <- ee$Image$constant(10)$multiply(gamma0$log10())

  # ratio_1 = gamma0dB$select('VV')$subtract(gamma0dB$select('VH'))

  # Volumetric Model
  nominator <- (ninetyRad$subtract(theta_iRad)$add(alpha_r))$tan()
  denominator <- (ninetyRad$subtract(theta_iRad))$tan()
  volModel <- (nominator$divide(denominator))$abs()

  # apply model
  gamma0_Volume <- gamma0$divide(volModel)
  gamma0_VolumeDB <- ee$Image$constant(10)$multiply(gamma0_Volume$log10())

  # we add a layover/shadow maskto the original implmentation
  # layover, where slope > radar viewing angle
  alpha_rDeg <- alpha_r$multiply(180 / pi)
  layover <- alpha_rDeg$lt(theta_i)$rename("layover")

  # shadow where LIA > 90
  shadow <- theta_liaDeg$lt(85)$
    rename("shadow")

  # calculate the ratio for RGB vis
  # ratio = gamma0_VolumeDB$select('VV')$subtract(gamma0_VolumeDB$select('VH'))

  output <- gamma0_VolumeDB$
    addBands(layover)

  output
}


#' Process Sentinel 1 image or collection
#'
#' @param x an ee image or collection
#' @param aoi.sf an sf/sfc object of the area of interest
#' @param mask
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
s1_process <- function(x, aoi.sf, mask=FALSE, ...) {
  UseMethod("s1_process")
}

#' @rdname s1_process
#'
#' @export
s1_process.ee.image.Image <- function(x, aoi.sf){

  aoi <- sf_ext_as_ee(aoi.sf)

  x <- s1_wrapper(x)
  x <- subset_bounds(x, aoi)

  x
}

#' @rdname s1_process
#'
#' @export
s1_process.ee.imagecollection.ImageCollection <- function(x, aoi.sf){
  aoi <- sf_ext_as_ee(aoi.sf)

  x <- x$map(combine_funcs)
  x <- x$map(subset_bounds, aoi)
  x
}


#' s1 processing wrapper
#'
#' @param x ee image
#'
#' @return ee image
s1_wrapper <- function(x){
    img <- terrain_flatten(x)|>
      toNatural()|>
      focalSpeckle()|>
      todB() |>
      s1_add_ratio()

  }
