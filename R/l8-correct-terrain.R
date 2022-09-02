

#' l8_terrain_correct
#'
#' correct illumination effect for Landsat8 collection.
#'
#' @param x either an earth engine Collection or Image
#' @param aoi
#' @param dem.src
#'
#' @return
#' @export
#'
#' @examples
l8_terrain_correct <- function(x, aoi,  dem.src=c("glo30", "srtm"), buf=10000) {

  dem <- dem_as_image(aoi, dem.src)

  combine_illumes <- function(img){
    img <- illuminationCondition(img, dem, buf)
    illuminationCorrection(img, dem, buf)
  }


  if (inherits(x, "ee.image.Image" )){
    return(combine_illumes(x))
  } else if (inherits(x, "ee.imagecollection.ImageCollection")){
    return(x$map(illuminationCondition))
  }

  # collection = collection$map(illuminationCondition, dem)
  #
  # collection = collection$map(illuminationCorrection, dem)
  #   return(collection)
}


illuminationCondition <- function (img, dem, buf=10000){

  # // Extract image metadata about solar position
  SZ_rad = ee$Image$constant(ee$Number(img$get('SOLAR_ZENITH_ANGLE')))$multiply(pi)$divide(180)$clip(img$geometry()$buffer(buf))
  SA_rad = ee$Image$constant(ee$Number(img$get('SOLAR_AZIMUTH_ANGLE'))$multiply(pi)$divide(180))$clip(img$geometry()$buffer(buf)) # CHECK PARENS HERE! DIFFERENT TO ABOVE.
  # // Creat terrain layers
  slp = ee$Terrain$slope(dem)$clip(img$geometry()$buffer(buf))
  slp_rad = ee$Terrain$slope(dem)$multiply(pi)$divide(180)$clip(img$geometry()$buffer(buf))
  asp_rad = ee$Terrain$aspect(dem)$multiply(pi)$divide(180)$clip(img$geometry()$buffer(buf))

  # // Calculate the Illumination Condition (IC)
  # // slope part of the illumination condition
  cosZ = SZ_rad$cos()
  cosS = slp_rad$cos()
  slope_illumination = cosS$expression("cosZ * cosS",
                                       list('cosZ'= cosZ,
                                            'cosS'= cosS$select('slope')))
  # // aspect part of the illumination condition
  sinZ = SZ_rad$sin()
  sinS = slp_rad$sin()
  cosAziDiff = (SA_rad$subtract(asp_rad))$cos()
  aspect_illumination = sinZ$expression("sinZ * sinS * cosAziDiff",
                                        list(
                                          'sinZ' = sinZ,
                                          'sinS' = sinS,
                                          'cosAziDiff' = cosAziDiff
                                        ))
  # // full illumination condition (IC)
  ic = slope_illumination$add(aspect_illumination)

  # // Add IC to original image
  img_plus_ic = ee$Image(
    img$
      addBands(ic$rename('IC'))$
      addBands(cosZ$rename('cosZ'))$
      addBands(cosS$rename('cosS'))$
      addBands(slp$rename('slope'))
  )
  return (img_plus_ic)
}



# ////////////////////////////////////////////////////////////////////////////////
#   // Function to apply the Sun-Canopy-Sensor + C (SCSc) correction method to each
# // image$ Function by Patrick Burns and Matt Macander
illuminationCorrection <- function (img, dem, buf=10000){
  props = img$toDictionary()
  st = img$get('system:time_start')

  img_plus_ic = img
  mask1 = img_plus_ic$select('nir')$gt(-0.1)

  mask2 <- function(ic){
    ma <- ic$select('slope')$gte(5)
    mb <-  ic$select('IC')$gte(0)
    mc <-  ic$select('B5')$gt(-0.1)
    ma$add(mb)$add(mc)$select('slope')$gte(1)
  }

  img_plus_ic_mask2 = ee$Image(img_plus_ic$updateMask(mask2(img_plus_ic)))

  # // Specify Bands to topographically correct
  bandList = c('B2','B3','B4','B5','B6','B7')
  compositeBands = img$bandNames()
  nonCorrectBands = img$select(compositeBands$removeAll(bandList))

  geom = ee$Geometry(img$get('system:footprint'))$bounds()$buffer(buf)

  apply_SCSccorr <- function (band){
    method = 'SCSc'
    out = img_plus_ic_mask2$select('IC', band)$reduceRegion(
      reducer= ee$Reducer$linearFit(), #// Compute coefficients: a(slope), b(offset), c(b/a)
      geometry= ee$Geometry(img$geometry()$buffer(-5000)), #// trim off the outer edges of the image for linear relationship
      scale= 300,
      maxPixels= buf
    )

    if (is.null(out) | out == "undefined") {
      #out === null || out === undefined
      return (img_plus_ic_mask2$select(band))
    } else {
      out_a = ee$Number(out$get('scale'))
      out_b = ee$Number(out$get('offset'))
      out_c = out_b$divide(out_a)
      # // Apply the SCSc correction
      SCSc_output = img_plus_ic_mask2$expression(
        "((image * (cosB * cosZ + cvalue)) / (ic + cvalue))",
        list(
          'image' = img_plus_ic_mask2$select(band),
          'ic' = img_plus_ic_mask2$select('IC'),
          'cosB' = img_plus_ic_mask2$select('cosS'),
          'cosZ' = img_plus_ic_mask2$select('cosZ'),
          'cvalue' = out_c
        )
      )

      return (SCSc_output)
    }

  }

  SCSccorr<- bandList|>
    lapply(apply_SCSccorr)

  img_SCSccorr = ee$Image(SCSccorr)$
    addBands(img_plus_ic$select('IC'))

  bandList_IC = ee$List(c(bandList, 'IC'))$
    flatten()

  img_SCSccorr = img_SCSccorr$unmask(img_plus_ic$select(bandList_IC))$select(bandList)

  return (
    img_SCSccorr$
      addBands(nonCorrectBands)$
      setMulti(props)$
      set('system:time_start', st)
  )
}
