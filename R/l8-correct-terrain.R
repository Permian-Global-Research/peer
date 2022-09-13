#' Terrain iluminaton correction for landsat 8.
#'
#' correct illumination effect for Landsat8 collection/images.
#'
#' @title l8_mask_clouds: cloud mask an l8 image or collection
#'
#' @param x An landsat8 ee image or image colletion.
#' @param aoi.sf An area of interest as an sf or sfc object.
#' @param buf default 10000. The buffer distance.
#' @param mask currently ignored.
#' @param ... Not currently used
#'
#' @export
l8_terrain_correct <- function(x, aoi.sf, buf=10000,
                               mask=FALSE, ...) {
  UseMethod("l8_terrain_correct")
}

#' @rdname l8_terrain_correct
#'
#' @export
l8_terrain_correct.ee.image.Image <- function(x, aoi.sf, buf=10000,
                                              mask=FALSE){

  dtm.aoi <- l8_terrain_dtm(x, aoi.sf)

  x <- illuminationCondition(x, dtm.aoi, buf)
  return(illuminationCorrection(x, dem, buf)# $ clip(dtm.aoi$aoi)
         )
}

#' @rdname l8_terrain_correct
#'
#' @export
l8_terrain_correct.ee.imagecollection.ImageCollection <- function(x, aoi.sf, buf=10000,
                                              mask=FALSE){

  dtm.aoi <- l8_terrain_dtm(x, aoi.sf)

  combine_illumes <- function(img){
    img <- illuminationCondition(img, dtm.aoi, buf)
    illuminationCorrection(img, dem, buf)
  }

  return(x$map(combine_illumes)) #$map(function(x) x$clip(aoi))
}




l8_RTC_msg <- function(){
  message(crayon::bgCyan(crayon::bgMagenta(
    "This Terrain Correction function may work less well where topographic changes
are extreme. Also, If you are getting strange looking results make sure you
have already run `peer::l8_mask_clouds"
  )))
}


#' l8_terrain_dtm
#'
#' correct illumination effect for Landsat8 collection.
#'
#' @param x either an earth engine Collection or Image
#' @param aoi.sf an sf or sfc object defining the aoi.
#'
#' @return A ee image of SRTM DEM.
l8_terrain_dtm <- function(x, aoi.sf) {
  dem.src ="srtm"

  aoi <- sf_ext_as_ee(aoi.sf)

  dem_as_image(aoi, dem.src)
}



#' illuminationCondition
#'
#' @param img
#' @param dem
#' @param buf
#'
#' @return ...
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



#' Correct Landsat 8
#'
#' Function to apply the Sun-Canopy-Sensor + C (SCSc) correction method to a
#' landsat 8 image. Function by Patrick Burns and Matt Macander.
#'
#' @param img
#' @param dem
#' @param buf
#'
#' @return
#' @export
#'
#' @examples
illuminationCorrection <- function (img, dem, buf=10000){
  props = img$toDictionary()
  st = img$get('system:time_start')

  img_plus_ic = img
  mask1 = img_plus_ic$select('B5')$gt(-0.1)

  mask2 <- function(ic){
    ma <- ic$select('slope')$gte(5)
    mb <-  ic$select('IC')$gte(0)
    mc <-  ic$select('B5')$gt(-0.1)
    ma$add(mb)$add(mc)$select('slope')$gte(1)
  }

  img_plus_ic_mask2 = ee$Image(img_plus_ic$updateMask(mask1)) #mask2(img_plus_ic)

  # // Specify Bands to topographically correct
  bandList = list('B2','B3','B4','B5','B6','B7')
  compositeBands = img$bandNames()
  nonCorrectBands = img$select(compositeBands$removeAll(bandList))

  geom = ee$Geometry(img$get('system:footprint'))$bounds()$buffer(buf)

  apply_SCSccorr <- ee_utils_pyfunc(function (band){
    # method = 'SCSc'

    out = img_plus_ic_mask2$select('IC', band)$reduceRegion(
      reducer= ee$Reducer$linearFit(), #// Compute coefficients: a(slope), b(offset), c(b/a)
      geometry= ee$Geometry(img$geometry()$buffer(-5000)), #// trim off the outer edges of the image for linear relationship
      scale= 300,
      maxPixels= 1000000000
    )

    #check if dictionary has null values
    notNullVals <- out$values()$filter(ee$Filter$notNull(list('item')))
    nVals <- out$size()
    nNotNullVals = notNullVals$size()

    # run the calculation.
    out_a = ee$Number(out$get('scale'))
    out_b = ee$Number(out$get('offset'))
    out_c = out_b$divide(out_a)
    # // Apply the SCSc correction
    SCSc_output = img_plus_ic_mask2$expression(
      "((image * (cosB * cosZ + cvalue)) / (ic + cvalue))",
      list(
        image = img_plus_ic_mask2$select(band),
        ic = img_plus_ic_mask2$select('IC'),
        cosB = img_plus_ic_mask2$select('cosS'),
        cosZ = img_plus_ic_mask2$select('cosZ'),
        cvalue = out_c
      )
    )


    # If NULLs introduced then return original (happens when all cells are NA -
    # i.e. 100% cloud cover.)
    SCSc_output = ee$Algorithms$If(nVals$neq(nNotNullVals),
                                   img_plus_ic_mask2$select(band),
                                   SCSc_output)

    return(SCSc_output)

  })

  # everything about this section feels not very R ish - not a fan but it works.
  # can we replace with $map?
  SCSccorr <- apply_SCSccorr(bandList[1])
  SCSccorr <- ee$Image(SCSccorr)


  for (i in bandList[2:length(bandList)]){
    x <- apply_SCSccorr(i)
    x <-  ee$Image(x)
    SCSccorr <- SCSccorr$addBands(x)
  }

  img_SCSccorr <- SCSccorr$
    addBands(img_plus_ic$select('IC'))

  # return(ee$Image(SCSccorr))
  bandList_IC <- ee$List(c(bandList, 'IC'))$
    flatten()

  img_SCSccorr = img_SCSccorr$unmask(img_plus_ic$select(bandList_IC))$select(bandList)

  img_SCSccorr <- img_SCSccorr$
  select(bandList)$
  addBands(nonCorrectBands)$
  setMulti(props)$
  set('system:time_start', st)

  return (
    ee$Image(img_SCSccorr)
  )
}
