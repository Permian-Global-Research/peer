l8_correct_terrain <-  function(collection, dem) {

  # ////////////////////////////////////////////////////////////////////////////////
  # // Function to calculate illumination condition (IC)$ Function by Patrick Burns and Matt Macander
  illuminationCondition <- function(img){

    # // Extract image metadata about solar position
     SZ_rad = ee$Image$constant(ee$Number(img$get('SOLAR_ZENITH_ANGLE')))$multiply(pi)$divide(180)$clip(img$geometry()$buffer(10000))
     SA_rad = ee$Image$constant(ee$Number(img$get('SOLAR_AZIMUTH_ANGLE'))$multiply(pi)$divide(180))$clip(img$geometry()$buffer(10000))
    # // Creat terrain layers
     slp = ee$Terrain$slope(dem)$clip(img$geometry()$buffer(10000))
     slp_rad = ee$Terrain$slope(dem)$multiply(pi)$divide(180)$clip(img$geometry()$buffer(10000))
     asp_rad = ee$Terrain$aspect(dem)$multiply(pi)$divide(180)$clip(img$geometry()$buffer(10000))

    # // Calculate the Illumination Condition (IC)
    # // slope part of the illumination conditions
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
                                              list('sinZ'= sinZ,
                                                'sinS'= sinS,
                                                'cosAziDiff'= cosAziDiff))
    # // full illumination condition (IC)
     ic = slope_illumination$add(aspect_illumination)

    # // Add IC to original image
     img_plus_ic = ee$Image(img$addBands(ic$rename('IC'))$addBands(cosZ$rename('cosZ'))$addBands(cosS$rename('cosS'))$addBands(slp$rename('slope')))

    return(img_plus_ic)
  }
  # ////////////////////////////////////////////////////////////////////////////////
  # // Function to apply the Sun-Canopy-Sensor + C (SCSc) correction method to each
  # // image$ Function by Patrick Burns and Matt Macander
  illuminationCorrection <- function(img){
     props = img$toDictionary()
     st = img$get('system:time_start')

     img_plus_ic = img
     mask1 = img_plus_ic$select('nir')$gt(-0.1)
     mask2 = img_plus_ic$select('slope')$gte(5)$
       and(img_plus_ic$select('IC')$gte(0))$
       and(img_plus_ic$select('nir')$gt(-0.1))
     img_plus_ic_mask2 = ee$Image(img_plus_ic$updateMask(mask2))

    # // Specify Bands to topographically correct
     bandList = c('blue','green','red','nir','swir1','swir2') # must be all right?
     compositeBands = img$bandNames()
     nonCorrectBands = img$select(compositeBands$removeAll(bandList))

     geom = ee$Geometry(img$get('system:footprint'))$bounds()$buffer(10000)

     apply_SCSccorr <- function(band){
       method = 'SCSc'
       out = img_plus_ic_mask2$select('IC', band)$reduceRegion(
         list(reducer= ee$Reducer$linearFit(), #// Compute coefficients: a(slope), b(offset), c(b/a)
        geometry= ee$Geometry(img$geometry()$buffer(-5000)), #// trim off the outer edges of the image for linear relationship
        scale= 300,
        maxPixels= 1000000000
      ))

      if (is.nul(out) | is.null(undefined)){ #(out === null || out === undefined ) ## not sure on this bit...
        return(img_plus_ic_mask2$select(band))
      }

      else{
         out_a = ee$Number(out$get('scale'))
         out_b = ee$Number(out$get('offset'))
         out_c = out_b$divide(out_a)
        # // Apply the SCSc correction
         SCSc_output = img_plus_ic_mask2$expression(
          "((image * (cosB * cosZ + cvalue)) / (ic + cvalue))",
          list(
            'image'= img_plus_ic_mask2$select(band),
            'ic'= img_plus_ic_mask2$select('IC'),
            'cosB'= img_plus_ic_mask2$select('cosS'),
            'cosZ'= img_plus_ic_mask2$select('cosZ'),
            'cvalue'= out_c
          ))

        return(SCSc_output)
      }

    }

     img_SCSccorr = ee$Image(bandList$map(apply_SCSccorr))$addBands(img_plus_ic$select('IC'))
     bandList_IC = ee$List(c(bandList, 'IC'))$flatten()
    img_SCSccorr = img_SCSccorr$unmask(img_plus_ic$select(bandList_IC))$select(bandList)

    return (img_SCSccorr$addBands(nonCorrectBands))$
      setMulti(props)$
      set('system:time_start',st)
  }

  collection = collection$map(illuminationCondition)
  collection = collection$map(illuminationCorrection)

  return(collection)

}
