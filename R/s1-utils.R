

# Script developed to obtain the Sentinel-1 annual composites for Kuamut project Area
# Sentinel-1







############/  MAIN ##########

  # Functions

# Terrain Correction

# Implementation by Andreas Vollrath (ESA), inspired by Johannes Reiche (Wageningen)
terrainCorrection <- function(img) {
   imgGeom <- img$geometry()
   srtm <- ee$Image('USGS/SRTMGL1_003')$clip(imgGeom) # 30m srtm
   sigma0Pow <- ee$Image$constant(10)$pow(img$divide(10.0))

  # Article ( numbers relate to chapters)
  # 2$1$1 Radar geometry
   theta_i = img$select('angle')
   phi_i = ee$Terrain$aspect(theta_i)$
     reduceRegion(ee$Reducer$mean(), theta_i$get('system:footprint'), 1000)$
     get('aspect')

  # 2$1$2 Terrain geometry
   alpha_s = ee$Terrain$slope(srtm)$select('slope')
   phi_s = ee$Terrain$aspect(srtm)$select('aspect')

  # 2$1$3 Model geometry
  # reduce to 3 angle
   phi_r = ee$Image$constant(phi_i)$subtract(phi_s)

  # convert all to radians
   phi_rRad = phi_r$multiply(pi / 180)
   alpha_sRad = alpha_s$multiply(pi / 180)
   theta_iRad = theta_i$multiply(pi / 180)
   ninetyRad = ee$Image$constant(90)$multiply(pi / 180)

  # slope steepness in range (eq$ 2)
   alpha_r = (alpha_sRad$tan()$multiply(phi_rRad$cos()))$atan()

  # slope steepness in azimuth (eq 3)
   alpha_az = (alpha_sRad$tan()$multiply(phi_rRad$sin()))$atan()

  # local incidence angle (eq$ 4)
   theta_lia = (alpha_az$cos()$multiply((theta_iRad$subtract(alpha_r))$cos()))$acos()
   theta_liaDeg = theta_lia$multiply(180 / pi)
  # 2$2
  # Gamma_nought_flat
   gamma0 = sigma0Pow$divide(theta_iRad$cos())
   gamma0dB = ee$Image$constant(10)$multiply(gamma0$log10())
   ratio_1 = gamma0dB$select('VV')$subtract(gamma0dB$select('VH'))

  # Volumetric Model
   nominator = (ninetyRad$subtract(theta_iRad)$add(alpha_r))$tan()
   denominator = (ninetyRad$subtract(theta_iRad))$tan()
   volModel = (nominator$divide(denominator))$abs()

  # apply model
   gamma0_Volume = gamma0$divide(volModel)
   gamma0_VolumeDB = ee$Image$constant(10)$multiply(gamma0_Volume$log10())

  # we add a layover/shadow maskto the original implmentation
  # layover, where slope > radar viewing angle
   alpha_rDeg = alpha_r$multiply(180 / pi)
   layover = alpha_rDeg$lt(theta_i)

  # shadow where LIA > 90
   shadow = theta_liaDeg$lt(85)

  # calculate the ratio for RGB vis
   ratio = gamma0_VolumeDB$select('VV')$subtract(gamma0_VolumeDB$select('VH'))

   output = gamma0_VolumeDB$
     addBands(ratio)$
     addBands(alpha_r)$
     addBands(phi_s)$
     addBands(theta_iRad)$
     addBands(layover)$
     addBands(shadow)$
     addBands(gamma0dB)$
     addBands(ratio_1)

  img$addBands(output$select(c('VV', 'VH'), c('VV', 'VH')),NULL,TRUE)
}



#Function to convert from dB
toNatural <- function (img){
  ee$Image(10.0)$
    pow(img$select(0)$
          divide(10.0))$
    copyProperties(img,c('system:time_start','system:time_end'))
}


todB <- function(img) {
  ee$Image(10.0)$
    multiply(img$select(0)$
               log10())$
    copyProperties(img,c('system:time_start','system:time_end'))
}


# Speckle filter - Reduce noise


 focalSpeckle = function(img,  SMOOTHING_RADIUS = 35){
   selected = img
   smoothed = selected$focal_mean(SMOOTHING_RADIUS,'circle','meters') #Apply a focal mean filter

   smoothed$copyProperties(img,c('system:time_start','system:time_end'))
}


# Filter and Stack Collection

 stackcol = function(img){

   VV = img$select('VV')
  # smoothVV = VV$focal_mean(SMOOTHING_RADIUS,'circle','meters')

   VH = img$select('VH')
  # smoothVH = VH$focal_mean(SMOOTHING_RADIUS,'circle','meters')

   crossPol = VH$subtract(VV)
  # smoothcrossPol = crossPol$focal_mean(SMOOTHING_RADIUS,'circle','meters')

   stack = VV$addBands(VH)
  stack = stack$addBands(crossPol)
  stack = stack$rename(c('VV','VH','Ratio'))

  stack$copyProperties(img,c('system:time_start','system:time_end'))
}


 #' s1_collect
 #'
 #' @param aoi
 #' @param start.date
 #' @param end.date
 #'
 #' @return
 #' @export
 #'
 #' @examples
 s1_collect <- function(aoi, start.date, end.date ){ #min.cloud=90

   aoi = rgee::sf_as_ee(aoi)

   #Pre-processing Functions

   # create function to crop with 'Subset' boundaries
   subset_bounds <- function(img) {
     # Crop by table extension
     img$clip(aoi)$
       copyProperties(img,c('system:time_start','system:time_end'))
   }

   # Load Sentinel-1 C-band SAR Ground Range collection (log scale, VV, descending)
   collection_Full = ee$ImageCollection('COPERNICUS/S1_GRD')$
     filter(ee$Filter$eq('instrumentMode', 'IW'))$
     filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VV'))$
     filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VH'))$
     #$filter(ee$Filter$eq('relativeOrbitNumber_start', 164))
     filter(ee$Filter$eq('orbitProperties_pass', 'DESCENDING'))$
     filterDate(start.date, end.date)$
     filterMetadata('resolution_meters', 'equals' , 10)$
     filterBounds(aoi)


   # Clipping Collections
   collection <- collection_Full$map(subset_bounds)

   return(collection)
 }

 s1_terrain_correct <- function(x){
   # print(collection, 'Collection_Clip')
   #Appy Terrain Correction

   combine_funcs <- function(img){
     img <- terrainCorrection(img)$
       subset_bounds()
     illuminationCorrection(img, dem, buf)
   }

   if (inherits(x, "ee.image.Image" )){
     return(combine_illumes(x))
   } else if (inherits(x, "ee.imagecollection.ImageCollection")){
     return(x$map(combine_illumes)) #$map(function(x) x$clip(aoi))
   }
   s1_TC <- collection$
     map(terrainCorrection)$
     map(subset_bounds)
   # Ratio calculations

   collection <- collection$map(stackcol)
   # print(collection, 'Collection after ratio')



   # print(s1_TC, 'Terrain Correction S1')

   firstTerrainCorrection <- s1_TC$
     filterDate('2016-01-01', '2016-12-31')$
     select('VV')$
     median()

   # Map$addLayer(firstTerrainCorrection,{min:-30,max:5},"Terrain corrected",0)
   #-----------------------------------------# Filter Collection

   # No TC / Terrain Correction
   # collection = collection
   collection = s1_TC
   # Select band
   band = 'VV'

   # Filter band collection
   collection = collection$select(band)

   #Convert to natural
   collection = collection$map(toNatural)

   #Speckle filter
   collection = collection$map(focalSpeckle)

   #Re-Convert to dB
   collection = collection$map(todB)

   #S1 Collection 2015

   s1.out = collection$filterDate('2015-01-01', '2015-12-31')$median()
   return(s1.out)
 }





# Map$addLayer(s12015, {min:-15,max:0}, 'Sentinel1  2015',0)

# #S1 Collection 2016
#
#  s12016 = collection$filterDate('2016-01-01', '2016-12-31')$median()
# # Map$addLayer(s12016, {min:-15,max:0}, 'Sentinel1  2016',0)
#
# #S1 Collection 2017
#
#  s12017 = collection$filterDate('2017-01-01', '2017-12-31')$median()
# # Map$addLayer(s12017, {min:-15,max:0}, 'Sentinel1  2017',0)
#
# #S1 Collection 2018
#
#  s12018 = collection$filterDate('2018-01-01', '2018-12-31')$median()
# # Map$addLayer(s12018, {min:-15,max:0}, 'Sentinel1  2018',0)
#
# #S1 Collection 2019
#
#  s12019 = collection$filterDate('2019-01-01', '2019-12-31')$median()
# # Map$addLayer(s12019, {min:-15,max:0}, 'Sentinel1  2019',0)
#
# #S1 Collection 2020
#
#  s12020 = collection$filterDate('2020-01-01', '2020-12-31')$median()
# # Map$addLayer(s12020, {min:-15,max:0}, 'Sentinel1  2020',0)
#
#
# #S1 Collection 2021
#
#  s12021 = collection$filterDate('2021-01-01', '2021-12-31')$median()
# # Map$addLayer(s12021, {min:-15,max:0}, 'Sentinel1  2021',0)
#
#
# #S1 Collection 2022
#
#  s12022 = collection$filterDate('2022-01-01', '2022-07-31')$median()
# Map$addLayer(s12022, {min:-15,max:0}, 'Sentinel1  2022',0)


#-------------------- DATA EXPORT (MAPS and Validation)

#
# # Export a cloud-optimized GeoTIFF$ (2016)
# Export$img$toDrive({
#   img: s12015,
#   description: 'S1_'+band+'_2015_RTC',
#   folder:'L2 - Radiometric Terrain Correction',
#   scale: 10,
#   region: area,
#   fileFormat: 'GeoTIFF',
#   maxPixels: 1e13
# })
#
# /*
#   # Export a cloud-optimized GeoTIFF$ (2016)
# Export$img$toDrive({
#   img: s12016,
#   description: 'S1_'+band+'_2016_RTC',
#   folder:'L2 - Radiometric Terrain Correction',
#   scale: 10,
#   region: area,
#   fileFormat: 'GeoTIFF',
#   maxPixels: 1e13
# })
#
# # Export a cloud-optimized GeoTIFF$ (2017)
# Export$img$toDrive({
#   img: s12017,
#   description: 'S1_'+band+'_2017_RTC',
#   folder:'L2 - Radiometric Terrain Correction',
#   scale: 10,
#   region: area,
#   fileFormat: 'GeoTIFF',
#   maxPixels: 1e13
# })
#
# # Export a cloud-optimized GeoTIFF$ (2018)
# Export$img$toDrive({
#   img: s12018,
#   description: 'S1_'+band+'_2018_RTC',
#   folder:'L2 - Radiometric Terrain Correction',
#   scale: 10,
#   region: area,
#   fileFormat: 'GeoTIFF',
#   maxPixels: 1e13
# })
#
# # Export a cloud-optimized GeoTIFF$ (2019)
# Export$img$toDrive({
#   img: s12019,
#   description: 'S1_'+band+'_2019_RTC',
#   folder:'L2 - Radiometric Terrain Correction',
#   scale: 10,
#   region: area,
#   fileFormat: 'GeoTIFF',
#   maxPixels: 1e13
# })
#
# # Export a cloud-optimized GeoTIFF$ (2020)
# Export$img$toDrive({
#   img: s12020,
#   description: 'S1_'+band+'_2020_RTC',
#   folder:'L2 - Radiometric Terrain Correction',
#   scale: 10,
#   region: area,
#   fileFormat: 'GeoTIFF',
#   maxPixels: 1e13
# })
#
# */
#
#   # Export a cloud-optimized GeoTIFF$ (2021)
# Export$img$toDrive({
#   img: s12021,
#   description: 'S1_'+band+'_2021_RTC',
#   folder:'L2 - Radiometric Terrain Correction',
#   scale: 10,
#   region: area,
#   fileFormat: 'GeoTIFF',
#   maxPixels: 1e13
# })
#
#
# # Export a cloud-optimized GeoTIFF$ (2022)
# Export$img$toDrive({
#   img: s12022,
#   description: 'S1_'+band+'_2022_RTC',
#   folder:'L2 - Radiometric Terrain Correction',
#   scale: 10,
#   region: area,
#   fileFormat: 'GeoTIFF',
#   maxPixels: 1e13
# })
#





