

library(sf)

aoi <- read_sf("data/Kuamut_boundaries_dissolved.fgb") %>%
  st_transform(4326)
aoi_ee <- sf_as_ee(aoi)
aoi_cent <- st_centroid(aoi) %>%
  st_coordinates()



subset_bounds <- function(img) {
  # // Crop by table extension
  img$clip(aoi_ee)$
  copyProperties(img,c('system:time_start','system:time_end'))
}


img_c <- ee$ImageCollection('LANDSAT/LC08/C01/T1_SR')$
  filterBounds(aoi_ee)$
  filter(ee$Filter$lt('CLOUD_COVER', 25))$
  map(subset_bounds)$
  filterDate('2015-01-01','2015-08-31')

ee_print(img_c)
# img_c2 <- img_c$filter(ee$Filter$lt('CLOUD_COVER', 25))#$
  # filterMetadata('index','equals','LC08_117056_20150711')
  # filter(ee$Filter('system:index', "LC08_117056_20150711"))

ee_print(img_c2)


getQABits <- function(img, start, end, newName) {
  # // Compute the bits we need to extract.
  pattern <- sapply(c(start:end), function(x) 2 ^x) |>
    sum()
  # // Return a single band image of the extracted QA bits, giving the band
  # // a new name.
  img$select(0, newName)$
   bitwiseAnd(pattern)$
  rightShift(start)
};

# // A function to mask out cloudy pixels.
cloud_shadows = function(img) {
  # // Select the QA band.
  QA <- img$select('pixel_qa')
  # // Get the internal_cloud_algorithm_flag bit.
  getQABits(QA, 3,3, 'cloud_shadows')$eq(0)
  # // Return an image masking out cloudy areas.
}

# // A function to mask out cloudy pixels.
clouds <- function(img) {
  # // Select the QA band.
  QA <- img$select('pixel_qa');
  # // Get the internal_cloud_algorithm_flag bit.
  getQABits(QA, 5,5, 'Cloud')$eq(0);
  # // Return an image masking out cloudy areas.
};

maskClouds <- function(img){
  cl.sh <- cloud_shadows(img)
  cl <- clouds(img)
  img = img$updateMask(cl.sh)
  img$updateMask(cl)
}

composite_free <- img_c$map(maskClouds)

img_free = img_c$median()

# img_c.info <- img_c$getInfo()
# img_c2.info <- img_c2$getInfo()

# i <- ee$Image('LANDSAT/LC08/C01/T1_TOA/LC08_044034_20140318')
#
# id <- i$getString('system:index')
# id$getInfo()


# img_c2 <-img_c2$map(function(img){
#   index <- img$getString('system:index')
#   return(index)
# })


# c2 <- img_c2$median()
  # reduceRegion(
  #   reducer = ee$Reducer$mean(),
  #   geometry = aoi_ee # Geometry is specified here.
  #   # scale = 30,
  #   # maxPixels = 1e10
  # )

vizParams <- list(
  bands = c('B4', 'B3', 'B2'),
  min = 0,
  max = 1000,
  gamma = c(0.95, 1.1, 1)
)

# Center the map and display the image.
Map$setCenter(lon = aoi_cent[1], lat = aoi_cent[2], zoom = 12) # San Francisco Bay
Map$addLayer(img_free$select("pixel_qa"), list(max = 1000, min = 0), name = "pixel_qa") +
Map$addLayer(img_free, vizParams, 'false color composite')

# srtm <- ee$Image("CGIAR/SRTM90_V4")
# ee_print(srtm)






img_c2 <-img_c$map(function(img){
  index <- img$getString('system:index')
  include <- img$
    getNumber("CLOUD_COVER")$lt(25)#$
    # and(index$match('LC08_117056_20150711')$
    #       size()$
    #       eq(0))
  return(ee$Algorithms$If(include, img, NA))
})

img_c2$getInfo()

l8_collection <- function(){

  sr14 <-  ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
  .filterBounds(Kuamut)
  .map(subset_bounds)
  .filterDate('2015-01-01','2015-08-31')
  # //.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',100))
  # //.sort('CLOUD_COVER',false);

  print(sr14)

  var sr14 = sr14.map(function (image) {
    var index = image.getString('system:index')
    var include = image.getNumber("CLOUD_COVER").lt(25)
    .and(index.match('LC08_117056_20150711').size().eq(0))


    return ee.Algorithms.If(include, image, null)
  }, true)
}


# Load an image.
landsat <- ee$Image('LANDSAT/LC08/C01/T1_TOA/LC08_044034_20140318')

# Define the visualization parameters.
vizParams <- list(
  bands = c('B5', 'B4', 'B3'),
  min = 0,
  max = 0.5,
  gamma = c(0.95, 1.1, 1)
)

# Center the map and display the image.
Map$setCenter(lon = -122.1899, lat = 37.5010, zoom = 10) # San Francisco Bay
Map$addLayer(landsat, vizParams, 'false color composite')
