

#Function to convert from dB
toNatural <- function (img, bands=c('VV', 'VH')){
  selected <- img$select(bands)

  adjusted <- ee$Image(10.0)$
    pow(selected$divide(10.0))

  no_speck_bands = anti_selection(img, bands)

  adjusted <- adjusted$addBands(no_speck_bands)$
    copyProperties(img,c('system:time_start','system:time_end'))

  return(ee$Image(adjusted))
}


todB <- function(img, bands=c('VV', 'VH')) {
  selected <- img$select(bands)

  adjusted <- ee$Image(10.0)$
    multiply(selected$log10())

  no_speck_bands = anti_selection(img, bands)

  adjusted <- adjusted$addBands(no_speck_bands)$
    copyProperties(img,c('system:time_start','system:time_end'))
  return(ee$Image(adjusted))
}


# Speckle filter - Reduce noise


#' focalSpeckle
#'
#'speckle filter for SAR imagery.
#'
#' @param img and ee image
#' @param SMOOTHING_RADIUS default 35.
#' @param bands The bands over which to undertake smoothing.
#' @param focal.type default 'mean'. The function to use over the focal area.
#'
#' @return
#' @export
#'
#' @examples
 focalSpeckle <- function(img,  SMOOTHING_RADIUS = 35, bands=c('VV', 'VH'),
                          focal.type=c("mean", "median", "mode", "min", "max")){

   selected <- img$select(bands)

   f_smooth <- function(x){
     if (focal.type[1] == "mean"){
       x <- x$focal_mean(SMOOTHING_RADIUS,'circle','meters')
     } else if (focal.type[1] == "median"){
       x <- x$focal_median(SMOOTHING_RADIUS,'circle','meters')
     } else if (focal.type[1] == "mode"){
       x <- x$focal_mode(SMOOTHING_RADIUS,'circle','meters')
     } else if (focal.type[1] == "max"){
       x <- x$focal_max(SMOOTHING_RADIUS,'circle','meters')
     } else if (focal.type[1] == "min"){
       x <- x$focal_min(SMOOTHING_RADIUS,'circle','meters')
     } else {
       stop(paste0("'", focal.type[1], "'", "must be one of:",
                   "'mean', 'median', 'mode', 'min', max"))
     }
     return(x)
   }

   smoothed <- f_smooth(selected) #Apply a focal mean filter

   no_speck_bands = anti_selection(img, bands)

   smoothed <- smoothed$addBands(no_speck_bands)$
     copyProperties(img,c('system:time_start','system:time_end'))

   return(ee$Image(smoothed))
}


#' s1_add_ratio
#'
#' Add band ratio values to image.
#'
#' @param img an ee image
#' @param method "DB" when image values have been converted to dB using, "DN"
#' otherwise.
#' @param bands default is c('VV', 'VH'). must be length 2. The bands to used to
#' calculate ratio. First value is either subtracted from or used to divide the
#' second dependng on the method arg.
#'
#' @return
#' @export
#'
#' @examples
 s1_add_ratio = function(img, method= c("DB", "DN"), bands=c('VV', 'VH')){

   VV <- img$select(bands[1])

   VH <- img$select(bands[2])

   if (method[1]=="DB"){
     crossPol = VH$subtract(VV)$rename('Ratio')
   } else if (method[1]=="DN"){
     crossPol = VH$divide(VV)$rename('Ratio')
   } else {
     stop("method not supported - must be either 'DB' or 'DN'")
   }

  stack <- img$addBands(crossPol)$
    copyProperties(img,c('system:time_start','system:time_end'))
  return(ee$Image(stack))
}


 #' s1_collect
 #'
 #' @param aoi
 #' @param start.date
 #' @param end.date
 #' @param orbit.pass Default 'DESCENDING'. Must be either
 #'
 #' @return
 #' @export
 #'
 #' @examples
 s1_collect <- function(aoi, start.date, end.date,
                        orbit.pass=c('DESCENDING', 'ASCENDING'),
                        dataset=c("S1_GRD", "S1_GRD_FLOAT")){

    if (!orbit.pass[1] %in% c('DESCENDING', 'ASCENDING')){
      stop(paste0("Orbit pass must be either 'DESCENDING' or 'ASCENDING' not '",
                  orbit.pass, "'"))
    }

   aoi = sf_ext_as_ee(aoi)

   collection <- ee$ImageCollection(paste0('COPERNICUS/', dataset[1]))$
     filterDate(start.date, end.date)

   collection_Full = collection$
     filter(ee$Filter$eq('instrumentMode', 'IW'))$
     filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VV'))$
     filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VH'))$
     #$filter(ee$Filter$eq('relativeOrbitNumber_start', 164))
     filter(ee$Filter$eq('orbitProperties_pass',  orbit.pass[1]))$
     # filterDate(start.date, end.date)$
     filterMetadata('resolution_meters', 'equals' , 10)$
     filterBounds(aoi)

   if (length(collection_Full$getInfo()$features)>0){
     return(collection_Full)
   } else {
     stop(paste0("The collection is empty - perhaps try a different ",
                 "orbital.pass value or change the start.date and ",
                 "end.date values"))
   }

 }





