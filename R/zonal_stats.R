
# NOT WORKING YET -  from https://developers.google.com/earth-engine/tutorials/community/extract-raster-values-for-points
zonalStats <- function (ic, fc) {
  # Initialize internal params dictionary$
  .params <- list(
    reducer= ee$Reducer$mean(),
    scale= NULL,
    crs= NULL,
    bands= NULL,
    bandsRename= NULL,
    imgProps= NULL,
    imgPropsRename= NULL,
    datetimeName= 'datetime',
    datetimeFormat= 'YYYY-MM-dd HH:mm:ss'
  )

  # Replace initialized params with provided params$
  # if (params) {
  #   for (param in params) {
  #     .params[param] = params[param] || .params[param]
  #   }
  # }

  # Set default parameters based on an image representative$
   imgRep = ic$first()
   nonSystemImgProps = ee$Feature(NULL)$
     copyProperties(imgRep)$propertyNames()
  if (!.params$bands) .params$bands = imgRep$bandNames()
  if (!.params$bandsRename) .params$bandsRename = .params$bands
  if (!.params$imgProps) .params$imgProps = nonSystemImgProps
  if (!.params$imgPropsRename) .params$imgPropsRename = .params$imgProps

  # Map the reduceRegions function over the image collection$
   results = ic$map(function(img) {
    # Select bands (optionally rename), set a datetime & timestamp property$
    img = ee$Image(img$select(.params$bands, .params$bandsRename))$
      set(.params$datetimeName, img$date()$format(.params$datetimeFormat))$
      set('timestamp', img$get('system:time_start'))

    # Define final image property dictionary to set in output features$
     propsFrom = ee$List(.params$imgProps)$
       cat(ee$List(c(.params$datetimeName, 'timestamp')))

     propsTo = ee$List(.params$imgPropsRename)$
       cat(ee$List(c(.params$datetimeName, 'timestamp')))

     imgProps = img$toDictionary(propsFrom)$
       rename(propsFrom, propsTo)

    # Subset points that intersect the given image$
     fcSub = fc$filterBounds(img$geometry())

    # Reduce the image by regions$
     out <- img$reduceRegions(
       list(
         collection = fcSub,
         reducer = .params$reducer,
         scale = .params$scale,
         crs = .params$crs
       )
     )$
       map(function(f) {# Add metadata to each feature$
      f$set(imgProps)
    })
  })$
     flatten()$
     filter(ee$Filter$notNull(.params$bandsRename))

  return(results)
}
