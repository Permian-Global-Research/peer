#' MapAll
#'
#' Map all the bands of an ee image
#'
#' @param img an ee image
#' @param pal default is `viridisLite::plasma(256)`. HEX code colours...
#'
#' @return
#' @export
#'
#' @examples
MapAll <- function(img, zoom=11, pal=viridisLite::plasma(256), opacity=1,
                   .n = 150){
  .names <- img_band_names(img)


  ee_fc <- rando_points(img, .n)

  s <- point_extract(img, ee_fc)
  mins <- dplyr::summarise_all(s, min)
  maxs <- dplyr::summarise_all(s, max)
  # mins <- rgeeExtra::ee_minValue(img$reproject(ee$Projection('EPSG:4326'), NULL, 1),
  #                                mode = "Points",
  #                                sample_size = 1000)
  # maxs <- rgeeExtra::ee_maxValue(img)

  .tab <- rbind(mins, maxs) |>
    as.data.frame()

  l <- lapply(seq_along(.tab), function(x){
    list(bands=.names[x],
         min=.tab[1, x],
         max=.tab[2, x],
         palette=pal,
         opacity = opacity)
  } )

  .cent <- img$geometry()$
    bounds()$
    centroid(10)$
    getInfo()$
    coordinates

  rgee::Map$setCenter(lon = .cent[1] , lat = .cent[2], zoom = zoom)

  m <- rgee::Map$addLayer(img, l[[1]],  l[[1]]$bands)

  if (length(l)>1){
    for (i in l[2:length(l)]){
      m <- m +
        rgee::Map$addLayer(img, i,  i$bands)
    }
  }


  return(m)
}



g.inf <- function(img){
  img$geometry()$getInfo()
}
