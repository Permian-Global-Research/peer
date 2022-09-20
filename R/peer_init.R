#' Initialize ee module to use {peer}
#'
#' A wrapper for rgee::ee_Initialize(). Loads rgee and initializes the ee
#' python module.
#'
#' @param ...
#'
#' @return
#' @export
#'
peer_init <- function(...){
  library(rgee)
  rgee::ee_Initialize(...)
}
