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
peer_init <- function(user = NULL,..., gcs=TRUE){
  rgee::ee_Initialize(user=user, gcs=gcs, ...)
}
