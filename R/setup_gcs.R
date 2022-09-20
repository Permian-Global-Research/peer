
#' Set up GOogle Cloud Storage
#'
#' @param SaK_file
#' @param user
#' @param create_bucket
#' @param bucket_name
#'
#' @return
#' @export
#'
#' @examples
setup_gcs <- function(sak_file, user, check_bucket=TRUE, bucket_name="rgee_store"){
  # Assign the SaK to a EE user.
  rgee::ee_utils_sak_copy(
    sakfile =  sak_file,
    users = user # Unlike GD, we can use the same SaK for multiple users.
  )

  if (isTRUE(check_bucket)){
    peer_init(user=user)
    # # Create your own container
    # project_id <- ee_get_earthengine_path() %>%
    #   list.files(., "\\.json$", full.names = TRUE) %>%
    #   jsonlite::read_json() %>%
    #   '$'(project_id) # Get the Project ID
    #
    # googleCloudStorageR::gcs_create_bucket(bucket_name, projectId = project_id)
    ee_utils_sak_validate(sak_file, bucket = bucket_name)
  }

}


#' Save an ee image to Google Cloud Stoarge
#'
#' @param x an ee image
#' @param filename filename without extension for GCS
#' @param scale
#' @param aoi
#' @param description
#' @param bucket
#' @param max.wait
#' @param task_time
#'
#' @return
#' @export
#'
#' @examples
save_to_gcs <- function(x, filename, scale, region, description, bucket= "rgee_store",
                        max.wait = 30, task_time=10, public = TRUE, add.vsi=TRUE){

  n.attempts <- round(max.wait*60)/task_time

  .d <- ee_image_to_gcs(
    x,
    scale = scale,
    region =  region,
    description = description,
    bucket = bucket,
    fileNamePrefix = filename,
    timePrefix = FALSE,
    fileFormat = "GEO_TIFF",
    formatOptions = list(cloudOptimized = TRUE)
  )

  .d$start()
  ee_monitoring(.d, task_time=task_time, max_attempts=n.attempts)

  .url <- paste0(
    googleCloudStorageR::gcs_download_url(object_name = filename,
                                                bucket = bucket,
                                                public = public),
    ".tif")

  if (isTRUE(add.vsi)){
    .url <- paste0("/vsicurl/",.url)
  }

  return(.url)

}
