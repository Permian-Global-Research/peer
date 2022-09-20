.onLoad <- function(libname, pkgname) {
  message(crayon::bgYellow(crayon::bold(crayon::black(
    "{peer} loaded...
Run `setup_gcs()` to set up your google cloud storage credentials atleast once
and then `peer_init()` at the start of each R session"
  ))))
  # peer_init()
}
