#' Set Path to Blast
#'
#' @param path a directory containing executable Blast program.
#' Installing it by conda is recommended.
#' @param store if `TRUE` (default), save the path to
#' config file so that you don't need to set this again when
#' start a new R session.
#'
#' @return Nothing
#' @export
#'
#' @examples
#' \donttest{
#' set_blast_path("/Users/wsx/anaconda3/bin/")
#' }
set_blast_path <- function(path, store = TRUE) {
  set_path(path, program = "blast", store = store)
}
