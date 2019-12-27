#' Set Path to Blast
#'
#' @param path a directory containing executable Blast program.
#' Installing it by conda is recommended.
#'
#' @return Nothing
#' @export
#'
#' @examples
#' \dontrun{
#' set_blast_path("/Users/wsx/anaconda3/bin/")
#' }
set_blast_path <- function(path) {
  if (.Platform$OS.type == "unix") {
    save_to_config("blastp", file.path(path, "blastp"))
    save_to_config("makeblastdb", file.path(path, "makeblastdb"))
  } else {
    save_to_config("blastp", file.path(path, "blastp.exe"))
    save_to_config("makeblastdb", file.path(path, "makeblastdb.exe"))
  }
}
