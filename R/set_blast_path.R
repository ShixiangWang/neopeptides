#' Set Path to Blast
#'
#' At default, configure file is '~/.neopeptide/config.yml', you
#' can set by `options(neopeptides.config="new_path.yml")`.
#' Please remember, if you change the configure file, you have
#' to set it per R session.
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
#' @seealso [find_path], [install_database]
set_blast_path <- function(path) {
  if (.Platform$OS.type == "unix") {
    save_to_config("blastp", file.path(path, "blastp"))
    save_to_config("makeblastdb", file.path(path, "makeblastdb"))
  } else {
    save_to_config("blastp", file.path(path, "blastp.exe"))
    save_to_config("makeblastdb", file.path(path, "makeblastdb.exe"))
  }
}
