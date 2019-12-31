#' Install IEDB and Proteome Reference Database
#'
#' At default, configure file is '~/.neopeptide/config.yml', you
#' can set by `options(neopeptides.config="new_path.yml")`.
#' Please remember, if you change the configure file, you have
#' to set it per R session.
#'
#' @param db_path a path to store the databases.
#' @param force if `TRUE`, force to re-download and re-install.
#' @param species can be 'human' and 'mouse'. Default select them all.
#' @param data_type can be 'IEDB' and 'Proteome'. Default select them all.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' \dontrun{
#' # Default all databases will be downloaded and installed.
#' install_database()
#' }
#' @seealso [set_blast_path], [calc_iedb_score], [calc_dissimilarity]
install_database <- function(db_path = "~/.neopeptide/db",
                             force = FALSE,
                             species = c("human", "mouse"),
                             data_type = c("IEDB", "Proteome")) {
  species <- match.arg(species, choices = c("human", "mouse"), several.ok = TRUE)
  data_type <- match.arg(data_type, choices = c("IEDB", "Proteome"), several.ok = TRUE)

  if (!dir.exists(db_path)) {
    dir.create(db_path, recursive = TRUE)
  }

  proteome_human <- "https://github.com/XSLiuLab/neopeptide_data/raw/master/Homo_sapiens.GRCh38.pep.all.fa.gz"
  proteome_mouse <- "https://github.com/XSLiuLab/neopeptide_data/raw/master/Mus_musculus.GRCm38.pep.all.fa.gz"
  iedb_human <- "https://github.com/XSLiuLab/neopeptide_data/raw/master/iedb.fasta.gz"
  iedb_mouse <- "https://github.com/XSLiuLab/neopeptide_data/raw/master/Mu_iedb.fasta.gz"

  proteome <- c(proteome_human, proteome_mouse)
  iedb <- c(iedb_human, iedb_mouse)
  names(proteome) <- names(iedb) <- c("human", "mouse")

  database <- c()
  if ("IEDB" %in% data_type) {
    database <- c(database, iedb[species])
  }
  if ("Proteome" %in% data_type) {
    database <- c(database, proteome[species])
  }

  # Download database
  local_file <- file.path(db_path, basename(database))
  local_fa_file <- sub("\\.gz", "", local_file)
  for (i in seq_along(database)) {
    if (!force & (file.exists(local_file[i]) | file.exists(local_fa_file[i]))) {
      message("=> File ", local_file[i], " or ", local_fa_file[i], " exist, skipping...")
    } else {
      message("=> Start downloading...")
      utils::download.file(url = database[i], destfile = local_file[i])
    }
  }

  # Unzip gz files
  for (i in seq_along(local_file)) {
    if (file.exists(local_file[i])) {
      message("=> Unzipping ", local_file[i])
      R.utils::gunzip(local_file[i], remove = TRUE)
    } else {
      if (file.exists(local_fa_file[i])) {
        message("=> File ", local_fa_file[i], " exist, skipping...")
      } else {
        stop("Some bad thing happened, please take a check or report to the developer!")
      }
    }
  }

  # Make blast database and store the path
  #
  # if files with a same prefix > 1, we think a database exists
  # otherwise, we make it with blast command
  for (i in seq_along(local_fa_file)) {
    files <- list.files(
      path = db_path,
      pattern = paste0("^", basename(local_fa_file[i])),
      full.names = TRUE
    )
    if (!force && length(files) > 1) {
      message("=> Database ", local_fa_file[i], " found, skipping...")
    } else {
      # Call makeblastdb
      message("=> Making blastp database ", local_fa_file[i], "...")
      cmd_blastdb(local_fa_file[i])
    }
    # Store database path to config file
    key <- paste0(
      "db_", names(database)[i], "_",
      ifelse(grepl("iedb", basename(local_fa_file[i])),
        "iedb", "proteome"
      )
    )
    save_to_config(key, local_fa_file[i])
  }
}
