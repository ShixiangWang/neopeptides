make_fasta_file <- function(x, path) {
  AA <- Biostrings::AAStringSet(x)
  Biostrings::writeXStringSet(AA, file = path, format = "fasta")
}

make_blastp_db <- function(db = c("human", "mouse"),
                           data_type = c("IEDB", "Genome")) {
  stopifnot(length(db) == 1, length(data_type) == 1)
  # Get an installed database or create a database based on a fasta file
  if (db %in% c("human", "mouse")) {
    # db_mouse_iedb
    # db_human_genome
    # db_mouse_genome
    # db_human_iedb
    db_key <- paste0("db_", db, "_", ifelse(data_type == "IEDB", "iedb", "genome"))
    load_from_config(db_key)
  } else {
    if (!file.exists(db)) {
      stop("File ", db, " does not exist!", call. = FALSE)
    }
    # Check if a database has been created
    files <- list.files(
      path = dirname(db),
      pattern = paste0("^", basename(db)),
      full.names = TRUE
    )
    if (length(files) == 1) {
      cmd_blastdb(db)
    }
    return(db)
  }
}
