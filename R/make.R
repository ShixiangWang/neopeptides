make_fasta_file = function(x, path) {
  AA <- Biostrings::AAStringSet(x)
  Biostrings::writeXStringSet(AA, file = path, format = "fasta")
}

make_blast_db = function() {
  if (is.null(db_file)) {
    if (db == "mouse") {
      db_file <- system.file(
        "extdata", "Mu_iedb.fasta",
        package = "neopeptides",
        mustWork = TRUE
      )
    }
    if (db == "human") {
      db_file <- system.file(
        "extdata", "iedb.fasta",
        package = "neopeptides",
        mustWork = TRUE
      )
    }
  } else {
    if (!file.exists(db_file)) {
      stop("File ", db_file, " does not exist!")
    }
  }

}
