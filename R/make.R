make_fasta_file = function(x, path) {
  AA <- Biostrings::AAStringSet(x)
  Biostrings::writeXStringSet(AA, file = path, format = "fasta")
}

make_blast_db = function(db = NULL,
                         data_type = c("IEDB", "ALL")) {
  #https://github.com/XSLiuLab/neopeptide_data/raw/master/Homo_sapiens.GRCh38.pep.all.fa.gz
  #https://github.com/XSLiuLab/neopeptide_data/raw/master/Mus_musculus.GRCm38.pep.all.fa.gz
  if (is.null(db)) {
    if (db == "mouse") {
      db <- system.file(
        "extdata", "Mu_iedb.fasta",
        package = "neopeptides",
        mustWork = TRUE
      )
    }
    if (db == "human") {
      db <- system.file(
        "extdata", "iedb.fasta",
        package = "neopeptides",
        mustWork = TRUE
      )
    }
  } else {
    if (!file.exists(db)) {
      stop("File ", db, " does not exist!")
    }
  }

}
