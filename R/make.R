make_fasta_file = function(x, path) {
  AA <- Biostrings::AAStringSet(x)
  Biostrings::writeXStringSet(AA, file = path, format = "fasta")
}
