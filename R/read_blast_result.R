read_blast_result <- function(path) {
  # TODO: read a list of files in a vector
  dt <- data.table::fread(path)
  colnames(dt) <-
    c(
      "id",
      "anno",
      "nmer",
      "q_start",
      "q_stop",
      "WT.peptide",
      "s_start",
      "s_end",
      "overlap_length",
      "mismatch_length",
      "pident",
      "evalue",
      "bitscore"
    )

  dt <- dt[, nmer := nmer %>% stringr::str_replace_all(pattern = "-|\\*", replacement = "")] %>%
    .[, WT.peptide := WT.peptide %>% stringr::str_replace_all(pattern = "-|\\*", replacement = "")] %>%
    .[!is.na(nmer) & !is.na(WT.peptide)]

  dt <- dt[nmer %like% "^[ARNDCQEGHILKMFPSTWYV]+$" & WT.peptide %like% "^[ARNDCQEGHILKMFPSTWYV]+$"]
}
