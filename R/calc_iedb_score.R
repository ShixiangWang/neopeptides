#' Calculate IEDB Score for Peptides
#'
#' @param pep a vector of peptides.
#' @param db an available database, can be one of "human" and "mouse".
#' @param db_file you can set this to specify the (fasta) database file
#' to be searched instead of using standard database.
#' @param tmp_dir path for storing temp files.
#' @param clean_tmp if `TRUE`, remove temp directory.
#'
#' @return Data table of IEDB score
#'
#' @import data.table
#' @importFrom data.table :=
#' @export
#' @examples
#' \donttest{
#' calc_iedb_score("MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDP")
#' calc_iedb_score("MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDP", db = "mouse")
#' }
calc_iedb_score <- function(pep,
                            db = c("human", "mouse"),
                            db_file = NULL,
                            tmp_dir = file.path(tempdir(), "neopeptides"),
                            clean_tmp = TRUE) {
  stopifnot(has_program("blastp"))
  db <- match.arg(db)

  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE)
  }

  if (clean_tmp) {
    on.exit({
      message("=> Removing temporary files...")
      try(
        unlink(tmp_dir, recursive = TRUE)
      )
    })
  }

  # generate fastas to query
  names(pep) <- 1:length(pep)
  sdt <- data.table::data.table(
    pep = pep %>% as.character(),
    nmer_id = names(pep)
  )

  AA <- Biostrings::AAStringSet(pep, use.names = TRUE)
  tmp_fasta <- file.path(tmp_dir, "iedb_score_fasta.fa")
  Biostrings::writeXStringSet(AA, file = tmp_fasta, format = "fasta")

  message("=> Running blastp for homology to IEDB antigens..")
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

  db <- db_file
  tmp_iedb_out <- file.path(tmp_dir, "blastp_iedbout.csv")

  cmd_blastp(
    tmp_fasta,
    db,
    tmp_iedb_out
  )

  blastdt <- list.files(
    path = tmp_dir,
    pattern = "iedbout\\.csv",
    full.names = TRUE
  )

  if (length(blastdt) == 0) {
    message("=> No blast output against database returned!")
    return(data.table::data.table(nmer = pep))
  }

  if (all(file.info(blastdt)$size == 0)) {
    message("=> No database matches found by blast!")
    return(data.table::data.table(nmer = pep, iedb_score = 0))
  }

  blastdt <- data.table::fread(blastdt)
  colnames(blastdt) <-
    c(
      "nmer_id",
      "IEDB_anno",
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

  blastdt <- blastdt[, nmer := nmer %>% stringr::str_replace_all(pattern = "-|\\*", replacement = "")] %>%
    .[, WT.peptide := WT.peptide %>% stringr::str_replace_all(pattern = "-|\\*", replacement = "")] %>%
    .[!is.na(nmer) & !is.na(WT.peptide)]

  blastdt <- blastdt[nmer %like% "^[ARNDCQEGHILKMFPSTWYV]+$" & WT.peptide %like% "^[ARNDCQEGHILKMFPSTWYV]+$"]

  if (nrow(blastdt) == 0) {
    message(paste("=> No IEDB matches found with cannonical AAs, can't compute IEDB score...."))
    return(data.table::data.table(nmer = pep))
  }

  message("=> Summing IEDB local alignments...")
  blastdt[, SW := SW_align(nmer, WT.peptide)]
  message("=> Done.")

  blastdt[, iedb_score := SW %>% modeleR(), by = "nmer_id"]

  # get full IEDB ref here
  fa <- Biostrings::readAAStringSet(db_file %>%
    stringr::str_replace("bdb$", "fasta") %>%
    stringr::str_replace("^-db\\ ", ""))
  f <- fa %>% as.character()
  names(f) <- names(fa)

  blastdt[, IEDB_anno := lapply(IEDB_anno, function(i) {
    mv <- f[which(stringr::str_detect(pattern = stringr::fixed(i), names(f)))]
    mv <- mv[which(stringr::str_detect(pattern = stringr::fixed(WT.peptide), mv))]
    return(paste(names(mv), WT.peptide, collapse = "|"))
  }), by = 1:nrow(blastdt)]

  anndt <- blastdt[, .SD %>% unique(), .SDcols = c(
    "nmer_id",
    "nmer",
    "IEDB_anno",
    "iedb_score",
    "SW"
  )]

  blastdt <- blastdt[, .SD %>% unique(), .SDcols = c("nmer_id", "iedb_score")]
  sdt[, nmer_id := as.character(nmer_id)]
  blastdt[, nmer_id := as.character(nmer_id)]

  sdt <- merge(sdt, blastdt, by = "nmer_id")
  sdt %>% data.table::setnames("pep", "nmer")

  anndt <- anndt[, msw := max(SW), by = "nmer_id"] %>%
    .[SW == msw] %>%
    .[, .SD %>% unique(), .SDcols = c("nmer_id", "IEDB_anno")]

  # merge equally good IEDB_annos into one
  anndt[, IEDB_anno := paste(IEDB_anno %>% unique(), collapse = "|"), by = "nmer_id"]
  anndt <- anndt %>% unique()
  anndt[, nmer_id := as.character(nmer_id)]
  sdt <- merge(sdt, anndt, by = "nmer_id", all.x = TRUE)
  sdt <- sdt[, .SD %>% unique(), .SDcols = c("nmer", "iedb_score", "IEDB_anno")]
  return(sdt)
}

# Internal function to Smith-Waterman align two vectors of peptides.
SW_align <- function(col1,
                     col2,
                     gap_open = -11L,
                     gap_extend = -1L) {
  al <- Biostrings::pairwiseAlignment(col1, col2,
    substitutionMatrix = "BLOSUM62",
    gapOpening = gap_open,
    gapExtension = gap_extend,
    type = "local",
    scoreOnly = TRUE
  )

  if (length(al) == 0) al <- as.numeric(NA)

  return(al)
}

modeleR <- function(als, a = 26, k = 4.86936, dislike = FALSE) {
  be <- -k * (a - als)
  sumexp <- sum(exp(be))
  Zk <- 1 + sumexp
  R <- sumexp / Zk
  if (dislike) {
    return(1 - R)
  } else {
    return(R)
  }
}

# makeblastdb -in Mu_iedb.fasta  -dbtype prot -out Mu_iedb.fasta -hash_index

utils::globalVariables(
  c(
    ".", "IEDB_anno", "SW", "WT.peptide", "iedb_score",
    "msw", "nmer", "nmer_id"
  )
)
