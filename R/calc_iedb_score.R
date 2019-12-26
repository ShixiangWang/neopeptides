#' Calculate IEDB Score for Peptides
#'
#' @param pep a vector of peptides.
#' @param db an available database, can be one of "human" and "mouse".
#' Or you can set this to specify the (fasta) database file
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
                            db = "human",
                            tmp_dir = file.path(tempdir(), "neopeptides"),
                            clean_tmp = TRUE) {
  stopifnot(has_program("blastp"),
            length(db) == 1)

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
    pep = pep,
    id = names(pep)
  )

  tmp_fasta <- file.path(tmp_dir, "iedb_score.fa")
  make_fasta_file(pep, tmp_fasta)

  message("=> Running blastp for homology to IEDB antigens..")
  tmp_iedb_out <- file.path(tmp_dir, "blastp_iedbout.csv")

  # TODO
  db <- make_blast_db()
  # Run blast
  cmd_blastp(tmp_fasta, db, tmp_iedb_out)

  blastdt <- list.files(
    path = tmp_dir,
    pattern = "blastp_iedbout\\.csv",
    full.names = TRUE
  )

  if (length(blastdt) == 0) {
    message("=> No blast output against database returned!")
    return(data.table::data.table(peptide = pep, iedb_score = 0))
  }

  if (all(file.info(blastdt)$size == 0)) {
    message("=> No database matches found by blast!")
    return(data.table::data.table(peptide = pep, iedb_score = 0))
  }

  blastdt = read_blast_result(blastdt)

  if (nrow(blastdt) == 0) {
    message(paste("=> No IEDB matches found with cannonical AAs, can't compute IEDB score...."))
    return(data.table::data.table(peptide = pep, iedb_score = 0))
  }

  message("=> Summing IEDB local alignments...")
  blastdt[, SW := SW_align(nmer, WT.peptide)]
  message("=> Done.")

  blastdt[, score := SW %>% modeleR(), by = "id"]

  # get full IEDB ref here
  fa <- Biostrings::readAAStringSet(db %>%
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
    "id",
    "nmer",
    "IEDB_anno",
    "score",
    "SW"
  )]

  blastdt <- blastdt[, .SD %>% unique(), .SDcols = c("id", "score")]
  sdt[, id := as.character(id)]
  blastdt[, id := as.character(id)]

  sdt <- merge(sdt, blastdt, by = "id")
  sdt %>% data.table::setnames("pep", "nmer")

  anndt <- anndt[, msw := max(SW), by = "id"] %>%
    .[SW == msw] %>%
    .[, .SD %>% unique(), .SDcols = c("id", "IEDB_anno")]

  # merge equally good IEDB_annos into one
  anndt[, IEDB_anno := paste(IEDB_anno %>% unique(), collapse = "|"), by = "id"]
  anndt <- anndt %>% unique()
  anndt[, id := as.character(id)]
  sdt <- merge(sdt, anndt, by = "id", all.x = TRUE)
  sdt <- sdt[, .SD %>% unique(), .SDcols = c("nmer", "score", "IEDB_anno")]
  return(sdt)
}

# Internal function to Smith-Waterman align two vectors of peptides.
SW_align <- function(pep1,
                     pep2,
                     gap_open = -11L,
                     gap_extend = -1L) {
  al <- Biostrings::pairwiseAlignment(pep1, pep2,
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


utils::globalVariables(
  c(
    ".", "IEDB_anno", "SW", "WT.peptide", "score",
    "msw", "nmer", "id"
  )
)
