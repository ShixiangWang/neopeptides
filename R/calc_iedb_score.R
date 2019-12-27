#' Calculate IEDB Score for Peptides
#'
#' @param pep a vector of peptides.
#' @param db an available database, can be one of "human" and "mouse".
#' Or you can set this to specify the (fasta) database file
#' to be searched instead of using standard database.
#' A blast database will be created if it does not exist.
#' @param fill a numeric value for filling default NA value when no
#' blast result.
#' @param threads number of threads to run.
#' @param tmp_dir path for storing temp files.
#' @param clean_tmp if `TRUE`, remove temp directory.
#'
#' @return Data table of IEDB scores.
#' - peptide - input peptide
#' - iedb_score - IEDB score
#' - annotation - IEDB annotation info
#' @import data.table
#' @importFrom data.table :=
#' @export
#' @examples
#' \donttest{
#' calc_iedb_score("AAAAAAAAA")
#' calc_iedb_score("MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDP")
#' calc_iedb_score("MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDP", db = "mouse")
#' calc_iedb_score(c("MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDP", "MTEYKLVVVG"))
#' }
#' @seealso [set_blast_path] for setting blast path, [install_database]
#' for installing database.
calc_iedb_score <- function(pep,
                            db = "human",
                            fill = NA_real_,
                            threads = parallel::detectCores(),
                            tmp_dir = file.path(tempdir(), "neopeptides"),
                            clean_tmp = TRUE) {
  stopifnot(length(db) == 1)
  old <- data.table::getDTthreads()
  data.table::setDTthreads(threads = threads)

  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE)
  }

  if (clean_tmp) {
    on.exit({
      message("=> Removing temporary files...")
      try(
        unlink(tmp_dir, recursive = TRUE)
      )
      data.table::setDTthreads(old)
    })
  }

  # generate fastas to query
  names(pep) <- 1:length(pep)
  sdt <- data.table::data.table(
    pep = pep,
    id = names(pep)
  )

  message("=> Running blastp for homology to IEDB antigens..")
  tmp_fasta <- file.path(tmp_dir, "iedb_score.fa")
  make_fasta_file(pep, tmp_fasta)
  tmp_iedb_out <- file.path(tmp_dir, "blastp_iedbout.csv")
  db <- make_blastp_db(db, data_type = "IEDB")
  # Run blast
  cmd_blastp(tmp_fasta, db, tmp_iedb_out, threads = threads)

  blastdt <- list.files(
    path = tmp_dir,
    pattern = "blastp_iedbout\\.csv",
    full.names = TRUE
  )

  if (length(blastdt) == 0) {
    message("=> No blast output against database returned!")
    return(data.table::data.table(peptide = pep, iedb_score = fill, annotation = NA_character_))
  }

  if (all(file.info(blastdt)$size == 0)) {
    message("=> No database matches found by blast!")
    return(data.table::data.table(peptide = pep, iedb_score = fill, annotation = NA_character_))
  }

  blastdt <- read_blast_result(blastdt)

  if (nrow(blastdt) == 0) {
    message(paste("=> No IEDB matches found with cannonical AAs, can't compute IEDB score...."))
    return(data.table::data.table(peptide = pep, iedb_score = fill, annotation = NA_character_))
  }

  message("=> Summing IEDB local alignments...")
  blastdt[, SW := SW_align(nmer, WT.peptide)]
  message("=> Done.")

  blastdt[, score := SW %>% modeleR(), by = "id"]

  # get full IEDB ref here
  fa <- Biostrings::readAAStringSet(db)
  f <- fa %>% as.character()
  names(f) <- names(fa)

  blastdt[, anno := lapply(anno, function(i) {
    mv <- f[which(stringr::str_detect(pattern = stringr::fixed(i), names(f)))]
    mv <- mv[which(stringr::str_detect(pattern = stringr::fixed(WT.peptide), mv))]
    return(paste(names(mv), WT.peptide, collapse = "|"))
  }), by = 1:nrow(blastdt)]

  anndt <- blastdt[, .SD %>% unique(), .SDcols = c(
    "id",
    "nmer",
    "anno",
    "score",
    "SW"
  )]

  blastdt <- blastdt[, .SD %>% unique(), .SDcols = c("id", "score")]
  sdt[, id := as.character(id)]
  blastdt[, id := as.character(id)]

  sdt <- merge(sdt, blastdt, by = "id", all.x = TRUE)

  anndt <- anndt[, msw := max(SW), by = "id"] %>%
    .[SW == msw] %>%
    .[, .SD %>% unique(), .SDcols = c("id", "anno")]

  # merge equally good annos into one
  anndt[, anno := paste(anno %>% unique(), collapse = "|"), by = "id"]
  anndt <- anndt %>% unique()
  anndt[, id := as.character(id)]
  sdt <- merge(sdt, anndt, by = "id", all.x = TRUE)
  sdt <- sdt[order(as.integer(id)), .SD %>% unique(),
    .SDcols = c("id", "pep", "score", "anno")
  ][
    , c("pep", "score", "anno"),
    with = FALSE
  ]

  colnames(sdt) <- c("peptide", "iedb_score", "annotation")
  return(sdt)
}


utils::globalVariables(
  c(
    ".", "anno", "SW", "WT.peptide", "score",
    "msw", "nmer", "id"
  )
)
