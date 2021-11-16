#' Calculate Dissimilarity Value to Reference Proteome for Peptides
#'
#' @inheritParams calc_iedb_score
#' @param k_val numeric. Steepness of sigmoidal curve at k.
#' Default 4.86936, the value used in the analysis of Van Allen,
#' Snyder, Rizvi, Riaz, and Hellmann datasets. See reference.
#' @param a_val numeric. Optionally can be "mean" to use mean alignment
#' for nmers passed. Horizontal displacement of partition function.
#' Default is 32, based on max_SW of 75 million 8-15mers from the
#' five clinical datasets against human, if using max_SW, use 52.
#' This value may not be meaningful for murine alignment so use with care.
#' See reference.
#'
#' @return Data table of dissimilarity values (to the non-mutated proteome).
#' - peptide - input peptide
#' - dissimilarity - dissimilarity value
#' @references Richman LP, Vonderheide RH, and Rech AJ.
#' "Neoantigen dissimilarity to the self-proteome predicts immunogenicity
#' and response to immune checkpoint blockade." Cell Systems 9, 375-382.E4, (2019).
#'
#' @import data.table
#' @importFrom data.table :=
#' @export
#' @examples
#' \dontrun{
#' calc_dissimilarity("AAAAAAAAA")
#' calc_dissimilarity("MRLVDRRWA")
#' calc_dissimilarity("VRLVDRRWA")
#' calc_dissimilarity("MTEYKLVVVGAGDVGKSALTI")
#' calc_dissimilarity("MTEYKLVVVGAGDVGKSALTI", db = "mouse")
#' calc_dissimilarity(c("MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDP", "MTEYKLVVVG"))
#' }
#' @seealso [set_blast_path] for setting blast path, [install_database]
#' for installing database.
calc_dissimilarity <- function(pep,
                               db = "human",
                               k_val = 4.86936,
                               a_val = 32,
                               fill = NA_real_,
                               use_blastp_short = TRUE,
                               threads = parallel::detectCores(),
                               tmp_dir = file.path(tempdir(), "neopeptides"),
                               clean_tmp = TRUE) {
  stopifnot(length(db) == 1, is.character(pep))
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

  message("=> Running blastp for homology to self antigens..")
  tmp_fasta <- file.path(tmp_dir, "dissimilarity.fa")
  make_fasta_file(pep, tmp_fasta)
  tmp_blastp_out <- file.path(tmp_dir, "blastp_self_out.csv")
  db <- make_blastp_db(db, data_type = "Proteome")
  # Run blast
  cmd_blastp(tmp_fasta, db, tmp_blastp_out, short_task = use_blastp_short, threads = threads)

  blastdt <- list.files(
    path = tmp_dir,
    pattern = "blastp_self_out\\.csv",
    full.names = TRUE
  )

  if (length(blastdt) == 0) {
    message("=> No blast output against self-proteome returned!")
    return(data.table::data.table(peptide = pep, dissimilarity = fill))
  }

  if (all(file.info(blastdt)$size == 0)) {
    message("=> No self-proteome matches found by blast!")
    return(data.table::data.table(peptide = pep, dissimilarity = fill))
  }

  blastdt <- read_blast_result(blastdt)

  if (nrow(blastdt) == 0) {
    message(paste("=> No self-proteome matches found with cannonical AAs, can't compute dissimilarity...."))
    return(data.table::data.table(peptide = pep, dissimilarity = fill))
  }

  # May memory intense
  message("=> Summing local alignments...")
  blastdt[, SW := SW_align(nmer, WT.peptide)]
  message("=> Done.")

  blastdt[, score := SW %>% modeleR(a = a_val, k = k_val, dislike = TRUE),
    by = "id"
  ]

  blastdt <- blastdt[, .SD %>% unique(), .SDcols = c("id", "score")]
  sdt[, id := as.character(id)]
  blastdt[, id := as.character(id)]

  sdt <- merge(sdt, blastdt, by = "id", all.x = TRUE)
  sdt <- sdt[order(as.integer(id)),
    .SD %>% unique(),
    .SDcols = c("id", "pep", "score")
  ][
    , c("pep", "score"),
    with = FALSE
  ]

  if (!is.na(fill)) {
    sdt$score[is.na(sdt$score)] <- fill
  }

  colnames(sdt) <- c("peptide", "dissimilarity")
  return(sdt)
}
