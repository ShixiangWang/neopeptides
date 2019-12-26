cmd_blastp <- function(query_file, db, out_file, threads = parallel::detectCores()) {
  # Run blastp-short
  # https://www.ncbi.nlm.nih.gov/books/NBK279684/
  # flags here taken from Lukza et al.:
  # -task blastp-short optimized blast for <30 AA, uses larger word sizes
  # -matrix use BLOSUM62 sub substitution matrix
  # -evalue expect value for saving hits
  # -gapopen, -gapextend, numeric cost to a gapped alignment and
  # -outfmt, output a csv with colums, seqids for query and database seuqnence, start and end of sequence match,
  # length of overlap, number of mismatches, percent identical, expected value, bitscore
  cmds <- paste(
    "blastp",
    "-query",
    query_file,
    "-db",
    db,
    "-evalue 100000000",
    "-matrix BLOSUM62",
    "-gapopen 11",
    "-gapextend 1",
    "-out",
    out_file,
    "-num_threads",
    threads,
    "-outfmt '10 qseqid sseqid qseq qstart qend sseq sstart send length mismatch pident evalue bitscore'"
  )
  system(cmds)
}

cmd_blastdb <- function(x, dbtype="prot") {
  # Create a database name as x
  stopifnot(has_program("makeblastdb"))
  # makeblastdb -in Mu_iedb.fasta -parse_seqids -hash_index -dbtype prot
  # -out Mu_iedb.fasta can be ignored
  cmds <- paste(
    "makeblastdb",
    "-in",
    x,
    "-parse_seqids",
    "-hash_index",
    "-dbtype",
    dbtype
  )
  system(cmds)
}
