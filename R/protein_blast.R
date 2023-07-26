#' @title protein_blast
#'
#' @description Performs protein BLAST using the specified BLAST database file and predicts matches for input sequences.
#'
#' @param seqs A character vector of input sequences to be used in the BLAST search.
#' @param blast_file A character vector specifying the path to the BLAST database file.
#' @param blast_type A character vector specifying the type of BLAST (default is "blastp").
#'
#' @return The function returns the BLAST results as a list, including the matches and information about the BLAST database used.
#'
#' @examples
#' # Perform protein BLAST using the "database.fasta" database file on input sequences
#' #blast_results <- protein_blast(seqs, blast_file = "database.fasta", blast_type = "blastp")
#'
#' @export

protein_blast <- function(seqs, blast_file, blast_type = "blastp") {
  # Check that blast_file is a character vector
  stopifnot(is.character(blast_file))

  # Create a BLAST database object
  db <- blast(db = blast_file, type = blast_type)

  # Run BLAST and predict matches for the input sequences
  blast_res <- predict(db, seqs)

  # Attach the filename of the BLAST database to the results
  blast_res$filename <- blast_file

  # Return the BLAST results
  return(blast_res)
}
