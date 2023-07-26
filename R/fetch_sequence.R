#' @title fetch_sequences
#'
#' @description Reads a FASTA file and fetches the sequences that match the SubjectID from the BLAST results.
#'
#' @param blast_results A list containing BLAST results, including the SubjectID and filename of the BLAST database.
#'
#' @importFrom Biostrings readAAStringSet
#'
#' @return The function returns the sequences that match the SubjectID from the BLAST results.
#'
#' @export

fetch_sequences <- function(blast_results) {
  # Read in the FASTA file
  cur_proteins <- Biostrings::readAAStringSet(blast_results$filename, format = "fasta")

  # Subset sequences by query name (SubjectID)
  cur_seqs <- cur_proteins[grepl(blast_results$SubjectID, names(cur_proteins))]

  # Return the matched sequences
  return(cur_seqs)
}
