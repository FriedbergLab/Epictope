#' @title shannon_reshape
#'
#' @description Calculates Shannon entropy for each position of the query sequence in the multiple sequence alignment (MSA).
#'
#' @param msa A multiple sequence alignment object (e.g., from the `msa` package).
#' @param query The index or name of the query sequence in the MSA.
#'
#' @return The function returns a data frame containing the positions, amino acids, and their corresponding Shannon entropy values for the query sequence.
#'
#' @importFrom Biostrings readAAStringSet
#' @export

shannon_reshape <- function(msa, query) {
  # Extract MSA as a matrix
  mat <- as.matrix(msa)

  # Determine the row for the query sequence
  chars <- mat[query,]

  # Calculate conservation on query positions
  ind <- which(chars != "-")

  # Calculate Shannon entropy for each column of the MSA
  score <- apply(mat, 2, shannon_entropy, nogap = FALSE)

  # Return Shannon entropy for only query positions
  temp <- score[ind]

  # Create a data frame to store the results
  res <- data.frame(position = 1:length(temp), aa = chars[ind], shannon = temp)

  # Return the result
  return(res)
}
