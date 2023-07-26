#' @title shannon_entropy
#'
#' @description Calculates the Shannon entropy of a sequence.
#'
#' @param seq A vector representing the sequence for which the entropy is to be calculated.
#' @param nogap A logical value indicating whether to exclude gaps (default is FALSE).
#'
#' @return The function returns the Shannon entropy of the input sequence.
#'
#' @examples
#' # Calculate the Shannon entropy of a sequence
#' entropy <- shannon_entropy(c("A", "A", "T", "C", "G"))
#'
#' @export

shannon_entropy <- function(seq, nogap = FALSE) {
  # Ensure that seq is a vector, not a list.
  seq <- unlist(seq)

  # Take the unique bases in the sequence.
  unique_base <- unique(seq)

  if (nogap) {
    # Exclude gaps if specified.
    unique_base <- unique_base[unique_base != "-"]
    adj_seq <- seq[seq != "-"]
  } else {
    adj_seq <- seq
  }

  # Get the length of the adjusted sequence.
  M <- length(adj_seq)

  # Initialize a vector to store the entropy for each base.
  entropy_list <- numeric(length(unique_base))

  for (i in seq_along(unique_base)) {
    # Take the current base.
    base <- unique_base[i]

    # Count the number of occurrences.
    n_i <- sum(adj_seq == base)

    # Calculate the probability.
    P_i <- n_i / M

    # Calculate the entropy.
    entropy_i <- P_i * log2(P_i)

    # Add entropy to the list.
    entropy_list[i] <- entropy_i
  }

  # Calculate the total entropy by summing the entropies for each base.
  res <- -sum(entropy_list)

  # Return the result.
  return(res)
}
