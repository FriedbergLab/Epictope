#' @title calculate_scores
#'
#' @description Calculates various scores for each position in the feature data frame.
#'
#' @param feature_df A data frame containing features for each position.
#'
#' @return The function returns a data frame with calculated scores for each position.
#'
#' @export

calculate_scores <- function(feature_df) {
  # Normalize entropy by dividing by max
  normalized_entropy <- feature_df$shannon / 4.321928

  # Convert secondary structure symbols to numeric value
  ss_score <- as.numeric(sapply(feature_df$ss, ss_convert))

  # Calculate relative solvent accessibility
  rsa <- feature_df$sasa / max_sasa[feature_df$aa]

  # Invert the anchor2 score
  inv_anchor2 <- 1 - feature_df$anchor2

  # Calculate the final score using a weighted sum of features
  tag_score <- (h_weight * normalized_entropy) + (ss_weight * ss_score) + (rsa_weight * rsa) + (br_weight * inv_anchor2)

  # Bind scores to the dataframe.
  res <- data.frame(
    position = 1:nrow(feature_df),
    normalized_entropy = normalized_entropy,
    ss_score = ss_score,
    rsa = rsa,
    inv_anchor2 = inv_anchor2,
    tag_score = tag_score
  )

  # Determine the minimum feature value and corresponding feature for each position.
  res$min <- apply(res, 1, FUN = function(.x) { return(min(.x)) })
  res$min_feature <- apply(res[c("normalized_entropy", "ss_score", "rsa", "inv_anchor2")], 1, FUN = function(.x) { colnames(res)[which(.x == min(.x))] })

  # Return the result.
  return(res)
}
