#' @title ss_convert
#'
#' @description Converts a Secondary Structure Symbol to its corresponding numeric value using the SS_key.
#'
#' @param .x A character vector representing the Secondary Structure Symbol to be converted.
#'
#' @return The function returns the numeric value corresponding to the Secondary Structure Symbol as defined in the SS_key.
#'
#' @examples
#' # Convert Secondary Structure Symbols to numeric values
#' ss_convert("H")  # Output: 1
#' ss_convert("E")  # Output: 2
#'
#' @export

ss_convert <- function(.x) {
  # Check if .x is a recognized Secondary Structure Symbol
  if (!grepl("[GHIECTBSP]", .x)) {
    # If not, print an error message and return NULL
    message("Character is not a recognized Secondary Structure Symbol.")
    return(NULL)
  } else {
    # Return the numeric value corresponding to the Secondary Structure Symbol as defined in the SS_key
    return(ss_key[.x])
  }
}
