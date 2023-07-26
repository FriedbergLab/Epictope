#' @title query_uniProt
#'
#' @description Queries the UniProt API and retrieves information based on the provided query and fields.
#'
#' @param query A character vector specifying the query terms to search for in UniProt.
#' @param fields A character vector specifying the fields to include in the query results (default is c("accession", "id")).
#' @param collapse A character vector specifying the separator for combining multiple query terms (default is " OR").
#'
#' @importFrom httr GET content
#'
#' @return A data frame containing the results of the UniProt query.
#'
#' @examples
#' # Query UniProt for information about proteins with accessions "P04114" and "P15056"
#' query_uniProt(query = c("P04114", "P15056"))
#'
#' @export

query_uniProt <- function(query = character(0L), fields = c("accession", "id"), collapse = " OR") {
  # Check that query and fields are both character vectors
  stopifnot(is.character(query), is.character(fields))

  # If the length of query is zero, throw an error message
  if (!length(query)) {
    stop("<internal> 'query' must be populated with valid inputs.")
  }

  # Send an HTTP GET request to the UniProt API, passing the query and field parameters
  # and specifying the desired response format
  resp <- httr::GET(
    url = "https://rest.uniprot.org/uniprotkb/search",
    query = list(query = paste(query, collapse = collapse), fields = paste(fields, collapse = ","), format = "tsv")
  )

  # Read the response text into a data frame using the read.delim function and return
  read.delim(text = httr::content(resp, encoding = "UTF-8"))
}
