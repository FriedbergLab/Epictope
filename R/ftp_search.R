#' @title ftp_search
#'
#' @description Searches the NCBI FTP directory for a a specified species folder
#'
#' @param species  A character vector of species names to search for
#' @param ftp A character vector of specifying the FTP site to search, defaults
#        to "https://ftp.ensembl.org/pub/release-108/fasta/"
#'
#' @return A named character vector of links to the species on the NCBI FTP site
#' @examples
#' ftp_link <- ftp_search("Mus_musculus", ftp = "https://ftp.ensembl.org/pub/release-108/fasta/")
#' @export
#' @importFrom rvest read_html html_attr html_text html_nodes
#'
####

ftp_search <- function(species, ftp = "https://ftp.ensembl.org/pub/release-108/fasta/") {
  # List of NCBI species
  webpage <- rvest::read_html(ftp)
  # Extract links
  links <- rvest::html_attr(rvest::html_nodes(webpage, "a"), "href")
  link_names <- rvest::html_text(rvest::html_nodes(webpage, "a"))
  names(links) <- link_names
  # Search for species
  species_regex <- paste0(species, "/$")
  species_query <- paste(paste0(species, "/$"), collapse = "|")
  species_links <- links[base::grepl(species_query, links, perl = TRUE)]
  not_found <- setdiff(species, gsub("/", "", species_links))
  if (length(not_found) > 0) {
    message(paste("The following species were not found:", paste(not_found, collapse = ", ")))
  }

  # Return links
  return(species_links)
}
