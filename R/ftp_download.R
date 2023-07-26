#' @title ftp_download
#'
#' @description Downloads the Ensembl CDS file for a specified species from the NCBI FTP site.
#'
#' @param link A character vector specifying the species link to download.
#' @param seq_type A character vector specifying the sequence type (default is "pep").
#'
#' @return The downloaded file will be saved to the specified destination folder.
#'
#' @importFrom rvest read_html html_attr html_nodes
#'
#' @examples
#' # Download the CDS file for Mus musculus (protein sequences)
#' #ftp_download("Mus_musculus/")
#' @export

ftp_download <- function(link, seq_type = "pep") {
  # Construct the URL for the CDS directory of the given link
  link_url <- paste0("https://ftp.ensembl.org/pub/release-108/fasta/", link, seq_type, "/")
  # Read webpage
  webpage <- rvest::read_html(link_url)
  # Extract links
  cds_link <- rvest::html_attr(rvest::html_nodes(webpage, "a"), "href")
  # Filter to ".fa.gz" extension
  filename <- cds_link[grepl("\\.all.fa\\.gz", cds_link)]
  # Download
  tryCatch(
    {
      utils::download.file(
        url = paste0(link_url, filename),
        destfile = paste0(getwd(), "/", filename),  # Set destination folder (change if needed)
        method = "wget",
        extra = "-nc"
      )
    },
    # If the download fails, search for the file in the specified directory and return if found
    error = function(cond) {
      message(cond$message)
      return(NULL)
    }
  )
}
