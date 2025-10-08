#' @title fetch_alphafold
#'
#' @description Downloads the Alphafold2 protein structure prediction in PDB format based on the provided protein ID.
#'
#' @param protein_id A character vector specifying the protein ID for which the prediction is to be fetched.
#'
#' @return If successful, the function returns the path to the downloaded Alphafold2 PDB file. If the file already
#' exists in the specified directory, a message is displayed, and the existing file path is returned. If the download
#' fails, an error message is displayed, and the function returns NULL.
#'
#' @examples
#' # Download the Alphafold2 PDB for a protein with ID "PDB_ID"
#' fetch_alphafold("P57102")
#'
#' @export

fetch_alphafold <- function(protein_id) {
  # Check that protein_id is a character vector
  stopifnot(is.character(protein_id))

  # Set the base URL for downloading Alphafold2 protein structure predictions
  base_url <- "https://alphafold.ebi.ac.uk/files/"

  # Set the version suffix
  version_suffix <- "-F1-model_v6.pdb"

  # Define output filename
  filename <- paste0("AF-", protein_id, version_suffix)

  # Add path to the output file
  output_file <- file.path(model_folder, filename)  # Set destination folder (change if needed)

  # Check if the file already exists in the specified directory
  if (file.exists(output_file)) {
    message("Alphafold2 PDB for ", protein_id, " already exists.")
    return(output_file)
  }

  # Define download URL
  download_url <- paste0(base_url, filename)

  # Try to download the Alphafold2 protein structure prediction using the specified URL and save it to the output file
  res <- tryCatch(
    {
      if (httr::http_error(GET(download_url))) {
        return(NA)
        } else {
        download.file(download_url, destfile = output_file, method = "curl", quiet = TRUE)  
        return(output_file)
        }
    },
    # If the download fails, display an error message
    error = function(cond) {
      message(paste("Error while downloading: ", protein_id))
      message(cond)
      return(NULL)
    }
  )

  return(res)
}

