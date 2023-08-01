#' @title make_protein_db
#'
#' @description Decompresses the input gzipped file and creates a protein BLASTable database using makeblastdb.
#'
#' @param gz_file A character vector specifying the input gzipped file to be used for creating the database.
#'
#' @examples
#' # Create a protein BLASTable database from the compressed FASTA file
#' #make_protein_db("Bos_taurus.ARS-UCD1.2.pep.all.fa.gz")
#'
#' @export

#' @title make_protein_db
#'
#' @description Decompresses the input gzipped file and creates a protein BLASTable database using makeblastdb.
#'
#' @param gz_file A character vector specifying the input gzipped file to be used for creating the database.
#'
#' @examples
#' # Create a protein BLASTable database from the compressed FASTA file
#' #make_protein_db("Bos_taurus.ARS-UCD1.2.pep.all.fa.gz")
#'
#' @importFrom R.utils gunzip
#' @export
#' 
make_protein_db <- function(gz_file) {
  # Remove the ".gz" extension from the input file name
  file <- gsub(".gz", "", gz_file)

  # Check if the unzipped file already exists
  if (!file.exists(file)) {
    # Decompress the input file
    R.utils::gunzip(gz_file, remove = FALSE)
  }

  # Create the protein BLASTable database
  system(command = paste("makeblastdb -in ", file, " -dbtype prot"))
}