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

make_protein_db <- function(gz_file) {
  # Decompress the input file
  system(paste0("gunzip -d -k ", gz_file))

  # Remove the ".gz" extension from the file name
  file <- tools::file_path_sans_ext(gz_file)

  # Create the protein BLASTable database
  system(command = paste("makeblastdb -in ", file, " -dbtype prot"))
}
