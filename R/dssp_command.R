#' @title dssp_command
#'
#' @description Runs the DSSP command with mkdssp to process a PDB file and generates a DSSP file.
#'
#' @param pdb_file A character vector specifying the path to the input PDB file.
#'
#' @return The function returns the name of the output DSSP file.
#'
#' @examples
#' # Generate the DSSP file from the input PDB file "protein.pdb"
#' #dssp_command("protein.pdb")
#'
#' @export

dssp_command <- function(pdb_file) {
  # Check that .x is a character vector
  stopifnot(is.character(pdb_file))

  # Define output filename by replacing the extension with ".dssp"
  output_file <- gsub("\\.pdb.gz|\\.pdb", ".dssp", pdb_file)

  # Run the DSSP command with mkdssp, return error message if an error occurs.
  system(command = paste0("mkdssp -i ", pdb_file, " -o ", output_file), intern = TRUE)

  # Return the name of the output file.
  return(output_file)
}
