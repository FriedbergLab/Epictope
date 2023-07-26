#' @title muscle
#'
#' @description Perform multiple sequence alignment using the MUSCLE program.
#'
#' @param x A set of biological sequences (RNAStringSet, DNAStringSet, or AAStringSet).
#'
#' @return The function returns the multiple sequence alignment result as an RNAStringSet, DNAStringSet, or AAStringSet, based on the input sequence type.
#'
#' @examples
#' # Create a set of DNA sequences
#' #dna_sequences <- DNAStringSet(c("ATCGATCG", "AT-AT-AT", "GCG-GCG"))
#' # Perform multiple sequence alignment
#' #alignment <- muscle(dna_sequences)
#'
#' @importFrom Biostrings readRNAMultipleAlignment readDNAMultipleAlignment readAAMultipleAlignment writeXStringSet
#' @export


#   https://github.com/mhahsler/rMSA
#   Copyright (C) 2012 Michael Hahsler and Anurag Nagar
muscle <- function(x) {
  # set up temp directory and files
  wd <- temp_folder
  dir <- getwd()
  temp_file <- basename(tempfile(tmpdir = wd))
  on.exit({
    file.remove(Sys.glob(paste(temp_file, ".*", sep = "")))
    setwd(dir)
  })
  infile <- file.path(wd, paste(temp_file, ".in", sep = ""))
  outfile <- file.path(wd, paste(temp_file, ".aln", sep = ""))

  # check sequence type
  reader <- if (is(x, "RNAStringSet")) {
    readRNAMultipleAlignment
  } else if (is(x, "DNAStringSet")) {
    readDNAMultipleAlignment
  } else if (is(x, "AAStringSet")) {
    readAAMultipleAlignment
  } else {
    stop("Unknown sequence type!")
  }

  # Write sequences to a temporary FASTA input file
  writeXStringSet(x, infile, append = FALSE, format = "fasta")

  # call muscle (MUSCLE executable needs to be installed and in the path!)
  system(paste(.findExecutable("muscle"), "-align", infile, "-output", outfile))

  # Read and process the output of MUSCLE
  rows <- scan(outfile, what = "", sep = "\n", strip.white = FALSE,
               quiet = TRUE, blank.lines.skip = FALSE)
  if (rows[[3L]] != "") {
    rows <- c(rows[1:2], "", rows[3:length(rows)])
  }
  cat(rows, file = outfile, sep = "\n")

  # Read the processed alignment using the appropriate reader
  result_alignment <- reader(outfile, format = "fasta")

  # Return the result
  return(result_alignment)
}
