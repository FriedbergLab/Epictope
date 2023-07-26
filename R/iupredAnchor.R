#' @title iupredAnchor
#'
#' @description Predicts intrinsic disorder and anchor residues in a protein using the IUPred2A web service.
#'
#' @param uniprotAccession A character string specifying the UniProt accession number of the protein.
#' @param proteinName A character string specifying the name of the protein (optional).
#'
#' @return The function returns a data frame containing the predicted intrinsic disorder and anchor residues at each position in the protein sequence. If \code{plotResults = TRUE}, it returns a ggplot object showing the prediction plot. If \code{plotResults = FALSE}, it returns the data frame only.
#'
#' @importFrom jsonlite fromJSON
#' @export


#   http://www.bioconductor.org/packages/release/bioc/html/idpr.html
#   Authors: William M. McFadden, Judith L. Yanowitz, Michael Buszczak
iupredAnchor <- function(uniprotAccession, proteinName = NA) {
  # Connecting to IUPred2A REST API
  iupredURL <- paste("https://iupred2a.elte.hu/iupred2a/anchor/",
                     uniprotAccession, ".json", sep = "")
  iupredJson <- jsonlite::fromJSON(iupredURL)

  # Reformatting data to be consistent in formatting across idpr
  iupredPrediction <- iupredJson$iupred2
  anchorPrediction <- iupredJson$anchor2
  iupredSequence <- unlist(strsplit(iupredJson$sequence, ""))
  iupredSequence <- unlist(iupredSequence)
  seqLength <- length(iupredSequence)
  iupredDF <- data.frame(Position = seq_len(seqLength),
                         AA = iupredSequence,
                         IUPred2 = iupredPrediction,
                         ANCHOR2 = anchorPrediction)
  colnames(iupredDF) <- tolower(colnames(iupredDF))
  return(iupredDF)
}
