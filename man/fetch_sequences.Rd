% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch_sequence.R
\name{fetch_sequences}
\alias{fetch_sequences}
\title{fetch_sequences}
\usage{
fetch_sequences(blast_results)
}
\arguments{
\item{blast_results}{A list containing BLAST results, including the SubjectID and filename of the BLAST database.}
}
\value{
The function returns the sequences that match the SubjectID from the BLAST results.
}
\description{
Reads a FASTA file and fetches the sequences that match the SubjectID from the BLAST results.
}
