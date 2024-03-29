% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blast.R
\name{blast}
\alias{blast}
\title{blast}
\usage{
blast(db = NULL, type = "blastp")
}
\arguments{
\item{db}{A character string specifying the path to the BLAST database.}

\item{type}{A character string specifying the type of BLAST (e.g., "blastp", "blastn", "blastx", etc.).}
}
\value{
The function returns an object of class "BLAST" containing information about the BLAST database.
}
\description{
Create a BLAST object for working with the BLAST database.
}
