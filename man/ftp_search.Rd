% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ftp_search.R
\name{ftp_search}
\alias{ftp_search}
\title{ftp_search}
\usage{
ftp_search(species, ftp = "https://ftp.ensembl.org/pub/release-108/fasta/")
}
\arguments{
\item{species}{A character vector of species names to search for}

\item{ftp}{A character vector of specifying the FTP site to search, defaults}
}
\value{
A named character vector of links to the species on the NCBI FTP site
}
\description{
Searches the NCBI FTP directory for a a specified species folder
}
\examples{
ftp_link <- ftp_search("Mus_musculus", ftp = "https://ftp.ensembl.org/pub/release-108/fasta/")
}
