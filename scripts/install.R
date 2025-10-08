#!/usr/bin/env Rscript
# Ensure R uses the conda environment library
.libPaths(file.path(Sys.getenv("CONDA_PREFIX"), "Lib", "R", "library"))

library(epictope)
rm(list = ls())
# load config and set up file structure.
setup_files()
check_config()
# download CDS files for species from NCBI FTP
cds_links <- ftp_search(species = species)
lapply(cds_links, ftp_download)

# convert cDNA to protein
# make protein fasta into blasteable database
cdna_files <- list.files(cds_folder, pattern = "\\.all.fa.gz$", full.names = TRUE, recursive = TRUE)
lapply(cdna_files, make_protein_db)


