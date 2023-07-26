library(epictope)
# download CDS files for species from NCBI FTP
cds_links <- ftp_search(species)
lapply(cds_links, ftp_download)

# convert cDNA to protein
# make protein fasta into blasteable database
cdna_files <- list.files(cds_folder, pattern = "all.fa$", full.names = TRUE, recursive = TRUE)
lapply(cdna_files, make_protein_db)

