source("code/ss_helper.R")
source("code/config.R")
# Vector of package names
pkg_names <- c("UniProt.ws","httr", "rMSA", "Biostrings", "rvest")

# Check if packages are already installed
pkg_installed <- pkg_names %in% installed.packages()

# Vector of packages that are not installed
pkg_to_install <- pkg_names[!pkg_installed]

# Install packages that are not already installed
if (length(pkg_to_install) > 0) {
  install.packages(pkg_to_install, 
    dependencies = TRUE,  
    repos = c("https://cran.r-project.org", "https://bioconductor.org/packages/3.12/bioc", "https://mhahsler.r-universe.dev"))
}
for (pkg in pkg_names){library(pkg, character.only = TRUE)}

# download CDS files for species from NCBI FTP
cds_links <- ftp_search(species)
lapply(cds_links, ftp_download)

# convert cDNA to protein
# make protein fasta into blasteable database
cdna_files <- list.files(cds_folder, pattern = "all.fa$", full.names = TRUE, recursive = TRUE)
lapply(cdna_files, make_protein_db)

