# Script Title: Epitope Scoring function
# Description: Calculates a epitope tag insertion score for a uniprot id
# Author: HC
# Date: 04/03/23

# load required packages
pkg_names <- c("UniProt.ws","httr", "rMSA", "Biostrings")
for (pkg in pkg_names){library(pkg, character.only = TRUE)}
rm(list = ls())

# helper functions
source("code/ss_helper.R")
# config file
source("code/config.R")

# download uniprot information 
uniprot_fields <- c("accession", "id", "gene_names", "xref_alphafolddb", "sequence", "organism_name", "organism_id")
uniprot_data <- query_uniProt(query = query, fields = uniprot_fields)

####
# download associated alphafold2 pdb
alphafold_file <- fetch_alphafold(gsub(";", "", uniprot_data$AlphaFoldDB))
# calculate dssp on alphafold pdb file
dssp_res <- dssp_command(alphafold_file)
# parse and read in dssp
dssp_df <- parse_dssp(dssp_res)

# retrieve iupred/anchor2 disordered binding regions
iupred_df <- idpr::iupredAnchor(query, plotResults = FALSE)
colnames(iupred_df) <- tolower(colnames(iupred_df))

# blast query aa sequence
seq <- Biostrings::AAStringSet(uniprot_data$Sequence, start=NA, end=NA, width=NA, use.names=TRUE)

# list of amino acid files to blast against
aa_files <-  list.files(cds_folder, pattern = "all.amino.fa$", full.names = TRUE, recursive = TRUE)
names(aa_files) <- aa_files

# blast
blast_results <- lapply(aa_files, function(.x){protein_blast(seq, .x)})

# take the highest blast match according to E score 
find_best_match <- function(.x){head(.x[base::order(.x$E),], 1)}
blast_best_match <- lapply(blast_results, find_best_match)
blast_seqs <- lapply(blast_best_match, fetch_sequences)

# fetch the AA sequence for blast best match
blast_seqs[[query]] <- seq
blast_stringset  <- Biostrings::AAStringSet(unlist(lapply(blast_seqs, function(.x){.x[[1]]})))

# multiple sequence alignment
msa_res <-  rMSA::muscle(blast_stringset)

# shannon entropy calculation
shannon_df <- shannon_reshape(msa_res, query)

# join tagging features in dataframe
features_df <- Reduce(function(x, y) merge(x, y, all=TRUE), list(shannon_df, dssp_df, iupred_df), accumulate=FALSE)
colnames(features_df)
# normalize features and calculate tagging score
norm_feats_df <- calculate_scores(features_df)

# write to file.
res_df <- merge(norm_feats_df, features_df)
write.csv(apply(res_df,2,as.character), file = paste0(outputFolder, "/", query, "_score.csv"))