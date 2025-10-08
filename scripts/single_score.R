#!/usr/bin/env Rscript
# Ensure R uses the conda environment library
.libPaths(file.path(Sys.getenv("CONDA_PREFIX"), "Lib", "R", "library"))

# load required packages
library(epictope)
rm(list = ls())

# helper functions
setup_files()
check_config()
args <- commandArgs(trailingOnly = TRUE); query <- args[1]
# download uniprot information 
uniprot_fields <- c("accession", "id", "gene_names", "xref_alphafolddb", "sequence", "organism_name", "organism_id")
uniprot_data <- query_uniProt(query = query, fields = uniprot_fields)

####
# download associated alphafold2 pdb
if(is.na(uniprot_data$AlphaFoldDB)){
  message("WARNING. NO CROSSREFERENCE ALPHAFOLD ENTRY FOUND; ATTEMPTING DIRECT LOOKUP. RESULTS MAY BE LOW CONFIDENCE.")
  uniprot_data$AlphaFoldDB <- query}
alphafold_file <- fetch_alphafold(gsub(";", "", uniprot_data$AlphaFoldDB))
if(is.na(alphafold_file)){stop("No alphafold file found for ", uniprot_id)}
# calculate dssp on alphafold pdb file
dssp_res <- dssp_command(alphafold_file)
# parse and read in dssp
dssp_df <- parse_dssp(dssp_res)

# retrieve iupred/anchor2 disordered binding regions
iupred_df <- iupredAnchor(query)

# blast query aa sequence
seq <- Biostrings::AAStringSet(uniprot_data$Sequence, start=NA, end=NA, width=NA, use.names=TRUE)

# list of amino acid files to blast against
aa_files <-  list.files(cds_folder, pattern = paste0(species, ".*\\.all.fa$", collapse = "|"), ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
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
msa_res <-  muscle(blast_stringset)
# write msa to file
# convert to XStringSet
msa_res <- Biostrings::AAStringSet(msa_res)
Biostrings::writeXStringSet(msa_res, file = paste0(outputFolder, "/", query, "_msa.fasta"))
# shannon entropy calculation
shannon_df <- shannon_reshape(msa_res, query)

# join tagging features in dataframe
features_df <- Reduce(function(x, y) merge(x, y, all=TRUE), list(shannon_df, dssp_df, iupred_df), accumulate=FALSE)
colnames(features_df)
# normalize features and calculate tagging score
norm_feats_df <- calculate_scores(features_df)

# write to file.
res_df <- merge(norm_feats_df, features_df)
write.csv(apply(res_df,2,as.character), file = paste0(outputFolder, "/", query, "_score.csv"), row.names = FALSE)
