library(seqinr)
library(tidyverse)
source("code/seqlogo_helper.R")

# Read in MSA Clustal files
clustal_folder <- "outputs/fastas/"
clustal_files <- list.files(clustal_folder, pattern = ".clw")
clustals = lapply(clustal_files, read_clustal)

clustal_dfs <- lapply(clustals, clustal_df)
names(clustal_dfs) = tolower(stringr::str_remove(clustal_files, ".clw"))

# Calculate shannon Entropy
entropy_dfs <- lapply(clustal_dfs, calculate_entropy)

# calculate information for each column in MSA
aligned_dfs <- lapply(entropy_dfs, align_scores)

# Read Alphafold Predicted Features
alphafold_folder <- "outputs/AlphaFoldPredictions_Features/"
alphafold_files <- list.files(alphafold_folder)

alphafold_data <- lapply(alphafold_files, function(.x){
	custom_read(.x, 
		path = alphafold_folder, 
		colnames = c("chain", "position", "residue", "sec_struct", "soluble_sa"))})

names(alphafold_data) <- tolower(str_remove(alphafold_files, ".txt"))

alpha_factors <- names(alphafold_data)[names(alphafold_data) %in% names(aligned_dfs)]

# Read Anchor Predicted Features
anchor_folder <- "outputs/AnchorPrediction_Features/"
anchor_files <- list.files(anchor_folder)

anchor_data <- lapply(anchor_files, function(.x){
	custom_read(.x,
		path = anchor_folder,
		colnames = c("position", "residue", "iupred2", "anchor2"))})

names(anchor_data) <- tolower(str_remove(anchor_files, ".txt"))


# temporary fix for tfs that do not have matching alphafold / anchor sequences.
a <- unlist(lapply(alphafold_data, nrow))
b <- unlist(lapply(anchor_data, nrow))
c <- names(a)[a != b-1 & a != b]
# "foxe1"  "foxg1c" "foxl2"  "foxp1b" "tbx1"

# Align AlphaFold and Anchor Predictions
new_factors <- alpha_factors[!(alpha_factors %in% c)]
pred_data <- lapply(new_factors, function(.x){merge_preds(.x, alphafold_data, anchor_data)})
names(pred_data) <- new_factors

# Merge Feature Predictions with MSA information
merged_data = lapply(new_factors, function(.x){merge_data(.x, pred_data, entropy_dfs, aligned_dfs)})
names(merged_data) = new_factors

# Calculate tagging score.
scored_data = lapply(merged_data, calculate_scores)
names(scored_data) = new_factors

# save scores to TSV
for( i in 1:length(scored_data)){
	write_tsv(scored_data[[i]],paste("outputs/scores/", names(merged_data[i]), ".tsv",  sep = ""))
}

#

