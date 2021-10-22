library(seqinr)
library(ggseqlogo)
library(tidyverse)
source("seqlogo_helper.R")


clustal_folder <- "outputs/fastas/"
clustal_files <- list.files(clustal_folder, pattern = ".clw")
clustals = lapply(clustal_files, read_clustal)

dfs <- lapply(clustals, clustal_df)
names(dfs) = stringr::str_remove(clustal_files, ".clw")

processed_dfs <- lapply(dfs, calculate_entropy)
processed_dfs2 <- lapply(processed_dfs,align_scores)

# Read Alphafold Predicted Features
alphafold_folder <- "outputs/AlphaFoldPredictions_Features"
alphafold_files <- list.files(alphafold_folder)

alphafold_data <- lapply(alphafold_files, read_alpha)
names(alphafold_data) <- str_remove(alphafold_files, ".txt")

factors <- names(alphafold_data)[names(alphafold_data) %in% names(processed_dfs2)]

merged_data = lapply(factors, function(.x){merge_alpha(.x, alphafold_data, processed_dfs, processed_dfs2)})
names(merged_data) = factors

scored_data = lapply(merged_data, calculate_scores)
names(scored_data) = names(merged_data)

for( i in 1:length(merged_data)){
	temp <- scored_data[[i]] %>%
		select(entropy, chain, sec_struct, soluble_sa, chain_no, struct_no, score) 
	write_tsv(temp,paste("outputs/scores/", names(merged_data[i]), ".tsv",  sep = ""))
}

