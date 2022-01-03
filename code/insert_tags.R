library(tidyverse)
library(seqinr)
library(Biostrings)

tfs <- c("hand2", "snai2", "smad5", "tcf21", "tbx20", "hey2")

#score data
score_folder <- "outputs/scores"
files = list.files(score_folder, "*.tsv")
score_data <- list()
for( i in 1:length(files)){
	score_data[[i]] = read_tsv(paste(score_folder, files[i], sep = "/"))
}
names(score_data) <- sapply(files, function(.x){strsplit(.x, "\\.")[[1]][1]})

#cDNA fastas
cDNA_folder <- "outputs/cDNA"
files = list.files(cDNA_folder, "*.fasta")
cDNA_data <- list()
for( i in 1:length(files)){
	cDNA_data[[i]] = read.fasta(paste(cDNA_folder, files[i], sep = "/"))[[1]] %>%
		paste0(collapse = "") %>%
		as.character() %>%
		toupper()
}
names(cDNA_data) <- sapply(files, function(.x){strsplit(.x, "\\.")[[1]][1]})
#https://resources.chromotek.com/blog/v5-tag-properties-origin-specifications-tips
V5_tag <- "GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG"
#https://signagen.com/blog/2016/01/20/common-epitope-tags/
HA_tag <- "TACCCATACGATGTTCCAGATTACGCT"

#Identify site locations to insert tag
site_locations <- list()
for(i in 1:length(tfs)){
	# get score data for tf of interest
	# sorte by desc min_value score
	temp <- score_data[[tfs[i]]] %>%
		arrange(desc(min_val))  
	# pulls residue positions (in AA) 
	ids = temp %>% pull(id)
	# take best insertion position
	res = ids[1]
	# look over the rest of the positions and check if they are within range of the tag insertion
	# add to list of positions if not
	for(j in 2:length(ids)){
		if(sum(abs(res - ids[j]) < (nchar(V5_tag)/3)) < 1){res <- c(res, ids[j])}
	}
	# filter score df to only positions that are acceptable tag positions
	ranks <- temp %>%
		rowid_to_column("rank") %>%
		rowwise() %>%
		mutate(tag = ifelse(id %in% res, TRUE, FALSE)) %>%
		arrange(tag, desc(min_val))
	# put tag sites into site_locations list.
	site_locations[[i]] <- ranks
}
names(site_locations) <- names(cDNA_data)

tag_list <- list()
for(i in 1:length(cDNA_data)){
	message(nchar(cDNA_data[[i]]))
	ids = site_locations[[names(cDNA_data)[i]]] %>% 
		filter(tag == TRUE) %>% 
		arrange(desc(min_val)) %>%
		pull(id)
	V5_list <- list()
	HA_list <- list()
	for(j in 1:length(ids)){
		temp_HA <- temp_V5 <- cDNA_data[[i]]
		pos = ids[j] * 3
		stringi::stri_sub(temp_V5, pos+1, pos) <- V5_tag
		stringi::stri_sub(temp_HA, pos+1, pos) <- HA_tag
		V5_list[[j]] <- temp_V5
		HA_list[[j]] <- temp_HA
	}
	names(V5_list) <- paste(names(cDNA_data)[i], "V5", as.character(ids*3))
	names(HA_list) <- paste(names(cDNA_data)[i], "HA", as.character(ids*3))
	tag_list[[i]] <- c(V5_list, HA_list)
}
names(tag_list) <- names(cDNA_data)

for(i in 1:length(tag_list)){
	write.fasta(tag_list[[i]], names(tag_list[[i]]), paste(names(tag_list)[i], "_cDNA.fasta", sep = ""))
}

for(i in 1:length(tag_list)){
	temp <- lapply(tag_list[[i]], function(.x){as.character(Biostrings::translate(DNAStringSet(.x)))})
	names_AA <- unlist(lapply(names(temp), function(.x){
		a <- unlist(str_split(.x, pattern = " "))
		b <- as.numeric(a[3])/3
		res <- paste(a[1], a[2], b)
		return(res)
	}))
	write.fasta(temp, names_AA, paste(names(tag_list)[i], "_AA.fasta", sep = ""))
}

for(i in 1:length(site_locations)){
	write_csv(site_locations[[i]], paste(names(site_locations)[i], "_tagged_scores.csv", sep = ""))
}

hand2_cDNA <- cDNA_data[[1]]
V5_tag <- "GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG"
lapply(temp, function(.x){grepl("GKPIPNPLLGLDST", .x)})
lapply(temp, function(.x){grepl("YPYDVPDYA", .x)})


lapply(names_AA, function(.x){
		a <- unlist(str_split(.x, pattern = " "))
		b <- as.numeric(a[3])/3
		res <- paste(a[1], a[2], b)
		return(res)
	})