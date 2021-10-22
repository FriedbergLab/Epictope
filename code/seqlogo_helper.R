read_clustal <- function(x){read.alignment(paste("outputs/fastas/", x, sep = ""), format = "clustal")}


clustal_df <- function(x){
	list <- x$seq
	df <- dplyr::bind_cols(lapply(list, function(.x){strsplit(.x, split = "")[[1]]}), 
		.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE))
	colnames(df) =  x$nam
	return(df)
}


shannon_entropy <- function(x, nogap = FALSE){
	x <- unlist(x)
	unique_base = unique(x)
	if(nogap == TRUE){unique_base = unique_base[unique_base != "-"]}
	M = length(x)
	entropy_list = list()
	for(i in 1:length(unique_base)){
		base = unique_base[i]
		n_i = table(x)[base]
		P_i = n_i/M
		entropy_i = P_i * log(P_i, base = 2)
		entropy_list[[i]] <- entropy_i
	}
	res <- -sum(unlist(entropy_list))
	return(res)
}

calculate_entropy <- function(.x){
	res <- .x %>%
		rowwise() %>%
		mutate(entropy = shannon_entropy(across(where(is.character)), nogap = TRUE)) %>%
		rowwise()
	return(res)
	}

align_scores <- function(.x){
	num = ncol(.x) - 1
	res <- .x %>% 
		rowid_to_column("id") %>%
		reshape2::melt(id.vars = c("id", "entropy")) %>%
		select(-variable) %>%
		rename("amino" = "value") %>% 
		group_by(id, entropy, amino) %>%
		summarise(value = n()) %>%
		ungroup() %>%
		pivot_wider(names_from = "amino", values_from = "value", values_fill = 0) %>%
		select(-c("-")) %>%
		mutate(info = log(20, base = 2) - entropy) %>%
		select(id, entropy, info, everything() ) %>%
		mutate_at(vars(-one_of("id", "entropy", "info")), ~./num/info)
	}

read_alpha <- function(.x){
	res <- read.table(paste(alphafold_folder, .x, sep = "/"), sep = "\t")
	colnames(res) <- c("chain", "position", "amino", "sec_struct", "soluble_sa")
	return(res)
}

merge_alpha <- function(.x, alphafold_data, processed_dfs, processed_dfs2){
	alpha <- alphafold_data[[.x]]
	msa <- processed_dfs[[.x]] %>%
		rowid_to_column("id") 
	entropy <- processed_dfs2[[.x]]

	msa_seq <- msa %>%
	filter(AlphaFold != "-") %>%
	pull(AlphaFold) %>%
	paste(collapse = "") %>%
	toupper()

	alpha_seq <- paste(alpha$amino, collapse = "")
	ind <- str_locate(alpha_seq, msa_seq)

	a <- alpha %>%
		rename("AlphaFold"= "amino") %>%
		filter(position %in% seq(ind[1], ind[2])) %>%
		select(-position) %>%
		rowid_to_column("position") 
	
	b <- msa %>%
		select(AlphaFold) %>%
		rowid_to_column("id") %>%
		filter(AlphaFold != "-") %>%
		rowid_to_column("position") %>%
		mutate(AlphaFold = toupper(AlphaFold)) %>%
		right_join(a) %>%
		select(id, chain, sec_struct, soluble_sa, AlphaFold) 

	res <- left_join(entropy, b) %>%
		select(id, entropy, info, chain, sec_struct, soluble_sa, AlphaFold, everything()) %>%
		filter(!is.na(chain)) 
}

calculate_scores <- function(.x){
	res <- .x %>%
		select(entropy, chain, sec_struct, soluble_sa, AlphaFold) %>%
		mutate(chain_no = ifelse(chain == "A", 0, 1), struct_no = ifelse(sec_struct == "-", 1, 0)) %>%
		rowwise() %>%
		mutate(score = sum(entropy, chain_no, struct_no, soluble_sa)) %>%
		rowwise()
	return(res)
}