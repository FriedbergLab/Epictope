read_clustal <- function(x){read.alignment(paste("outputs/fastas/", x, sep = ""), format = "clustal")}


clustal_df <- function(x){
	list <- x$seq
	df <- dplyr::bind_cols(lapply(list, function(.x){strsplit(.x, split = "")[[1]]}), 
		.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE))
	colnames(df) =  x$nam
	return(df)
}

# calculate shannon entropy of each column in MSA
shannon_entropy <- function(x, nogap = FALSE){
	x <- unlist(x)
	unique_base = unique(x)
	if(nogap == TRUE){unique_base = unique_base[unique_base != "-"]; adj_x = x[x != "-"]}
	M = length(adj_x)
	entropy_list = list()
	for(i in 1:length(unique_base)){
		base = unique_base[i]
		n_i = table(adj_x)[base]
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
	options(dplyr.summarise.inform = FALSE)
	res <- .x %>% 
		rowid_to_column("id") %>%
		reshape2::melt(id.vars = c("id", "entropy")) %>%
		select(-variable) %>%
		dplyr::rename("residue"= "value") %>% 
		group_by(id, entropy, residue) %>%
		summarise(value = n()) %>%
		ungroup() %>%
		pivot_wider(names_from = "residue", values_from = "value", values_fill = 0) %>%
		select(-c("-")) %>%
		rowwise() %>%
		mutate(info = log(20, base = 2) - entropy) %>%
		rowwise() %>%
		select(id, entropy, info, everything() ) %>%
		mutate_at(vars(-one_of("id", "entropy", "info")), ~./num/info)
	options(dplyr.summarise.inform = TRUE)
	return(res)
	}

custom_read <- function(.x, path = "", colnames = ""){
	res <- read.table(paste(path, .x, sep = ""), sep = "\t")
	colnames(res) <- colnames
	return(res)
}

calculate_scores <- function(.x){
	res <- .x %>%
		mutate(chain_no = ifelse(chain == "A", 0, 1), struct_no = ifelse(sec_struct == "-", 1, 0.5)) %>%
		mutate(struct_no = ifelse(grepl("[GHIE]", sec_struct), 0, struct_no)) %>%
		rowwise() %>%
		mutate(normalized_entropy = entropy / 4.321928) %>%
		mutate(score = sum(1.5 * normalized_entropy, chain_no, struct_no, soluble_sa, 1-anchor2)) %>%
		rowwise() %>%
		mutate(inv_anchor2 = 1-anchor2) %>%
		select(-anchor2) %>%
		ungroup()

	min <- res %>% 
		select(id, residue, normalized_entropy, struct_no, soluble_sa, inv_anchor2, score) %>%
		group_by(id, residue) %>%
		pivot_longer(cols = c("normalized_entropy", "struct_no", "soluble_sa", "inv_anchor2", "score"), names_to = "var", values_to = "vals") %>%
		slice(which.min(vals)) %>%
		right_join(res, by = c("id", "residue")) %>%
		rename("minimum_var" = var, "min_val" = vals) 
	return(min)
}

merge_preds <- function(.x, alphafold_data, anchor_data){
	alpha <- alphafold_data[[.x]]
	anchor <- anchor_data[[.x]]

	alpha_seq <- paste(alpha[["residue"]], collapse = "")
	anchor_seq <- paste(anchor[["residue"]], collapse = "")
	ind <- str_locate(alpha_seq, anchor_seq)
	if(is.na(sum(ind))){ind = c(1, nchar(alpha_seq))}

	alpha_adj <- alpha %>%
		filter(position %in% seq(ind[1], ind[2])) %>%
		select(-position) %>%
		rowid_to_column("position") 

	res <- left_join(alpha_adj, anchor, by = c("position", "residue"))
	return(res)
}


merge_data <- function(.x, pred_data, entropy_dfs, aligned_dfs){
	pred = pred_data[[.x]]
	msa <- entropy_dfs[[.x]] %>% rowid_to_column("id") 
	#entropy <- aligned_dfs[[.x]]

	msa_seq <- msa %>%
		dplyr::rename("AlphaFold" = contains("ENSDART")) %>%
		filter(AlphaFold != "-") %>%
		pull(AlphaFold) %>%
		paste(collapse = "") %>%
		toupper()

	pred_seq <- paste(pred$residue, collapse = "")
	ind <- str_locate(pred_seq, msa_seq)
	if(is.na(sum(ind))){ind = c(1, nchar(pred_seq))}

	adj_pred <- pred %>%
		filter(position %in% seq(ind[1], ind[2])) %>%
		select(-position) %>%
		rowid_to_column("position") 
	
	adj_msa <- msa %>%
		dplyr::rename("AlphaFold" = contains("ENSDART")) %>%
		select(AlphaFold, entropy) %>%
		rowid_to_column("id") %>%
		filter(AlphaFold != "-") %>%
		rowid_to_column("position") %>%
		mutate(AlphaFold = toupper(AlphaFold)) %>% 
		dplyr::rename("residue" = "AlphaFold") %>%
		right_join(adj_pred, by = c("position", "residue")) %>%
		select(chain, residue, sec_struct, soluble_sa, iupred2, anchor2, entropy) %>% 
		rowid_to_column("id") %>%
		filter(!is.na(chain))
	return(adj_msa)
}
