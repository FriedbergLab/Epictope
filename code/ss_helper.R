# searches the NCBI FTP site for a given species
ftp_search <- function(species){
    # list of ncbi species
    webpage <- read_html("https://ftp.ensembl.org/pub/release-108/fasta/")
    # extract links
    links <- html_nodes(webpage, "a") %>% html_attr("href")
    link_names <- html_nodes(webpage, "a") %>% html_text()
    names(links) <- link_names
    # search for species
    species_regex <- paste0(species, "/$")
    species_query <- paste(paste0(species, "/$"), collapse = "|")
    species_links <- links[base::grepl(species_query, links, perl = TRUE)]
    not_found <- setdiff(species, gsub("/", "", species_links))
    if (length(not_found) > 0) {
        message(paste("The following species were not found:", paste(not_found, collapse = ", ")))
    }
    # return links
    return(species_links)
}

# downloads the CDS fasta files from the NCBI FTP site
ftp_download <- function(link){
  # construct the URL for the CDS directory of the given link
  link_url <- paste0("https://ftp.ensembl.org/pub/release-108/fasta/", link, "cds/")
  # read webpage
  webpage <- rvest::read_html(link_url)
  # extract links
  cds_link <- rvest::html_nodes(webpage, "a") %>% rvest::html_attr("href")
  # filter to ".fa.gz" extension
  filename <- cds_link[grepl("\\.fa\\.gz", cds_link)]
  # download
  tryCatch(
        {
            download.file(url = paste0(link_url, filename), destfile = paste0(cds_folder,"/", filename), method = "wget", extra = "-nc")
        },
        # if the download fails, search for the file in the specified directory and return if found
        error = function(cond) {
            return(NULL)
        }
    )
}


# make protein blastable database from cds file
make_protein_db <- function(file){
    # read in fasta file
    fasta <- Biostrings::readDNAStringSet(file, format="fasta")
    # replace cds
    names(fasta) <- gsub(" cds.*", "", names(fasta), perl = TRUE)
    # convert to AA
    prot <- Biostrings::translate(fasta, if.fuzzy.codon = "solve")
    # write to file 
    amino <-  gsub(".fa$", ".amino.fa", file)
    Biostrings::writeXStringSet(prot, amino, append = TRUE)
    # make protein database of AA file.
    system(command = paste("makeblastdb -in ", amino, " -dbtype prot"))
    #rBLAST::makeblastdb(cur_amino, dbtype = "prot", args = "")
}


# define a function to access the UniProt API
query_uniProt <- function (query = character(0L), fields = c("accession", "id"), 
    collapse = " OR ") 
{
    # check that query and fields are both character vectors
    stopifnot(is.character(query), is.character(fields))
    
    # if the length of query is zero, throw an error message
    if (!length(query)) 
        stop("<internal> 'qlist' must be populated with queries")
    
    # send an HTTP GET request to the UniProt API, passing the query and field parameters
    # and specifying the desired response format
    resp <- httr::GET(paste0("https://rest.uniprot.org/", "uniprotkb/search"), 
        query = list(query = paste(query, collapse = collapse), 
            fields = paste(fields, collapse = ","), format = "tsv"))
    
    # read the response text into a data frame using the read.delim function and return 
    read.delim(text = httr::content(resp, encoding = "UTF-8"))
}

# define a function to fetch Alphafold2 protein structure predictions for a given protein
fetch_alphafold <- function(.x) {
    # Check that .x is a character vector
    stopifnot(is.character(.x))
    # set the base URL for downloading Alphafold2 protein structure predictions
    base_url <- "https://alphafold.ebi.ac.uk/files/"
    # set the version suffix 
    version_suffix <- "-F1-model_v4.pdb"
    # define output filename
    filename <- paste0("AF-", .x, version_suffix)
    # add path 
    output_file <- paste0(model_folder, "/", filename)
    # define download url
    download_url <- paste0(base_url, filename)
    
    # try to download the Alphafold2 protein structure prediction using the specified URL and save it to the output file
    res <- tryCatch(
        {
            download.file(download_url, destfile = output_file, method = "wget", extra = "-nc")
            message("Downloading Alphafold2 PDB for ", .x)
            return(output_file)
        },
        
        # if the download fails, search for the file in the specified directory and return if found
        error = function(cond) {
            message(paste("Error while downloading: ", .x))
            message(cond)
            message("\nSearching ", model_folder, " for existing ", filename, ".")
            cur_file <- list.files(path = model_folder, pattern = filename, recursive = TRUE, full.names = TRUE)
            if (!is.na(cur_file)) {
                message(cur_file, " found, returning.")
            }
            return(cur_file)
        }
    )
    return(res)
}

# run dssp command from R
dssp_command <- function(.x){
    # Check that .x is a character vector
    stopifnot(is.character(.x))
    # define output filename
    output_file <- gsub("\\.pdb.gz|\\.pdb", ".dssp", .x)
    # run dssp command with mkdssp, return error message if error.
    system(command = paste0("mkdssp -i ", .x, " -o ", output_file), intern = TRUE)
    # return name of output file.
    return(output_file)
}

# parse dssp file
parse_dssp <- function(file, keepfiles = FALSE){
  ## https://rdrr.io/github/jcaledo/ptm_0.1.1/src/R/dssp.R
  ## --------------- Reading the dssp file ------------------ ##
  con <- file(file, 'r')

  counter <- 0
  resnum <- c()
  respdb <- c()
  chain <- c()
  aa <- c()
  ss <- c()
  sasa <- c()
  phi <- c()
  psi <- c()

  while(TRUE){
    line <- readLines(con, n = 1)
    counter <- counter + 1

    if (counter == 1){
      l <- strsplit(line, split = "")[[1]]
      l <- paste(l, collapse = "")
      if ("have bz2" %in% l){
        first_valid_line <- 29 # dssp file coming from the API
      } else {
        first_valid_line <- 28 # dssp file coming from the sync
      }
    }

    if (counter > first_valid_line & length(line) != 0){
      a <- strsplit(line, split = "")[[1]]
      resnum <- c(resnum, paste(a[1:5], collapse = ""))
      respdb <- c(respdb, paste(a[6:10], collapse = ""))
      chain <- c(chain, paste(a[11:12], collapse = ""))
      aa <- c(aa, paste(a[13:14], collapse = ""))
      ss <- c(ss, paste(a[15:17], collapse = ""))
      sasa <- c(sasa, paste(a[36:38], collapse = ""))
      phi <- c(phi, paste(a[104:109], collapse = ""))
      psi <- c(psi, paste(a[110:115], collapse = ""))
    }
    if (length(line) == 0){
      break
    }
  }
  close(con)

  ## ------ Setting the variable types ------------- ##
  resnum <- as.numeric(resnum)
  respdb <- as.numeric(respdb)
  chain <- gsub(" ", "", chain)
  aa <- gsub(" ", "", aa)
  ss <- gsub("   ", "C", ss)
  ss <- gsub(" ", "", ss)

  ## -------- Building the dataframe ---------------- ##
  df <- as.data.frame(matrix(c(resnum, respdb, chain, aa,
                               ss, sasa, phi, psi), ncol = 8),
                      stringsAsFactors = FALSE)

  colnames(df) <- c('resnum', 'position', 'chain', 'aa', 'ss',
                    'sasa', 'phi', 'psi')

  df$resnum <- as.numeric(df$resnum)
  df$position <- as.character(df$position)
  df$sasa <- as.numeric(df$sasa)
  df$phi <- as.numeric(df$phi)
  df$psi <- as.numeric(df$psi)

  if (keepfiles == TRUE){
    save(df, file = paste(file, ".Rda", sep = ""))
  } else {
    file.remove(file)
  }

  ## --------------- Remove empty lines between chains ------------- ##
  badlines <- c()
  for (i in 1:nrow(df)){
    if (df$aa[i] == '!' | df$aa[i] == 'X'){
      badlines <- c(badlines, i)
    }
  }
  if (length(badlines) != 0){
    df <- df[-badlines,]
    df$resnum <- 1:nrow(df)
  }
  return(df)
}

# run blast on input sequence against blast file 
protein_blast <- function(seqs, blast_file){
    # Check that blast_file is a character vector
    stopifnot(is.character(blast_file))
    db <- rBLAST::blast(db = blast_file, type = "blastp") 
    blast_res <- predict(db, seqs) 
    blast_res$filename = blast_file
    return(blast_res)
}

# fetch query sequence from file
fetch_sequences <- function(.x){
    # read in fasta file
    cur_proteins <- Biostrings::readAAStringSet(.x$filename, format = "fasta")
    # subset by query name
    cur_seqs <- cur_proteins[.x$SubjectID]
    return(cur_seqs)
}

# calculate shannon entropy
shannon_entropy <- function(x, nogap = FALSE){
    # ensure that x is a vector, not a list.
    x <- unlist(x) 
    # take the unique bases in the sequence.
    unique_base = unique(x)
    if(nogap == TRUE){
        unique_base = unique_base[unique_base != "-"] # exclude gaps if specified.
        adj_x = x[x != "-"]
    }else{
        adj_x = x
    }
    # get the length of the adjusted sequence.
    M = length(adj_x) 
    # initialize a list to store the entropy for each base.
    entropy_list = list() 
    for(i in 1:length(unique_base)){
    # take the current base.
    base = unique_base[i] 
    # count the number of occurrences.
    n_i = table(adj_x)[base] 
    # calculate the probability.
    P_i = n_i/M 
    # calculate the entropy.
    entropy_i = P_i * log(P_i, base = 2)
    # add entropy to list. 
    entropy_list[[i]] <- entropy_i 
  }
  # calculate the total entropy by summing the entropies for each base.
  res <- -sum(unlist(entropy_list)) 
  # return the result.
  return(res) 
}

# calculate shannon entropy on each column of query MSA
shannon_reshape <- function(msa, query){
  # extract msa as matrix 
  mat <- as.matrix(msa@unmasked)
  # determine row for query sequence
  chars <- mat[query,]
  # calculate conservation on query positions
  ind <- which(chars != "-")
  # calculate shannon entropy
  score <- apply(mat, 2, shannon_entropy, nogap = FALSE)
  # return shannon entropy for only query positions
  temp <- score[ind]
  # write result to dataframe.
  res <- data.frame(position = 1:length(temp), aa = chars[ind], shannon = temp)
  return(res)
}


# define a function called ss_convert that takes one argument, .x
ss_convert <- function(.x) {
  # check if .x is a recognized Secondary Structure Symbol
  if (!grepl("[GHIECTBS]", .x)) {
    # if not, print an error message and return NULL
    message("Character is not a recognized Secondary Structure Symbol.")
    return(NULL)
  } else {
    return(ss_key[.x])
  }
}

# normalize the tagging features and calculate final score
calculate_scores <- function(.x){

    # normalie entropy by dividing by max
    normalized_entropy <- .x$shannon / 4.321928
    # convert secondary structure symbols to numeric value
    ss_score <- as.numeric(sapply(.x$ss, ss_convert))
    # calculate relative solvent accesibility
    rsa <- .x$sasa / max_sasa[.x$aa]
    # invert the anchor2 score
    inv_anchor2 <- 1-.x$anchor2

    # calculate final score
    tag_score <- (h_weight * normalized_entropy)  + (ss_weight * ss_score) + (rsa_weight * rsa) + (br_weight * inv_anchor2)

    # bind scores to dataframe.
    res <- data.frame(
        position = 1:nrow(.x),
        normalized_entropy = normalized_entropy, 
        ss_score = ss_score, 
        rsa = rsa, 
        inv_anchor2 = inv_anchor2, 
        tag_score = tag_score) 
    # determine minimum feature.
    res$min <- apply(res, 1, FUN = function(.x){return(min(.x))})
    res$min_feature <- apply(res, 1, FUN = function(.x){colnames(res)[which(.x == min(.x))]})
    
    return(res)
}

