#' @title parse_dssp
#'
#' @description Parses a DSSP file and returns the data in a data frame.
#'
#' @param file A character vector specifying the path to the input DSSP file.
#' @param keepfiles A logical value indicating whether to keep the parsed data in an Rda file (default is FALSE).
#'
#' @return The function returns a data frame containing the parsed DSSP data.
#'
#' @examples
#' #parse_dssp("protein.dssp", keepfiles = TRUE)
#' @export

# https://rdrr.io/github/jcaledo/ptm_0.1.1/src/R/dssp.R
parse_dssp <- function(file, keepfiles = TRUE) {
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

  while (TRUE) {
    line <- readLines(con, n = 1)
    counter <- counter + 1

    if (counter == 1) {
      l <- strsplit(line, split = "")[[1]]
      l <- paste(l, collapse = "")
      if ("have bz2" %in% l) {
        first_valid_line <- 28  # DSSP file coming from the API
      } else {
        first_valid_line <- 27  # DSSP file coming from the sync
      }
    }

    if (counter > first_valid_line & length(line) != 0) {
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
    if (length(line) == 0) {
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

  if (keepfiles == TRUE) {
    save(df, file = paste(file, ".Rda", sep = ""))
  } else {
    file.remove(file)
  }

  ## --------------- Remove empty lines between chains ------------- ##
  badlines <- c()
  for (i in 1:nrow(df)) {
    if (df$aa[i] == '!' | df$aa[i] == 'X') {
      badlines <- c(badlines, i)
    }
  }
  if (length(badlines) != 0) {
    df <- df[-badlines, ]
    df$resnum <- 1:nrow(df)
  }
  return(df)
}
