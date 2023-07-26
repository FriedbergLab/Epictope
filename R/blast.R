#' @title blast
#'
#' @description Create a BLAST object for working with the BLAST database.
#'
#' @param db A character string specifying the path to the BLAST database.
#' @param type A character string specifying the type of BLAST (e.g., "blastp", "blastn", "blastx", etc.).
#'
#' @return The function returns an object of class "BLAST" containing information about the BLAST database.
#' @export


#   https://github.com/mhahsler/rBLAST
#   Copyright (C) 2012 Michael Hahsler and Anurag Nagar
blast <- function(db = NULL, type = "blastp") {
  # Check if db is specified
  if (is.null(db))
    stop("No BLAST database specified!")

  db <- file.path(normalizePath(dirname(db)), basename(db))

  if (length(Sys.glob(paste(db, "*", sep = ""))) < 1)
    stop("BLAST database does not exist! (tried to open: ", db, ")")

  # Check for spaces in input
  if (length(grep(" ", db)) > 0)
    stop(
      "Database name or path contains spaces. Rename or move the database to remove spaces (current path: ",
      db,
      ")"
    )

  # Check if executable is available
  .findExecutable(type)

  # Check database
  status <- try(system2(
    .findExecutable("blastdbcmd"),
    args = c("-db", db, "-info"),
    stdout = FALSE
  ))

  if (status != 0)
    stop("Problem loading the database! (trying to execute: blastdbcmd)")

  structure(list(db = db, type = type), class = "BLAST")
}

print.BLAST <- function(x, info = TRUE, ...) {
  cat("BLAST Database\nLocation:", x$db, "\n")
  cat("BLAST Type:", x$type, "\n")

  if (info) {
    out <- system2(
      .findExecutable("blastdbcmd"),
      args = c("-db", x$db,
               "-info"),
      stdout = TRUE
    )
    cat(paste(out, collapse = "\n"))
    cat("\n")
  }
}

predict.BLAST <-
  function(object, newdata, silent = FALSE, BLAST_args = "", custom_format = "", ...) {
    # assigne db, executable, and data to new local variables
    db <- object$db
    exe <- object$type
    x <- newdata

    # get temp files and change working directory
    wd <- tempdir()
    dir <- getwd()
    temp_file <- basename(tempfile(tmpdir = wd))
    on.exit({
      file.remove(Sys.glob(paste(temp_file, "*", sep = "")))
      setwd(dir)
    })
    setwd(wd)
    # assign input and output files
    infile <- paste(temp_file, ".fasta", sep = "")
    outfile <- paste(temp_file, "_BLAST_out.txt", sep = "")
    # write sequences to file
    writeXStringSet(x, infile, append = FALSE, format = "fasta")
    # run BLAST with system call
    system2(
      command = .findExecutable(exe),
      args = c("-db", db, "-query", infile, "-out", outfile, '-outfmt "10', custom_format, '"', BLAST_args)
    )

    # rdp output column names
    if (custom_format == "") {
      c_names <- c("QueryID", "SubjectID", "Perc.Ident","Alignment.Length","Mismatches","Gap.Openings",
                   "Q.start","Q.end","S.start","S.end","E","Bits")
    } else{
      c_names <- unlist(strsplit(custom_format, split = " +"))
    }

    # read and parse BLAST output
    if (is(try(cl_tab <-
               read.table(outfile, sep = ",", quote = ""), silent = silent)
           , "try-error")) {
      warning("BLAST did not return a match!")
      cl_tab <- data.frame(matrix(ncol = length(c_names), nrow = 0))
    }
    colnames(cl_tab) <- c_names
    return(cl_tab)
  }


