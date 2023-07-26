#' @title findExecutable
#'
#' @description Find the location of an executable in the system's PATH environment variable.
#'
#' @param exe A character string specifying the name of the executable to be found.
#' @param interactive A logical value indicating whether to stop with an error message if the executable is not found.
#'
#' @return The function returns the path to the executable or an empty character vector if the executable is not found and interactive is set to FALSE.
#'
#' @examples
#' # Find the location of the "muscle" executable
#' #muscle_path <- .findExecutable("muscle")
#' #muscle_path
#' @export

.findExecutable <- function(exe, interactive = TRUE) {
  path <- Sys.which(exe)

  if (all(path == "")) {
    if (interactive)
      stop(
        "Executable for ",
        paste(exe, collapse = " or "),
        " not found! Please make sure that the software is correctly installed and, if necessary, path variables are set.",
        call. = FALSE
      )
    return(character(0))
  }

  path[which(path != "")[1]]
}
