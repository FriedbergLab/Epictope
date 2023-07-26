#' @title setup_files
#'
#' @description Creates the required folder structure..
#'
#' @examples
#' setup_files()
#' @export

setup_files <- function() {
  # folder structure
  dataFolder <- "data"
  outputFolder <- "outputs"
  model_folder <- file.path(dataFolder, "models")
  cds_folder <- file.path(dataFolder, "CDS")
  temp_folder <- file.path(dataFolder, "temp")

  folders <- c(dataFolder, outputFolder, model_folder, cds_folder, temp_folder)

  for (folder in folders) {
    if (!dir.exists(folder)) {
      dir.create(folder)
    }
  }

  assign("dataFolder", dataFolder, envir = .GlobalEnv)
  assign("outputFolder", outputFolder, envir = .GlobalEnv)
  assign("model_folder", model_folder, envir = .GlobalEnv)
  assign("cds_folder", cds_folder, envir = .GlobalEnv)
  assign("temp_folder", temp_folder, envir = .GlobalEnv)
}
