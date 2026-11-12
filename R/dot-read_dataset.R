#' Read a Dataset from a File
#'
#' @description
#' Internal helper of \code{\link{processing_dataset}} that reads a CSV,
#' Excel, SAS or SAS transport file according to its extension.
#'
#' @param path Single character string with the path of the file.
#'
#' @return A data frame with the contents of the file.
#'
#' @keywords internal
#' @noRd
.read_dataset <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }
  ext <- tolower(sub("^.*\\.", "", basename(path)))
  if (!ext %in% c("csv", "xlsx", "xls", "sas7bdat", "xpt")) {
    stop("Unsupported file format: '", ext, "'. Supported formats are ",
         "csv, xlsx, xls, sas7bdat and xpt.", call. = FALSE)
  }
  if (ext %in% c("xlsx", "xls") &&
      !requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is needed to read .", ext, " files.",
         call. = FALSE)
  }
  if (ext %in% c("sas7bdat", "xpt") &&
      !requireNamespace("haven", quietly = TRUE)) {
    stop("Package 'haven' is needed to read .", ext, " files.",
         call. = FALSE)
  }
  out <- switch(
    ext,
    csv = utils::read.csv(path, check.names = FALSE,
                          stringsAsFactors = FALSE),
    xlsx = ,
    xls = readxl::read_excel(path),
    sas7bdat = haven::read_sas(path),
    xpt = haven::read_xpt(path)
  )
  as.data.frame(out)
}
