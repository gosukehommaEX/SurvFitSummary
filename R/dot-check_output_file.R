#' Check the Path of an Output File
#'
#' @description
#' Internal helper that checks that \code{file} is a single path with the
#' given extension in an existing folder.
#'
#' @param file Path of the file to write.
#' @param extension Extension without the dot, for example \code{"xlsx"}.
#'
#' @return \code{NULL}, invisibly; an error if the path is not valid.
#'
#' @keywords internal
#' @noRd
.check_output_file <- function(file, extension) {
  if (!is.character(file) || length(file) != 1L || is.na(file) ||
      !grepl(paste0("\\.", extension, "$"), file, ignore.case = TRUE)) {
    stop("'file' must be the path of a .", extension, " file.",
         call. = FALSE)
  }
  if (!dir.exists(dirname(file))) {
    stop("The folder of 'file' does not exist: ", dirname(file),
         call. = FALSE)
  }
  invisible(NULL)
}
