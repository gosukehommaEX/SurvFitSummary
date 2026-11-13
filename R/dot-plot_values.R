#' Colours or Line Types Given to a Plot
#'
#' @description
#' Internal helper of the plots: checks the colours or line types given in
#' \code{colours} or \code{linetypes} and puts them in place of the defaults.
#' A named vector sets the values of the entries it names, and the other
#' entries keep their defaults; an unnamed vector gives one value per entry,
#' in the order of the entries. Line types are the names used by ggplot2
#' (\code{"solid"}, \code{"dashed"}, \code{"dotted"}, \code{"dotdash"},
#' \code{"longdash"}, \code{"twodash"}, \code{"blank"}) or strings of 2, 4,
#' 6 or 8 hexadecimal digits (for example \code{"44"}).
#'
#' @param values \code{NULL} or the values given by the user.
#' @param entries Names of the entries (arms, curves or elements of the
#'   plot), in the order of the legend.
#' @param defaults Default values, one per entry.
#' @param kind \code{"colour"} or \code{"linetype"}.
#' @param arg Name of the argument, for the error messages.
#'
#' @return A character vector of the values, named by the entries.
#'
#' @keywords internal
#' @noRd
.plot_values <- function(values, entries, defaults, kind, arg) {
  out <- stats::setNames(as.character(defaults), entries)
  if (is.null(values)) return(out)
  what <- if (kind == "colour") "colours" else "line types"
  if (!is.character(values) || length(values) == 0L || anyNA(values)) {
    stop("'", arg, "' must be NULL or a character vector of ", what, ".",
         call. = FALSE)
  }
  valid <- if (kind == "colour") {
    vapply(values, function(v) {
      tryCatch({
        grDevices::col2rgb(v)
        TRUE
      }, error = function(e) FALSE)
    }, logical(1))
  } else {
    values %in% c("blank", "solid", "dashed", "dotted", "dotdash",
                  "longdash", "twodash") |
      grepl("^([1-9A-Fa-f]{2}){1,4}$", values)
  }
  if (!all(valid)) {
    stop("'", arg, "' contains invalid ", what, ": ",
         paste(unique(values[!valid]), collapse = ", "), call. = FALSE)
  }
  given <- names(values)
  if (is.null(given)) {
    if (length(values) != length(entries)) {
      stop("An unnamed '", arg, "' must have ", length(entries),
           " values, one for each of: ", paste(entries, collapse = ", "),
           call. = FALSE)
    }
    out[] <- unname(values)
    return(out)
  }
  if (anyNA(given) || any(given == "") || anyDuplicated(given) > 0L ||
      !all(given %in% entries)) {
    stop("The names of '", arg, "' must be distinct and among: ",
         paste(entries, collapse = ", "), call. = FALSE)
  }
  out[given] <- unname(values)
  out
}
