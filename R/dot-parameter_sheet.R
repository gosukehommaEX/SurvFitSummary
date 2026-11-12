#' Parameters with Their Covariance Matrix in One Table
#'
#' @description
#' Internal helper of \code{\link{write_survival_workbook}}: the parameter
#' estimates of \code{\link{parameter_table}} on one scale, with the rows
#' of the variance-covariance matrix of each model of
#' \code{\link{covariance_table}} on the same scale in the columns
#' \code{cov_1}, \code{cov_2}, ... (column k is the covariance with the
#' k-th parameter of the model, so the diagonal holds the variances), as in
#' the parameter workbook used by health economic models.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param scale \code{"natural"} or \code{"transformed"}.
#'
#' @return A data frame with the columns \code{fit_type}, \code{arm},
#'   \code{model}, \code{parameter}, \code{est}, \code{lcl}, \code{ucl},
#'   \code{se} and \code{cov_1} to \code{cov_k}, where k is the largest
#'   number of parameters of a model (\code{NA} beyond the parameters of a
#'   model).
#'
#' @keywords internal
#' @noRd
.parameter_sheet <- function(x, scale) {
  par <- parameter_table(x, scale = scale)
  par$scale <- NULL
  if (nrow(par) == 0L) return(par)
  cov <- covariance_table(x, type = "covariance", scale = scale)
  key <- paste(par$fit_type, par$arm, par$model, sep = "\t")
  cov_key <- paste(cov$fit_type, cov$arm, cov$model, sep = "\t")
  n_max <- max(table(key))
  values <- matrix(NA_real_, nrow(par), n_max)
  for (k in unique(key)) {
    rows <- which(key == k)
    pars <- par$parameter[rows]
    v <- cov[cov_key == k, , drop = FALSE]
    wanted <- paste(rep(pars, each = length(pars)),
                    rep(pars, times = length(pars)), sep = "\t")
    found <- v$value[match(wanted, paste(v$row, v$col, sep = "\t"))]
    values[rows, seq_along(rows)] <- matrix(found, length(pars),
                                            length(pars), byrow = TRUE)
  }
  colnames(values) <- paste0("cov_", seq_len(n_max))
  out <- cbind(par, as.data.frame(values))
  rownames(out) <- NULL
  out
}
