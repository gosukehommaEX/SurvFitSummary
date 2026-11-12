#' Prepare a Dataset for Survival Analysis
#'
#' @description
#' Takes a data frame or reads a file, renames the columns of the treatment
#' arm, the survival time and the event or censoring indicator to
#' \code{ARM}, \code{SURVTIME} and \code{EVENT}, keeps the requested
#' covariate and stratification columns, and checks the values. The
#' treatment arm is returned as a factor whose first level is the control
#' arm, so that every function of the package uses the same reference arm.
#'
#' @param dataset A data frame, or a single character string with the path
#'   of a CSV (\code{.csv}), Excel (\code{.xlsx}, \code{.xls}), SAS
#'   (\code{.sas7bdat}) or SAS transport (\code{.xpt}) file. Excel files need
#'   the package readxl and SAS files the package haven.
#' @param column_arm Name of the column with the treatment arm.
#' @param column_survtime Name of the column with the survival time. The
#'   times must be positive.
#' @param column_cnsr Name of the column with the censoring indicator
#'   (1 = censored, 0 = event). The default is \code{NULL}.
#' @param column_event Name of the column with the event indicator
#'   (1 = event, 0 = censored). The default is \code{NULL}. When both
#'   \code{column_cnsr} and \code{column_event} are given,
#'   \code{column_event} is used.
#' @param column_covariates Character vector with the names of the columns
#'   kept as covariates for adjusted hazard ratios. The default is
#'   \code{NULL} (none).
#' @param column_strata Character vector with the names of the columns kept
#'   as stratification factors. The default is \code{NULL} (none).
#' @param control_arm Value of the treatment arm used as the control
#'   (reference) arm. The default is \code{NULL}, which uses the first level
#'   when the arm column is a factor and the first value in sorted order
#'   otherwise.
#'
#' @return A data frame with the columns \code{ARM} (factor with the control
#'   arm as the first level), \code{SURVTIME} (numeric), \code{EVENT}
#'   (integer, 1 = event, 0 = censored), followed by the covariate and
#'   stratification columns with their original names. Rows with a missing
#'   value in any of these columns are removed with a warning.
#'
#' @examples
#' raw <- generate_dummy_survival_data(
#'   arm = c("Control", "Treatment"), n = c(50, 50),
#'   hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 1
#' )
#' dat <- processing_dataset(raw, column_arm = "ARM",
#'                           column_survtime = "SURVTIME",
#'                           column_cnsr = "CNSR", control_arm = "Control")
#' head(dat)
#' levels(dat$ARM)
#'
#' @export
processing_dataset <- function(dataset,
                               column_arm,
                               column_survtime,
                               column_cnsr = NULL,
                               column_event = NULL,
                               column_covariates = NULL,
                               column_strata = NULL,
                               control_arm = NULL) {
  if (is.null(column_cnsr) && is.null(column_event)) {
    stop("Either 'column_cnsr' or 'column_event' must be specified.",
         call. = FALSE)
  }
  if (is.character(dataset) && length(dataset) == 1L) {
    dataset <- .read_dataset(dataset)
  }
  if (!is.data.frame(dataset)) {
    stop("'dataset' must be a data frame or the path of a file.",
         call. = FALSE)
  }
  dataset <- as.data.frame(dataset)

  # Columns used, in the order of the output
  column_indicator <- if (!is.null(column_event)) column_event else column_cnsr
  extra <- c(column_covariates, column_strata)
  used <- c(column_arm, column_survtime, column_indicator, extra)
  absent <- setdiff(used, names(dataset))
  if (length(absent) > 0L) {
    stop("Column(s) not found in 'dataset': ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(used) > 0L) {
    stop("The arm, time, indicator, covariate and strata columns must be ",
         "different columns.", call. = FALSE)
  }
  if (any(extra %in% c("ARM", "SURVTIME", "EVENT"))) {
    stop("Covariate and strata columns cannot be named ARM, SURVTIME or ",
         "EVENT.", call. = FALSE)
  }

  # Remove rows with missing values in the columns used
  complete <- stats::complete.cases(dataset[used])
  if (!all(complete)) {
    warning(sum(!complete), " row(s) with missing values were removed.",
            call. = FALSE)
    dataset <- dataset[complete, , drop = FALSE]
  }

  arm <- dataset[[column_arm]]
  time <- as.vector(unclass(dataset[[column_survtime]]))
  indicator <- as.vector(unclass(dataset[[column_indicator]]))
  if (!is.numeric(time)) {
    stop("The survival time column must be numeric.", call. = FALSE)
  }
  if (any(time <= 0)) {
    stop("Survival times must be positive: ", sum(time <= 0),
         " row(s) have a time of 0 or less.", call. = FALSE)
  }
  if (is.logical(indicator)) {
    indicator <- as.integer(indicator)
  }
  if (!is.numeric(indicator) || !all(indicator %in% c(0, 1))) {
    stop("The event or censoring indicator must contain only 0 and 1.",
         call. = FALSE)
  }
  event <- if (!is.null(column_event)) indicator else 1 - indicator

  # Treatment arm as a factor with the control arm first
  if (is.factor(arm)) {
    arm_levels <- levels(droplevels(arm))
  } else {
    arm <- as.vector(unclass(arm))
    arm_levels <- as.character(sort(unique(arm)))
  }
  arm <- as.character(arm)
  if (!is.null(control_arm)) {
    control_arm <- as.character(control_arm)
    if (length(control_arm) != 1L || !control_arm %in% arm_levels) {
      stop("'control_arm' must be one of: ",
           paste(arm_levels, collapse = ", "), call. = FALSE)
    }
    arm_levels <- c(control_arm, setdiff(arm_levels, control_arm))
  }

  out <- data.frame(
    ARM = factor(arm, levels = arm_levels),
    SURVTIME = as.numeric(time),
    EVENT = as.integer(event)
  )
  for (column in extra) {
    x <- dataset[[column]]
    # Drop SAS labels and formats, keep factors as they are
    out[[column]] <- if (is.factor(x)) x else as.vector(unclass(x))
  }
  rownames(out) <- NULL
  out
}
