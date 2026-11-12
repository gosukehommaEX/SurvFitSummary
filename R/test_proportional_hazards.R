#' Grambsch-Therneau Test of Proportional Hazards
#'
#' @description
#' For each treatment arm other than the control arm (the first level of
#' \code{ARM}), fits a Cox model to the data of the two arms and tests the
#' proportional hazards assumption with the score test of Grambsch and
#' Therneau (1994) as implemented in \code{survival::cox.zph()}.
#'
#' @param data Data frame returned by \code{\link{processing_dataset}}, with
#'   two or more treatment arms.
#' @param transform Transformation of time for the test: \code{"km"} (the
#'   default of \code{survival::cox.zph()}), \code{"rank"},
#'   \code{"identity"} or \code{"log"}.
#'
#' @return A data frame with the columns \code{comparison}, \code{log_hr}
#'   (Cox estimate of the log hazard ratio), \code{chisq}, \code{df} and
#'   \code{p_value}.
#'
#' @references
#' Grambsch PM, Therneau TM (1994). Proportional hazards tests and
#' diagnostics based on weighted residuals. Biometrika 81(3):515-526.
#' \doi{10.1093/biomet/81.3.515}
#'
#' @examples
#' raw <- generate_dummy_survival_data(
#'   arm = c("Control", "Treatment"), n = c(100, 100),
#'   hazards = log(2) / c(12, 18), dropout_per_year = 0.1, seed = 123
#' )
#' dat <- processing_dataset(raw, column_arm = "ARM",
#'                           column_survtime = "SURVTIME",
#'                           column_cnsr = "CNSR", control_arm = "Control")
#' test_proportional_hazards(dat)
#'
#' @export
test_proportional_hazards <- function(data,
                                      transform = c("km", "rank", "identity",
                                                    "log")) {
  .check_survival_data(data, min_arms = 2L)
  transform <- match.arg(transform)
  data$ARM <- droplevels(data$ARM)
  arm_levels <- levels(data$ARM)
  rows <- lapply(arm_levels[-1L], function(arm) {
    d <- data[data$ARM %in% c(arm_levels[1L], arm), , drop = FALSE]
    d$ARM <- droplevels(d$ARM)
    fit <- survival::coxph(survival::Surv(SURVTIME, EVENT) ~ ARM, data = d)
    zph <- survival::cox.zph(fit, transform = transform)
    data.frame(
      comparison = paste(arm, "vs", arm_levels[1L]),
      log_hr = unname(stats::coef(fit)[1L]),
      chisq = unname(zph$table["ARM", "chisq"]),
      df = unname(zph$table["ARM", "df"]),
      p_value = unname(zph$table["ARM", "p"])
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
