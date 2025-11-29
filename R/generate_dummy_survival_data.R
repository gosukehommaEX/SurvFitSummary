#' Generate Dummy Survival Data
#'
#' This function generates simulated time-to-event data with treatment arms.
#' It is useful for testing survival analysis workflows and parametric survival
#' modeling functions.
#'
#' @param arm Character vector specifying treatment arm names.
#'   Example: c("Treatment A", "Treatment B", "Control")
#' @param n Integer vector specifying number of subjects per arm.
#'   Length must match length of \code{arm}. Example: c(50, 100, 150)
#' @param hazards Numeric vector specifying hazard rates for each ARM.
#'   Length must match length of \code{arm}. Higher hazard = shorter survival.
#'   Example: log(2) / c(18, 15, 10) for median survival times of 18, 15, 10 months
#' @param dropout_per_year Numeric value between 0 and 1 specifying annual dropout rate.
#'   Example: 0.1 means 10\% dropout per year
#' @param seed Integer specifying random seed for reproducibility.
#'   Use the same seed to generate identical datasets across runs
#'
#' @return A data frame (tibble) containing the following columns:
#'   \item{SUBJID}{Subject identifier (character)}
#'   \item{ARM}{Treatment arm (factor)}
#'   \item{SURVTIME}{Observed survival time in months (numeric)}
#'   \item{CNSR}{Censoring indicator: 1 = censored, 0 = event (numeric)}
#'
#' @details
#' The function generates survival times from exponential distributions with
#' specified hazard rates for each treatment arm. Dropout times are also
#' generated exponentially based on the annual dropout rate. The observed
#' survival time is the minimum of event time and dropout time.
#'
#' @importFrom dplyr select mutate
#' @importFrom tibble tibble
#' @importFrom stats rexp
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Two-arm trial
#' data1 <- generate_dummy_survival_data(
#'   arm = c("Treatment", "Control"),
#'   n = c(100, 100),
#'   hazards = log(2) / c(20, 15),  # Median survival: 20 vs 15 months
#'   dropout_per_year = 0.1,
#'   seed = 123
#' )
#'
#' # Example 2: Three-arm trial
#' data2 <- generate_dummy_survival_data(
#'   arm = c("Trt_A", "Trt_B", "Control"),
#'   n = c(50, 100, 150),
#'   hazards = log(2) / c(18, 15, 10),
#'   dropout_per_year = 0.1,
#'   seed = 456
#' )
#'
#' # Example 3: Single arm study
#' data3 <- generate_dummy_survival_data(
#'   arm = c("Active"),
#'   n = c(75),
#'   hazards = log(2) / 24,
#'   dropout_per_year = 0.15,
#'   seed = 789
#' )
#' }
generate_dummy_survival_data <- function(arm,
                                         n,
                                         hazards,
                                         dropout_per_year,
                                         seed) {

  # Validate inputs
  if (length(arm) != length(n)) {
    stop("Length of 'arm' must match length of 'n'")
  }

  if (length(arm) != length(hazards)) {
    stop("Length of 'arm' must match length of 'hazards'")
  }

  # Set seed number
  set.seed(seed)

  # Total number of patients
  total_n <- sum(n)

  # Data generation
  dataset <- tibble::tibble(
    SUBJID = paste0('ID', seq(total_n)),
    ARM = rep(arm, n),
    EVENTTIME = unlist(lapply(seq(arm), function(i) stats::rexp(n[i], hazards[i]))),
    DROPTIME = stats::rexp(total_n, rate = -log(1 - dropout_per_year) / 12),
    SURVTIME = pmin(EVENTTIME, DROPTIME),
    CNSR = as.numeric(DROPTIME < EVENTTIME)
  ) %>%
    dplyr::select(SUBJID, ARM, SURVTIME, CNSR)

  # Output
  return(dataset)
}
