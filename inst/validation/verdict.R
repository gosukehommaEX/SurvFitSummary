# Verdicts of the validation and reproduction scripts.
#
# Sourced by inst/validation/validate_*.R and inst/reproduce/reproduce_*.R
# (run from the root of the package after devtools::load_all()).
#
# verdict_table() compares the values of the package with reference values
# (published, closed form or another package): PASS if the absolute
# difference is within the tolerance, EXPLAINED if not but the item is
# registered with its reason in inst/validation/explained.csv and the
# difference is at most the max_difference registered with it (the largest
# difference the reason accounts for), FAIL otherwise, and INFO for values
# reported without a reference (a missing value or tolerance with a
# reference FAILs). (Without this bound, a registered item was
# EXPLAINED whatever its difference, also in the mutation run of
# 2026-11-13.)
#
# Mutation run: with the environment variable SURVFITSUMMARY_MUTATION set to
# "1", the values of the package are shifted by 0.01 * max(1, |value|)
# before the comparison, and the results are written to results/mutation/.
# Every comparison must then FAIL, which shows that the comparisons are not
# vacuous (inst/validation/run_validation.R checks this). This holds if
# every tolerance is at most 0.005 * max(1, |reference|), half the shift,
# and every max_difference below 0.005, half the smallest shift, which
# verdict_table() checks, so that a comparison that cannot detect an error
# of that size stops the normal run. (A shift of 1% + 0.01 left values near
# -1 almost unchanged, so that 3 comparisons passed in the mutation run of
# 2026-11-13; tolerances in standard errors larger than the shift let 71
# comparisons of validate_python_reference.R pass in that of the same day.)

mutation_run <- function() {
  identical(Sys.getenv("SURVFITSUMMARY_MUTATION"), "1")
}

verdict_table <- function(script, item, reference, value, tolerance) {
  n <- length(item)
  reference <- rep_len(as.numeric(reference), n)
  value <- rep_len(as.numeric(value), n)
  tolerance <- rep_len(as.numeric(tolerance), n)
  too_wide <- !is.na(reference) & !is.na(tolerance) &
    tolerance > 0.005 * pmax(1, abs(reference))
  if (any(too_wide)) {
    stop(script, ": the tolerance of ", sum(too_wide), " comparison(s) ",
         "exceeds half the shift of the mutation run, 0.005 * max(1, ",
         "|reference|), for example '", item[too_wide][1L], "'.",
         call. = FALSE)
  }
  if (mutation_run()) {
    value <- value + 0.01 * pmax(1, abs(value))
  }
  explained <- utils::read.csv(file.path("inst", "validation",
                                         "explained.csv"),
                               stringsAsFactors = FALSE)
  if (!is.numeric(explained$max_difference) ||
      anyNA(explained$max_difference) ||
      any(explained$max_difference <= 0 |
            explained$max_difference >= 0.005)) {
    stop("explained.csv: max_difference must be a number in (0, 0.005).",
         call. = FALSE)
  }
  explained <- explained[explained$script == script, , drop = FALSE]
  bound <- explained$max_difference[match(item, explained$item)]
  difference <- abs(value - reference)
  verdict <- ifelse(is.na(reference), "INFO",
                    ifelse(!is.na(difference) & !is.na(tolerance) &
                             difference <= tolerance,
                           "PASS",
                           ifelse(!is.na(bound) & !is.na(difference) &
                                    difference <= bound,
                                  "EXPLAINED", "FAIL")))
  data.frame(script = script, item = item, reference = reference,
             value = value, difference = difference, tolerance = tolerance,
             max_difference = bound, verdict = verdict,
             stringsAsFactors = FALSE)
}

write_verdicts <- function(result, folder) {
  if (mutation_run()) folder <- file.path(folder, "mutation")
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  script <- result$script[1L]
  path <- file.path(folder, paste0(script, ".csv"))
  utils::write.csv(result, path, row.names = FALSE)
  print(table(factor(result$verdict,
                     levels = c("PASS", "EXPLAINED", "FAIL", "INFO"))))
  invisible(path)
}

# Half a unit of the last printed decimal: tolerance of a published value
printed_tolerance <- function(decimals) {
  0.5 * 10^(-decimals) + 1e-9
}
