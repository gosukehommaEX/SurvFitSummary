#' Data Overview and Hazard Ratios in the Word Report
#'
#' @description
#' Internal helper: the start of the Results section of the report,
#' with the data summary, the Kaplan-Meier curves (all arms in
#' one plot and each arm alone, each with and/or without confidence
#' intervals), the Kaplan-Meier table of each arm and the hazard ratios.
#'
#' @param x Object returned by \code{\link{fit_survival_models}}.
#' @param s List returned by \code{.report_settings()}.
#' @param conf_int Logical vector: versions of the figures without
#'   (\code{FALSE}) and/or with (\code{TRUE}) confidence intervals.
#'
#' @return A list of elements made by \code{.report_block()}.
#'
#' @keywords internal
#' @noRd
.report_data_results <- function(x, s, conf_int) {
  b <- .report_block
  level <- paste0(format(100 * x$conf_level), "%")
  ci_text <- function(ci) {
    if (ci) paste0(", with ", level, " confidence intervals") else ""
  }
  two_arms <- length(x$arms) >= 2L
  sm <- summarise_data(x$data, conf_level = x$conf_level)
  total <- if (two_arms) sm[sm$ARM == "Overall", ] else sm[1L, ]
  arms <- sm[sm$ARM != "Overall", , drop = FALSE]
  median_text <- paste0(
    .format_number(arms$median, 1L, na = "NR"), " (",
    .format_number(arms$median_lcl, 1L, na = "NR"), ", ",
    .format_number(arms$median_ucl, 1L, na = "NR"), ")"
  )
  summary_table <- data.frame(
    "Treatment" = as.character(arms$ARM),
    "N patients" = as.character(arms$n),
    "N events" = as.character(arms$events),
    "N censored" = as.character(arms$censored),
    "Median" = median_text,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  names(summary_table)[5L] <- paste0("Median time to event, ", s$unit,
                                     " (", level, " CI)")
  blocks <- list(
    b("heading", "Results", level = 1L),
    b("heading", "Data overview", level = 2L),
    b("heading", "Data Summary", level = 3L),
    b("paragraph", paste0(
      "Total observations: ", total$n, ". Total events: ", total$events,
      ". Total censored: ", total$censored, "."
    )),
    b("table", "Data summary", table = summary_table),
    b("paragraph", "Abbreviations: CI, confidence interval; NR, not reached."),
    b("heading", "Kaplan-Meier curves and tables", level = 3L)
  )
  km_args <- function(arms, ci) {
    list(fun = plot_km,
         args = list(data = x$data, arms = arms, conf_int = ci,
                     conf_level = x$conf_level, time_unit = x$time_unit))
  }
  if (two_arms) {
    for (ci in conf_int) {
      blocks <- c(blocks, list(b(
        "figure", paste0("Kaplan-Meier curves of all arms", ci_text(ci)),
        plot = km_args(NULL, ci), height = 6.5
      )))
    }
  }
  for (a in x$arms) {
    for (ci in conf_int) {
      blocks <- c(blocks, list(b(
        "figure", paste0("Kaplan-Meier curve of ", a, ci_text(ci)),
        plot = km_args(a, ci), height = 6
      )))
    }
  }
  km <- km_table(x$data, interval = s$interval, conf_level = x$conf_level)
  for (a in x$arms) {
    k <- km[km$ARM == a, , drop = FALSE]
    table <- data.frame(
      "Time" = as.character(round(k$time, 4L)),
      "Number at risk" = as.character(k$n_risk),
      "Number of events" = as.character(k$n_event),
      "Number of censors" = as.character(k$n_censor),
      "Survival probability (SE)" = paste0(.format_number(k$surv, 3L), " (",
                                           .format_number(k$se, 3L), ")"),
      "CI" = paste0(.format_number(k$lcl, 3L), ", ",
                    .format_number(k$ucl, 3L)),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    names(table)[1L] <- paste0("Time (", s$unit, ")")
    names(table)[6L] <- paste("Survival probability", level, "CI")
    blocks <- c(blocks, list(
      b("table", paste("Kaplan-Meier estimates of", a), table = table)
    ))
  }
  if (two_arms) {
    hr <- s$hazard_ratios
    hr_table <- data.frame(
      "Model" = hr$model,
      "Comparison" = hr$comparison,
      "Covariates" = hr$covariates,
      "Strata" = hr$strata,
      "HR" = paste0(.format_number(hr$hr, 2L), " (",
                    .format_number(hr$lcl, 2L), ", ",
                    .format_number(hr$ucl, 2L), ")"),
      "p-value" = .format_p(hr$p_value),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    names(hr_table)[5L] <- paste0("HR (", level, " CI)")
    blocks <- c(blocks, list(
      b("heading", "Hazard ratios", level = 3L),
      b("table", "Hazard ratios from Cox proportional hazards models",
        table = hr_table),
      b("paragraph", paste(
        "Abbreviations: CI, confidence interval; HR, hazard ratio (each arm",
        "versus the control arm)."
      ))
    ))
  }
  blocks
}
