# Global variables to avoid R CMD check NOTEs
# This file declares all global variables used in non-standard evaluation (NSE)
# contexts throughout the package, particularly in dplyr and ggplot2 chains

utils::globalVariables(c(
  # General data frame columns
  ".",
  "ARM",
  "SURVTIME",
  "EVENT",
  "CNSR",

  # Dummy data generation
  "SUBJID",
  "EVENTTIME",
  "DROPTIME",

  # Statistical model outputs
  "term",
  "estimate",
  "std.error",
  "statistic",
  "p.value",
  "conf.low",
  "conf.high",

  # Survival metrics
  "est",
  "lcl",
  "ucl",
  "Median Survival Time",
  "Mean Survival Time",
  "Restricted Mean Survival Time",

  # Hazard ratio outputs
  "HR",
  "L95%",
  "U95%",

  # Covariance matrix
  "row_param",

  # Plotting variables - time-related
  "time",
  "Time",
  "time_value",
  "time_scale_label",
  "log_time",

  # Plotting variables - survival-related
  "surv",
  "hazard",
  "lower",
  "upper",
  "lower_ci",
  "upper_ci",

  # Plotting variables - diagnostics
  "Residuals",
  "Comparison",
  "Annotation",
  "Annot_Color",
  "x_pos",
  "y_pos",

  # Q-Q plot
  "Control_Quantile",
  "Treatment_Quantile",
  "x",
  "y",
  "xend",
  "yend",

  # PH assumption plots
  "arm",
  "cloglog",
  "probit",

  # Type indicator
  "type"
))

# Import stats functions that were flagged
#' @importFrom stats coef relevel time
NULL
