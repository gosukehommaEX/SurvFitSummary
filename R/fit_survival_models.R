#' Fit Standard Parametric and Spline Survival Models
#'
#' @description
#' Fits parametric survival models with \code{flexsurv::flexsurvreg()} and
#' Royston-Parmar spline models with \code{flexsurv::flexsurvspline()} to
#' each treatment arm separately (independent fits) and to all arms together
#' with the arm as a covariate (dependent fits, labelled \code{"Both arms"}).
#' A model that fails is kept in the list of models with its error message
#' and the other models are still fitted; warnings of flexsurv are recorded
#' in the same way instead of being printed. The returned object is used by
#' \code{\link{parameter_table}}, \code{\link{covariance_table}},
#' \code{\link{gof_table}}, \code{\link{predict_survival}},
#' \code{\link{survival_summary}} and \code{\link{spline_knot_table}}.
#'
#' @details
#' The distributions, with their parameters in the order and on the natural
#' scale of flexsurv, are:
#' \itemize{
#'   \item \code{"exp"} (Exponential): rate.
#'   \item \code{"weibull"} (Weibull, accelerated failure time form): shape,
#'     scale; \code{S(t) = exp(-(t / scale)^shape)}.
#'   \item \code{"weibullPH"} (Weibull, proportional hazards form): shape,
#'     scale; \code{S(t) = exp(-scale * t^shape)}.
#'   \item \code{"rayleigh"} (Rayleigh): scale; the Weibull with the shape
#'     fixed at 2, so that \code{S(t) = exp(-(t / scale)^2)}. The fixed
#'     shape is not counted as a parameter.
#'   \item \code{"lnorm"} (Lognormal): meanlog, sdlog.
#'   \item \code{"llogis"} (Log-logistic): shape, scale;
#'     \code{S(t) = 1 / (1 + (t / scale)^shape)}.
#'   \item \code{"gompertz"} (Gompertz): shape, rate; hazard
#'     \code{rate * exp(shape * t)}.
#'   \item \code{"gamma"} (Gamma): shape, rate.
#'   \item \code{"gengamma"} (Generalized Gamma, Prentice 1974): mu, sigma, Q.
#'   \item \code{"gengamma.orig"} (Generalized Gamma, Stacy 1962): shape,
#'     scale, k. Equivalent to \code{"gengamma"} with Q > 0 only.
#'   \item \code{"genf"} (Generalized F, Prentice 1975): mu, sigma, Q, P.
#'   \item \code{"genf.orig"} (Generalized F, original form): mu, sigma, s1,
#'     s2.
#' }
#' The default is the seven standard parametric models of NICE DSU
#' Technical Support Document 14. In dependent fits the arm is a covariate
#' of the location parameter on its transformed scale: the rate of the
#' exponential, Gompertz and gamma, the scale of the Weibull, Weibull (PH),
#' Rayleigh, log-logistic and original generalized gamma, meanlog of the
#' lognormal, and mu of the generalized gamma and generalized F.
#'
#' The spline models (Royston and Parmar 2002) model the log cumulative
#' hazard (\code{"hazard"}), the log cumulative odds (\code{"odds"}) or the
#' probit of the survival (\code{"normal"}) as a natural cubic spline of log
#' time with parameters \code{gamma0}, \code{gamma1}, ...; in dependent fits
#' the arm is a covariate of \code{gamma0}. The boundary knots are the
#' smallest and largest log event times of the data fitted (both arms
#' together for dependent fits). By default the internal knots are at
#' equally spaced quantiles of the log event times: the median for one
#' knot, the 1/3 and 2/3 quantiles for two knots and the quartiles for three
#' knots.
#'
#' The likelihood is maximised by \code{stats::optim()} (BFGS) with the
#' settings in \code{control}; the default tolerance is smaller than that of
#' \code{optim()} so that the estimates are accurate to well below their
#' standard errors.
#'
#' @param data Data frame returned by \code{\link{processing_dataset}}.
#' @param distributions Character vector of distribution codes (see
#'   Details), or \code{NULL} for no parametric model. The default is the
#'   seven standard models \code{c("exp", "weibull", "lnorm", "llogis",
#'   "gompertz", "gengamma", "gamma")}.
#' @param spline_scales Character vector with the scales of the spline
#'   models, a subset of \code{c("hazard", "odds", "normal")} (the default),
#'   or \code{NULL} for no spline model.
#' @param spline_knots Numbers of internal knots of the spline models, whole
#'   numbers from 1 to 3. The default is \code{1:3}. Ignored when
#'   \code{knots} is given.
#' @param knots \code{NULL} (the default) for knots at quantiles of the log
#'   event times, or a list of increasing vectors of internal knots on the
#'   time scale, for example \code{list(12, c(6, 18))}; each vector gives
#'   one spline model per scale, and the vectors must have different
#'   lengths. The knots must lie strictly between the smallest and largest
#'   event times of each fit.
#' @param fit_type Character vector with \code{"independent"} (each arm
#'   separately) and/or \code{"dependent"} (all arms with the arm as a
#'   covariate). The default is both. With a single arm only independent
#'   models are fitted, with a warning if dependent models were requested.
#' @param time_unit Time unit of \code{SURVTIME}: \code{"months"} (the
#'   default), \code{"days"}, \code{"weeks"} or \code{"years"}. It is stored
#'   in the result for the labels and time conversions of the outputs.
#' @param conf_level Confidence level, a single number in (0, 1). The
#'   default is 0.95.
#' @param control List of settings passed to \code{stats::optim()}. The
#'   default is \code{list(reltol = 1e-10, maxit = 1000)}.
#'
#' @return An object of class \code{survfitsummary}, a list with
#'   \describe{
#'     \item{models}{Data frame with one row per model: \code{model_id},
#'       \code{fit_type} (\code{"independent"} or \code{"dependent"}),
#'       \code{arm} (the arm, or \code{"Both arms"} for dependent fits),
#'       \code{family} (\code{"parametric"} or \code{"spline"}),
#'       \code{distribution} (code, \code{NA} for splines), \code{scale} and
#'       \code{n_knots} (\code{NA} for parametric models), \code{model}
#'       (print-ready label such as \code{"Log-logistic"} or
#'       \code{"Hazards, 2 knots"}), \code{n} and \code{events} (of the data
#'       fitted), \code{status} (\code{"ok"}, \code{"warning"} or
#'       \code{"error"}) and \code{message} (warnings or error, \code{NA}
#'       if none).}
#'     \item{fits}{List of flexsurvreg objects in the order of
#'       \code{models}, \code{NULL} for models that failed.}
#'     \item{knots}{Data frame of the internal knots of the spline models:
#'       \code{fit_type}, \code{arm}, \code{n_knots}, \code{knot} (index),
#'       \code{quantile} and \code{log_time}.}
#'     \item{boundary_knots}{Data frame with \code{fit_type}, \code{arm},
#'       \code{lower} and \code{upper} (log time).}
#'     \item{data}{The columns \code{ARM}, \code{SURVTIME} and \code{EVENT}
#'       of \code{data}.}
#'     \item{arms}{Treatment arms, the control arm first.}
#'     \item{time_unit, conf_level}{The arguments.}
#'   }
#'
#' @references
#' Jackson C (2016). flexsurv: A platform for parametric survival modeling
#' in R. Journal of Statistical Software 70(8):1-33.
#' \doi{10.18637/jss.v070.i08}
#'
#' Latimer NR (2011, updated 2013). NICE DSU Technical Support Document 14:
#' Survival analysis for economic evaluations alongside clinical trials -
#' extrapolation with patient-level data.
#'
#' Royston P, Parmar MKB (2002). Flexible parametric proportional-hazards
#' and proportional-odds models for censored survival data, with
#' application to prognostic modelling and estimation of treatment
#' effects. Statistics in Medicine 21(15):2175-2197.
#' \doi{10.1002/sim.1203}
#'
#' @examples
#' dat <- processing_dataset(
#'   system.file("extdata", "fit_example.csv", package = "SurvFitSummary"),
#'   column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
#'   control_arm = "Control"
#' )
#' fits <- fit_survival_models(dat, distributions = c("exp", "weibull"),
#'                             spline_scales = "hazard", spline_knots = 1)
#' fits
#' fits$models[, c("fit_type", "arm", "model", "status")]
#'
#' @export
fit_survival_models <- function(data,
                                distributions = c("exp", "weibull", "lnorm",
                                                  "llogis", "gompertz",
                                                  "gengamma", "gamma"),
                                spline_scales = c("hazard", "odds",
                                                  "normal"),
                                spline_knots = 1:3,
                                knots = NULL,
                                fit_type = c("independent", "dependent"),
                                time_unit = c("months", "days", "weeks",
                                              "years"),
                                conf_level = 0.95,
                                control = list(reltol = 1e-10,
                                               maxit = 1000)) {
  .check_survival_data(data)
  info <- .distribution_info()
  if (!is.null(distributions)) {
    if (!is.character(distributions) || length(distributions) == 0L ||
        !all(distributions %in% info$code)) {
      stop("'distributions' must be NULL or a subset of: ",
           paste(info$code, collapse = ", "), call. = FALSE)
    }
    distributions <- unique(distributions)
  }
  scale_labels <- c(hazard = "Hazards", odds = "Odds", normal = "Normal")
  if (!is.null(spline_scales)) {
    if (!is.character(spline_scales) || length(spline_scales) == 0L ||
        !all(spline_scales %in% names(scale_labels))) {
      stop("'spline_scales' must be NULL or a subset of: hazard, odds, ",
           "normal", call. = FALSE)
    }
    spline_scales <- unique(spline_scales)
  }
  if (is.null(distributions) && is.null(spline_scales)) {
    stop("No model to fit: 'distributions' and 'spline_scales' are both ",
         "NULL.", call. = FALSE)
  }
  if (!is.null(knots)) {
    valid <- is.list(knots) && length(knots) > 0L &&
      all(vapply(knots, function(k) {
        is.numeric(k) && length(k) > 0L && all(is.finite(k)) &&
          all(k > 0) && !is.unsorted(k, strictly = TRUE)
      }, logical(1)))
    if (!valid) {
      stop("'knots' must be NULL or a list of increasing vectors of ",
           "positive knot times.", call. = FALSE)
    }
    if (anyDuplicated(lengths(knots)) > 0L) {
      stop("The vectors of 'knots' must have different numbers of knots.",
           call. = FALSE)
    }
    spline_knots <- lengths(knots)
  } else if (!is.null(spline_scales)) {
    if (!is.numeric(spline_knots) || length(spline_knots) == 0L ||
        any(!is.finite(spline_knots)) ||
        any(spline_knots != round(spline_knots)) ||
        any(spline_knots < 1 | spline_knots > 3)) {
      stop("'spline_knots' must contain whole numbers from 1 to 3.",
           call. = FALSE)
    }
    spline_knots <- sort(unique(as.integer(spline_knots)))
  }
  if (!is.character(fit_type) || length(fit_type) == 0L ||
      !all(fit_type %in% c("independent", "dependent"))) {
    stop("'fit_type' must contain \"independent\" and/or \"dependent\".",
         call. = FALSE)
  }
  fit_type <- intersect(c("independent", "dependent"), fit_type)
  time_unit <- match.arg(time_unit)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.list(control)) {
    stop("'control' must be a list.", call. = FALSE)
  }

  data <- data[c("ARM", "SURVTIME", "EVENT")]
  data$ARM <- droplevels(data$ARM)
  arms <- levels(data$ARM)
  if ("dependent" %in% fit_type && length(arms) < 2L) {
    warning("Dependent models need two or more treatment arms; only ",
            "independent models are fitted.", call. = FALSE)
    fit_type <- "independent"
  }

  # Data sets fitted: each arm, then all arms together
  groups <- list()
  if ("independent" %in% fit_type) {
    for (a in arms) {
      groups[[length(groups) + 1L]] <- list(
        fit_type = "independent", arm = a,
        data = data[data$ARM == a, , drop = FALSE]
      )
    }
  }
  if ("dependent" %in% fit_type) {
    groups[[length(groups) + 1L]] <- list(fit_type = "dependent",
                                          arm = "Both arms", data = data)
  }

  # User knots must lie inside the range of the event times of every fit
  if (!is.null(knots) && !is.null(spline_scales)) {
    for (g in groups) {
      event_times <- g$data$SURVTIME[g$data$EVENT == 1L]
      if (length(unique(event_times)) < 2L) next
      if (any(unlist(knots) <= min(event_times) |
              unlist(knots) >= max(event_times))) {
        stop("The knots must lie strictly between the smallest and ",
             "largest event times of each fit; for ", g$arm, " these are ",
             format(min(event_times)), " and ", format(max(event_times)),
             ".", call. = FALSE)
      }
    }
  }

  model_rows <- list()
  fits <- list()
  knot_rows <- list()
  boundary_rows <- list()
  add_model <- function(g, family, distribution, scale, n_knots, label,
                        result) {
    k <- length(model_rows) + 1L
    model_rows[[k]] <<- data.frame(
      model_id = k, fit_type = g$fit_type, arm = g$arm, family = family,
      distribution = distribution, scale = scale, n_knots = n_knots,
      model = label, n = nrow(g$data), events = sum(g$data$EVENT),
      status = result$status, message = result$message,
      stringsAsFactors = FALSE
    )
    fits[k] <<- list(result$fit)
  }
  for (g in groups) {
    dependent <- g$fit_type == "dependent"
    for (d in distributions) {
      result <- .fit_one_model(g$data, dependent = dependent,
                               distribution = d, conf_level = conf_level,
                               control = control)
      add_model(g, "parametric", d, NA_character_, NA_integer_,
                info$label[info$code == d], result)
    }
    if (is.null(spline_scales)) next
    knot_sets <- lapply(seq_along(spline_knots), function(j) {
      tryCatch(
        .spline_knots(g$data$SURVTIME, g$data$EVENT,
                      n_knots = spline_knots[j],
                      knot_times = if (is.null(knots)) NULL else knots[[j]]),
        error = function(e) e
      )
    })
    for (j in seq_along(knot_sets)) {
      ks <- knot_sets[[j]]
      if (inherits(ks, "error")) next
      knot_rows[[length(knot_rows) + 1L]] <- data.frame(
        fit_type = g$fit_type, arm = g$arm, n_knots = spline_knots[j],
        knot = seq_along(ks$internal), quantile = ks$quantile,
        log_time = ks$internal, stringsAsFactors = FALSE
      )
      if (length(boundary_rows) == 0L ||
          !identical(boundary_rows[[length(boundary_rows)]]$arm, g$arm)) {
        boundary_rows[[length(boundary_rows) + 1L]] <- data.frame(
          fit_type = g$fit_type, arm = g$arm, lower = ks$boundary[1L],
          upper = ks$boundary[2L], stringsAsFactors = FALSE
        )
      }
    }
    for (s in spline_scales) {
      for (j in seq_along(spline_knots)) {
        k <- spline_knots[j]
        ks <- knot_sets[[j]]
        result <- if (inherits(ks, "error")) {
          list(fit = NULL, status = "error", message = conditionMessage(ks))
        } else {
          .fit_one_model(g$data, dependent = dependent,
                         distribution = NA_character_, scale = s,
                         knots = ks, conf_level = conf_level,
                         control = control)
        }
        add_model(g, "spline", NA_character_, s, k,
                  paste0(scale_labels[[s]], ", ", k,
                         if (k == 1L) " knot" else " knots"),
                  result)
      }
    }
  }

  empty_knots <- data.frame(fit_type = character(0), arm = character(0),
                            n_knots = integer(0), knot = integer(0),
                            quantile = numeric(0), log_time = numeric(0),
                            stringsAsFactors = FALSE)
  empty_boundary <- data.frame(fit_type = character(0),
                               arm = character(0), lower = numeric(0),
                               upper = numeric(0), stringsAsFactors = FALSE)
  models <- do.call(rbind, model_rows)
  rownames(models) <- NULL
  knot_table <- if (length(knot_rows) > 0L) {
    do.call(rbind, knot_rows)
  } else {
    empty_knots
  }
  rownames(knot_table) <- NULL
  boundary_table <- if (length(boundary_rows) > 0L) {
    do.call(rbind, boundary_rows)
  } else {
    empty_boundary
  }
  rownames(boundary_table) <- NULL
  structure(
    list(models = models, fits = fits, knots = knot_table,
         boundary_knots = boundary_table, data = data, arms = arms,
         time_unit = time_unit, conf_level = conf_level),
    class = "survfitsummary"
  )
}
