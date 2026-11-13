# Spline models of the package against rstpm2::stpm2().
#
# rstpm2 implements the Royston-Parmar models independently of flexsurv,
# with a natural spline basis (nsx) instead of the truncated power basis of
# flexsurv. With the same internal and boundary knots (log event times),
# both span the same restricted cubic splines of log time, so the maximized
# log-likelihood and the predicted survival must agree, although the spline
# coefficients differ. The scales correspond to link.type = "PH" (hazard),
# "PO" (odds) and "probit" (normal). For the arms fitted together, the arm
# shifts the spline (gamma0 in flexsurv, the coefficient of the arm in
# rstpm2), so the effect of the arm is compared too.
# Data: inst/extdata/fit_example.csv (all models), and flexsurv::bc (time
# in years, groups Good, Medium, Poor) for the models of all groups with 1
# knot on the hazard and odds scales (sp3 and sp1 of Jackson 2016, PH(2)
# and PO(2) of Royston and Parmar 2002), to check the estimates against
# which the published values are compared. Tolerance: 1e-3 for the
# log-likelihood and the effect of the arm, 1e-4 for the survival
# probabilities (both programs maximize numerically).
# rstpm2 is used only here, not by the package. Its fits are S4 objects
# (class "stpm2", based on "mle2"), whose methods are not found by
# stats::logLik(), stats::coef() and stats::predict(): the S4 generics
# are called directly (stats4::logLik(), stats4::coef() and the generic
# predict() of rstpm2).
#
# Run from the root of the package after devtools::load_all():
#   source(file.path("inst", "validation", "validate_splines_rstpm2.R"))
# Output: inst/validation/results/validate_splines_rstpm2.csv

source(file.path("inst", "validation", "verdict.R"))
script <- "validate_splines_rstpm2"
predict_stpm2 <- get("predict", envir = asNamespace("rstpm2"))
# stpm2() calls gsm() by evaluating its call in the environment of the
# caller, so the functions of rstpm2 must be found from here: rstpm2 is
# attached, and gsm is also bound here
suppressPackageStartupMessages(library(rstpm2))
gsm <- get("gsm", envir = asNamespace("rstpm2"))

dat <- processing_dataset(
  file.path("inst", "extdata", "fit_example.csv"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(dat, distributions = NULL,
                            spline_scales = c("hazard", "odds", "normal"),
                            spline_knots = 1:3)
gof <- gof_table(fits)
knots <- spline_knot_table(fits)
par <- parameter_table(fits, scale = "natural")
times <- c(3, 6, 12, 24, 36, 60, 120)
pred <- predict_survival(fits, times)
links <- c(hazard = "PH", odds = "PO", normal = "probit")
labels <- c(hazard = "Hazards", odds = "Odds", normal = "Normal")
arms <- levels(dat$ARM)

groups <- c(lapply(arms, function(a) {
  list(fit_type = "independent", arm = a, data = dat[dat$ARM == a, ],
       formula = survival::Surv(SURVTIME, EVENT) ~ 1)
}), list(list(fit_type = "dependent", arm = "Both arms", data = dat,
              formula = survival::Surv(SURVTIME, EVENT) ~ ARM)))

result <- NULL
for (g in groups) {
  boundary <- range(log(g$data$SURVTIME[g$data$EVENT == 1]))
  for (scale in names(links)) {
    for (n_knots in 1:3) {
      model <- paste0(labels[[scale]], ", ", n_knots,
                      if (n_knots == 1L) " knot" else " knots")
      internal <- knots$log_time[knots$fit_type == g$fit_type &
                                   knots$arm == g$arm &
                                   knots$n_knots == n_knots]
      smooth <- stats::as.formula(paste0(
        "~ rstpm2::nsx(log(SURVTIME), knots = c(",
        paste(sprintf("%.15g", internal), collapse = ", "),
        "), Boundary.knots = c(",
        paste(sprintf("%.15g", boundary), collapse = ", "), "))"
      ))
      ref <- rstpm2::stpm2(g$formula, data = g$data, smooth.formula = smooth,
                           link.type = links[[scale]])
      prefix <- paste0(g$arm, ", ", model, ": ")
      ll <- gof$loglik[gof$fit_type == g$fit_type & gof$arm == g$arm &
                         gof$model == model]
      result <- rbind(result, verdict_table(
        script, paste0(prefix, "log-likelihood"),
        as.numeric(stats4::logLik(ref)), ll, 1e-3
      ))
      for (a in if (g$fit_type == "dependent") arms else g$arm) {
        newdata <- data.frame(SURVTIME = times,
                              ARM = factor(a, levels = arms))
        s_ref <- as.numeric(predict_stpm2(ref, newdata = newdata,
                                          type = "surv"))
        rows <- pred$fit_type == g$fit_type & pred$arm == a &
          pred$model == model
        result <- rbind(result, verdict_table(
          script, paste0(prefix, a, " S(", times, ")"), s_ref,
          pred$est[rows][match(times, pred$time[rows])], 1e-4
        ))
      }
      if (g$fit_type == "dependent") {
        effect <- paste0("ARM", arms[2L])
        result <- rbind(result, verdict_table(
          script, paste0(prefix, "effect of ", arms[2L]),
          unname(stats4::coef(ref)[effect]),
          par$est[par$fit_type == "dependent" & par$model == model &
                    par$parameter == effect],
          1e-3
        ))
      }
    }
  }
}

# flexsurv::bc: models of all groups with 1 knot
bc <- processing_dataset(
  data.frame(GROUP = flexsurv::bc$group, YEARS = flexsurv::bc$recyrs,
             EVENT = flexsurv::bc$censrec),
  column_arm = "GROUP", column_survtime = "YEARS", column_event = "EVENT",
  control_arm = "Good"
)
bc_fits <- fit_survival_models(bc, distributions = NULL,
                               spline_scales = c("hazard", "odds"),
                               spline_knots = 1L, fit_type = "dependent",
                               time_unit = "years")
bc_gof <- gof_table(bc_fits)
bc_par <- parameter_table(bc_fits, scale = "natural")
bc_knots <- spline_knot_table(bc_fits)
bc_pred <- predict_survival(bc_fits, c(1, 2, 5))
boundary <- range(log(bc$SURVTIME[bc$EVENT == 1]))
smooth <- stats::as.formula(paste0(
  "~ rstpm2::nsx(log(SURVTIME), knots = ",
  sprintf("%.15g", bc_knots$log_time[bc_knots$n_knots == 1]),
  ", Boundary.knots = c(",
  paste(sprintf("%.15g", boundary), collapse = ", "), "))"
))
for (scale in c("hazard", "odds")) {
  model <- paste0(labels[[scale]], ", 1 knot")
  ref <- rstpm2::stpm2(survival::Surv(SURVTIME, EVENT) ~ ARM, data = bc,
                       smooth.formula = smooth, link.type = links[[scale]])
  prefix <- paste0("bc, Both arms, ", model, ": ")
  effects <- paste0("ARM", levels(bc$ARM)[-1L])
  rows <- bc_par$model == model
  result <- rbind(
    result,
    verdict_table(script, paste0(prefix, "log-likelihood"),
                  as.numeric(stats4::logLik(ref)),
                  bc_gof$loglik[bc_gof$model == model], 1e-3),
    verdict_table(script, paste0(prefix, "effect of ", effects),
                  unname(stats4::coef(ref)[effects]),
                  bc_par$est[rows][match(effects, bc_par$parameter[rows])],
                  1e-4)
  )
  for (a in levels(bc$ARM)) {
    newdata <- data.frame(SURVTIME = c(1, 2, 5),
                          ARM = factor(a, levels = levels(bc$ARM)))
    p_rows <- bc_pred$model == model & bc_pred$arm == a
    result <- rbind(result, verdict_table(
      script, paste0(prefix, a, " S(", c(1, 2, 5), ")"),
      as.numeric(predict_stpm2(ref, newdata = newdata, type = "surv")),
      bc_pred$est[p_rows][match(c(1, 2, 5), bc_pred$time[p_rows])], 1e-4
    ))
  }
}

write_verdicts(result, file.path("inst", "validation", "results"))
