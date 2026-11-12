# The tables of the sheets are those of the functions of the package (each
# compared with inst/validation/python in its own test file); these tests
# check that each table is written to its sheet unchanged, and the layout
# of the parameter sheets (variances on the diagonal of the covariance
# columns, from the standard errors).

workbook_file <- local({
  path <- NULL
  function() {
    if (is.null(path)) {
      path <<- tempfile(fileext = ".xlsx")
      write_survival_workbook(get_fits("report"), path, conf_int = FALSE)
    }
    path
  }
})

test_that("each table of the package is written to its sheet", {
  x <- get_fits("report")
  path <- workbook_file()
  expect_identical(
    openxlsx::getSheetNames(path),
    c("Settings", "Data summary", "KM table", "Hazard ratios", "PH test",
      "Model status", "Fit statistics", "Parameters",
      "Parameters (transformed)", "Covariance", "Cholesky",
      "Survival summary", "Landmark survival", "Spline knots")
  )
  expect_same_sheet(read_sheet(path, "Data summary"), summarise_data(x$data))
  expect_same_sheet(read_sheet(path, "KM table"),
                    km_table(x$data, interval = 12))
  expect_same_sheet(read_sheet(path, "Hazard ratios"),
                    estimate_hazard_ratio(x$data))
  expect_same_sheet(read_sheet(path, "PH test"),
                    test_proportional_hazards(x$data))
  expect_same_sheet(read_sheet(path, "Fit statistics"), gof_table(x))
  expect_same_sheet(read_sheet(path, "Covariance"),
                    covariance_table(x, scale = c("transformed", "natural")))
  expect_same_sheet(read_sheet(path, "Cholesky"),
                    covariance_table(x, type = "cholesky",
                                     scale = c("transformed", "natural")))
  expect_same_sheet(read_sheet(path, "Spline knots"), spline_knot_table(x))
  status <- read_sheet(path, "Model status")
  expect_identical(status$model, x$models$model)
  expect_identical(status$status, x$models$status)
})

test_that("predictions at the reporting and landmark times", {
  x <- get_fits("report")
  path <- workbook_file()
  times <- c(0, 6, seq(12, 120, by = 12))
  expect_same_sheet(read_sheet(path, "Landmark survival"),
                    predict_survival(x, times))
  expect_same_sheet(read_sheet(path, "Survival summary"),
                    survival_summary(x, time_horizon = 120))
  settings <- read_sheet(path, "Settings")
  value <- function(item) settings$value[settings$item == item]
  expect_identical(value("Time unit"), "months")
  expect_identical(value("Reporting interval"), "12")
  expect_identical(value("Landmark times"), "6, 12, 24, 60, 120")
  expect_identical(value("Time horizon"), "120")
  expect_identical(value("Confidence intervals of the predictions"), "No")
})

test_that("parameters with the rows of their covariance matrix", {
  x <- get_fits("report")
  path <- workbook_file()
  for (scale in c("natural", "transformed")) {
    sheet <- if (scale == "natural") {
      "Parameters"
    } else {
      "Parameters (transformed)"
    }
    p <- read_sheet(path, sheet)
    par <- parameter_table(x, scale = scale)
    expect_identical(p$parameter, par$parameter)
    expect_equal(p$est, par$est, tolerance = 1e-12)
    expect_equal(p$se, par$se, tolerance = 1e-12)
    cov <- covariance_table(x, scale = scale)
    key <- paste(p$fit_type, p$arm, p$model)
    cov_key <- paste(cov$fit_type, cov$arm, cov$model)
    for (k in unique(key)) {
      rows <- which(key == k)
      for (i in seq_along(rows)) {
        for (j in seq_along(rows)) {
          expected <- cov$value[cov_key == k &
                                  cov$row == p$parameter[rows[i]] &
                                  cov$col == p$parameter[rows[j]]]
          expect_equal(p[[paste0("cov_", j)]][rows[i]], expected,
                       tolerance = 1e-12)
        }
        # Variance on the diagonal: the square of the standard error
        expect_equal(p[[paste0("cov_", i)]][rows[i]], p$se[rows[i]]^2,
                     tolerance = 1e-10)
      }
      # Columns beyond the parameters of the model are empty
      if (length(rows) < 3L) {
        expect_true(all(is.na(p$cov_3[rows])))
      }
    }
  }
})

test_that("settings given by the user, one arm and confidence intervals", {
  d <- make_fit_data()
  d <- d[d$ARM == "Control", ]
  x <- fit_survival_models(d, distributions = "exp", spline_scales = NULL,
                           fit_type = "independent")
  path <- withr::local_tempfile(fileext = ".xlsx")
  write_survival_workbook(x, path, interval = 6, landmarks = c(48, 36),
                          time_horizon = 60, n_sim = 20)
  expect_false(any(c("Hazard ratios", "PH test", "Spline knots") %in%
                     openxlsx::getSheetNames(path)))
  expect_same_sheet(read_sheet(path, "KM table"),
                    km_table(x$data, interval = 6))
  pred <- read_sheet(path, "Landmark survival")
  expect_equal(pred$time, seq(0, 60, by = 6))
  expect_true(all(pred$lcl[-1] < pred$est[-1] & pred$est[-1] < pred$ucl[-1]))
  summary <- read_sheet(path, "Survival summary")
  expect_false(is.na(summary$rmst_lcl))
})

test_that("hazard ratios given by the user are written as given", {
  x <- get_fits("report")
  hr <- estimate_hazard_ratio(x$data, ties = "breslow")
  path <- withr::local_tempfile(fileext = ".xlsx")
  write_survival_workbook(x, path, hazard_ratios = hr, conf_int = FALSE)
  expect_same_sheet(read_sheet(path, "Hazard ratios"), hr)
})

test_that("invalid inputs are rejected", {
  x <- get_fits("report")
  path <- tempfile(fileext = ".xlsx")
  expect_error(write_survival_workbook(list(), path), "fit_survival_models")
  expect_error(write_survival_workbook(x, tempfile(fileext = ".csv")),
               "'file'")
  expect_error(write_survival_workbook(x, file.path(tempfile(), "a.xlsx")),
               "folder")
  expect_error(write_survival_workbook(x, path, conf_int = NA), "'conf_int'")
  expect_error(write_survival_workbook(x, path, interval = -1), "'interval'")
  expect_error(write_survival_workbook(x, path, landmarks = 0), "'landmarks'")
  expect_error(write_survival_workbook(x, path, time_horizon = 0),
               "'time_horizon'")
  expect_error(write_survival_workbook(x, path, hazard_ratios = data.frame()),
               "'hazard_ratios'")
  expect_false(file.exists(path))
})
