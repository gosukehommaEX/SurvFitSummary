# The report is first built as a list of elements (.report_blocks()), then
# written with officer (.render_report()). The structure is checked against
# the report template; the numbers in the tables against the data and the
# values computed independently by inst/validation/python/expected_fits.py
# (to within the rounding of the report); the Word files are read back
# with officer.

report_overview <- list(data_label = NULL, outcome = NULL, author = NULL,
                        version = NULL, report_date = as.Date("2026-11-12"))

report_plan <- local({
  blocks <- NULL
  function() {
    if (is.null(blocks)) {
      x <- get_fits("report")
      s <- .report_settings(x, NULL, NULL, NULL, NULL)
      blocks <<- .report_blocks(x, s, report_overview, c(FALSE, TRUE),
                                n_sim = 1000, seed = 1)
    }
    blocks
  }
})

block_field <- function(blocks, field) {
  vapply(blocks, function(b) {
    if (is.null(b[[field]])) NA_character_ else as.character(b[[field]])
  }, character(1))
}

test_that("reporting times of the template by default", {
  x <- get_fits("report")
  s <- .report_settings(x, NULL, NULL, NULL, NULL)
  expect_identical(s$interval, 12)
  expect_identical(s$landmarks$time, c(6, 12, 24, 60, 120))
  expect_identical(s$landmarks$label,
                   c("6 months", "1 year", "2 years", "5 years", "10 years"))
  expect_identical(s$time_horizon, 120)
  expect_equal(.landmark_times("days")$time,
               c(182.625, 365.25, 730.5, 1826.25, 3652.5))
  expect_equal(.landmark_times("weeks")$time[2], 365.25 / 7)
  expect_equal(.landmark_times("years")$time, c(0.5, 1, 2, 5, 10))
  s <- .report_settings(x, 6, c(48, 36), NULL, NULL)
  expect_identical(s$interval, 6)
  expect_identical(s$landmarks$label, c("36 months", "48 months"))
  expect_identical(s$time_horizon, 48)
})

test_that("the headings follow the report template", {
  blocks <- report_plan()
  type <- block_field(blocks, "type")
  level <- block_field(blocks, "level")
  text <- block_field(blocks, "text")
  heading <- function(l) type == "heading" & level == l
  expect_identical(text[heading(1)], c("Overview", "Methods", "Results"))
  expect_identical(
    text[heading(2)],
    c("Kaplan-Meier analysis", "Hazard Ratios", "Diagnostics",
      "Standard parametric models", "Spline-based extrapolation",
      "Software", "Data overview", "Diagnostics",
      "Standard parametric models", "Spline-based extrapolation")
  )
  index <- seq_along(blocks)
  results <- which(text == "Results")
  parametric <- which(text == "Standard parametric models" & index > results)
  splines <- which(text == "Spline-based extrapolation" & index > results)
  distributions <- c("Distribution: Exponential", "Distribution: Weibull")
  expect_identical(
    text[heading(3) & index > parametric & index < splines],
    c("Control", distributions, "Treatment", distributions, "Both arms",
      distributions)
  )
  expect_identical(text[heading(3) & index > splines],
                   c("Control", "Treatment", "Both arms"))
  expect_identical(
    unique(text[heading(4)]),
    c("Smoothed empirical hazard plots", "Kaplan-Meier curves")
  )
  expect_identical(text[type == "title"], "Survival Analysis Report")
  expect_identical(text[type == "subtitle"], "Generated on: 2026-11-12")
  expect_true(any(grepl("^Data: \\[user to complete", text)))
})

test_that("figures with and without confidence intervals", {
  blocks <- report_plan()
  figures <- blocks[block_field(blocks, "type") == "figure"]
  # Kaplan-Meier curves and smoothed hazards: all arms and each arm, without
  # and with CI; distribution diagnostics, Schoenfeld residuals and QQ plot;
  # standard models: for each arm 2 overlaid plots and 4 plots per model,
  # for both arms 2 overlaid plots per arm and 4 plots per model; spline
  # models: 2 overlaid plots per arm, for each arm and for both arms
  n_models <- 2
  expected <- 2 * (1 + 2) * 2 + 3 + 2 * (2 + 4 * n_models) +
    (2 * 2 + 4 * n_models) + 2 * 2 + 2 * 2
  expect_equal(length(figures), expected)
  captions <- block_field(figures, "text")
  with_ci <- grepl("with 95% confidence intervals$", captions)
  ci_arg <- vapply(figures, function(b) isTRUE(b$plot$args$conf_int),
                   logical(1))
  expect_identical(with_ci, ci_arg)
  # Kaplan-Meier and hazards: all arms and each arm; standard models: 2
  # plots per model of each arm and of both arms
  expect_equal(sum(with_ci), (1 + 2) * 2 + 2 * n_models * 2 + n_models * 2)
  # The model plots of both arms show the models fitted to both arms
  both <- grepl("fitted to both arms", captions)
  fit_type <- vapply(figures, function(b) {
    if (is.null(b$plot$args$fit_type)) NA_character_ else
      b$plot$args$fit_type
  }, character(1))
  expect_true(all(fit_type[both] == "dependent"))
  expect_true(all(fit_type[!both & !is.na(fit_type)] == "independent"))
  overlay <- figures[[which(captions == paste(
    "Kaplan-Meier curve of Control and the survival curves of the standard",
    "parametric models"
  ))]]
  expect_identical(overlay$plot$fun, plot_survival_fits)
  expect_identical(overlay$plot$args$models, c("Exponential", "Weibull"))
  expect_identical(overlay$plot$args$time_horizon, 120)
})

test_that("data summary and Kaplan-Meier tables from the data", {
  blocks <- report_plan()
  tables <- blocks[block_field(blocks, "type") == "table"]
  captions <- block_field(tables, "text")
  d <- make_fit_data()
  summary <- tables[[which(captions == "Data summary")]]$table
  expect_identical(summary$Treatment, c("Control", "Treatment"))
  expect_identical(summary$`N patients`, c("300", "300"))
  expect_identical(summary$`N events`,
                   as.character(as.vector(tapply(d$EVENT, d$ARM, sum))))
  km <- tables[[which(captions == "Kaplan-Meier estimates of Control")]]
  km <- km$table
  times <- seq(0, 60, by = 12)
  expect_identical(km$`Time (months)`, as.character(times))
  control <- d$SURVTIME[d$ARM == "Control"]
  expect_identical(km$`Number at risk`,
                   as.character(vapply(times, function(t) {
                     sum(control >= t)
                   }, integer(1))))
  e <- read_expected("predictions")
  e <- e[e$fit_type == "independent" & e$arm == "Control" &
           e$model == "Exponential" & e$time %in% times, ]
  shown <- as.numeric(sub(" .*", "", km$`Survival probability (SE)`))
  expect_lte(max(abs(shown[match(e$time, times)] - e$km)), 5e-4 + 1e-10)
})

test_that("fit statistics and predictions agree with the reference", {
  blocks <- report_plan()
  tables <- blocks[block_field(blocks, "type") == "table"]
  captions <- block_field(tables, "text")
  gof <- tables[[which(captions == paste(
    "Goodness of fit statistics of the standard parametric models",
    "(Control), in ascending order of AIC"
  ))]]$table
  expect_false(is.unsorted(as.numeric(gof$AIC)))
  e <- read_expected("gof")
  e <- e[match(paste("independent", "Control", gof$Distribution),
               paste(e$fit_type, e$arm, e$model)), ]
  expect_lte(max(abs(as.numeric(gof$AIC) - e$AIC)), 0.05 + TOL_LOGLIK)
  expect_lte(max(abs(as.numeric(gof$BIC) - e$BIC)), 0.05 + TOL_LOGLIK)

  landmarks <- c("6 months", "1 year", "2 years", "5 years", "10 years")
  pred <- read_expected("predictions")
  summary <- read_expected("summary")
  for (a in c("Control", "Treatment")) {
    row <- tables[[which(captions == paste0(
      "Tests of fit and predicted survival probability (Both arms, ", a, ")"
    ))]]$table
    expect_identical(names(row),
                     c("Model", "AIC", "BIC", "Median (months)",
                       "Mean (months)", paste("Survival at", landmarks)))
    expect_identical(row$Model, "Hazards, 1 knot")
    p <- pred[pred$fit_type == "dependent" & pred$arm == a &
                pred$model == "Hazards, 1 knot", ]
    expected <- p$surv[match(c(6, 12, 24, 60, 120), p$time)]
    shown <- as.numeric(unlist(row[paste("Survival at", landmarks)]))
    expect_lte(max(abs(shown - expected)), 5e-4 + TOL_SURV)
    s <- summary[summary$fit_type == "dependent" & summary$arm == a &
                   summary$model == "Hazards, 1 knot", ]
    expect_lte(abs(as.numeric(row$`Median (months)`) - s$median),
               0.05 + TOL_TIME * s$median)
    expect_lte(abs(as.numeric(row$`Mean (months)`) - s$mean),
               0.05 + TOL_TIME * s$mean)
  }
  knots <- tables[[which(captions == "Placement of knots (Control)")]]$table
  k <- read_expected("knots")
  k <- k[k$fit_type == "independent" & k$arm == "Control" &
           k$n_knots == 1, ]
  expect_identical(knots$n_risk, as.character(k$n_risk))
  expect_lte(abs(as.numeric(knots$surv_est) - k$surv_est), 5e-4 + TOL_SURV)
})

test_that("the Word file holds the text, tables and figures in order", {
  x <- get_fits("report")
  s <- .report_settings(x, NULL, NULL, NULL, NULL)
  blocks <- c(.report_front(x, s, report_overview),
              .report_data_results(x, s, FALSE))
  path <- withr::local_tempfile(fileext = ".docx")
  .render_report(blocks, path, dpi = 30)
  content <- docx_content(path)
  type <- block_field(blocks, "type")
  text <- block_field(blocks, "text")
  number <- ave(seq_along(type), type, FUN = seq_along)
  expected <- ifelse(type == "table", paste0("Table ", number, ": ", text),
                     ifelse(type == "figure",
                            paste0("Figure ", number, ": ", text), text))
  expect_identical(docx_paragraphs(content), expected)
  expected_tables <- lapply(blocks[type == "table"], block_cells)
  expect_identical(docx_table_cells(path), expected_tables)
  expect_length(docx_media(path), sum(type == "figure"))
  # Figures within the text width of the A4 page (6.3 inches, in EMU)
  extent <- xml2::xml_find_all(docx_xml(path), "//wp:extent", docx_ns)
  expect_length(extent, sum(type == "figure"))
  expect_true(all(as.numeric(xml2::xml_attr(extent, "cx")) <= 6.3 * 914400))
  styles <- content$style_name[content$text %in% c("Overview", "Results")]
  expect_identical(styles, c("heading 1", "heading 1"))
})

test_that("a complete report for one arm with two models", {
  d <- make_fit_data()
  d <- d[d$ARM == "Control", ]
  x <- fit_survival_models(d, distributions = c("exp", "weibull"),
                           spline_scales = NULL, fit_type = "independent")
  path <- withr::local_tempfile(fileext = ".docx")
  out <- write_survival_report(x, path, data_label = "Example trial",
                               outcome = "Overall survival",
                               author = "Statistics", version = "1.0",
                               conf_int = FALSE, dpi = 30)
  expect_identical(out, path)
  text <- docx_paragraphs(docx_content(path))
  expect_true(all(c("Data: Example trial", "Outcome: Overall survival",
                    "Prepared by: Statistics", "Version: 1.0",
                    "Distribution: Exponential") %in% text))
  expect_false(any(c("Hazard Ratios", "Proportional hazards assessment",
                     "Spline-based extrapolation") %in% text))
  expect_false(any(grepl("could not be drawn", text)))
  # Kaplan-Meier curve, smoothed hazard, diagnostics, 2 overlaid plots and
  # the hazard and survival plots of each model
  expect_length(docx_media(path), 9)
})

test_that("a figure that cannot be drawn is reported in the text", {
  blocks <- list(
    .report_block("paragraph", "Before"),
    .report_block("figure", "Broken figure",
                  plot = list(fun = function() stop("no data"),
                              args = list())),
    .report_block("paragraph", "After")
  )
  path <- withr::local_tempfile(fileext = ".docx")
  .render_report(blocks, path, dpi = 30)
  expect_identical(
    docx_paragraphs(docx_content(path)),
    c("Before",
      "[The figure 'Broken figure' could not be drawn: no data]", "After")
  )
  expect_length(docx_media(path), 0)
})

test_that("styles of the document and formatting of numbers", {
  styles <- data.frame(style_type = c("paragraph", "paragraph", "table"),
                       style_name = c("Normal", "Heading 2", "Table Grid"),
                       stringsAsFactors = FALSE)
  expect_identical(.docx_style(styles, "heading 2"), "Heading 2")
  expect_identical(.docx_style(styles, c("Caption", "Normal")), "Normal")
  expect_identical(.docx_style(styles, c("table_template", "Table Grid"),
                               type = "table"), "Table Grid")
  expect_null(.docx_style(styles, "Title"))
  expect_identical(.format_number(c(-0.0001, 1.25, 12.3456, NA), 3),
                   c("0.000", "1.250", "12.346", "NE"))
  expect_identical(.format_number(NA_real_, 1, na = "NR"), "NR")
  expect_identical(.format_p(c(0.0004, 0.56214, 1, NA)),
                   c("<0.001", "0.562", "1.000", "NE"))
  expect_identical(.and_list(c("a", "b", "c")), "a, b and c")
  expect_identical(.and_list("a"), "a")
})

test_that("the template of the package gives the styles of the report", {
  blocks <- list(
    .report_block("title", "Survival Analysis Report"),
    .report_block("subtitle", "Generated on: 2026-11-13"),
    .report_block("heading", "Overview", level = 1L),
    .report_block("heading", "Smoothed empirical hazard plots", level = 4L),
    .report_block("table", "Example",
                  table = data.frame(Model = "Weibull", AIC = "1.0"))
  )
  path <- withr::local_tempfile(fileext = ".docx")
  .render_report(blocks, path)
  content <- docx_content(path)
  style_of <- function(text) {
    content$style_name[content$content_type == "paragraph" &
                         content$text %in% text]
  }
  expect_identical(style_of("Survival Analysis Report"), "Title")
  expect_identical(style_of("Generated on: 2026-11-13"), "Subtitle")
  expect_identical(style_of("Overview"), "heading 1")
  expect_identical(style_of("Smoothed empirical hazard plots"), "heading 4")
  expect_identical(style_of("Table 1: Example"), "Table Caption")
  doc <- docx_xml(path)
  table_style <- xml2::xml_find_all(doc, "//w:tbl/w:tblPr/w:tblStyle",
                                    docx_ns)
  expect_identical(xml2::xml_attr(table_style, "w:val", docx_ns),
                   "ReportTable")
  cell_style <- xml2::xml_find_all(doc, "//w:tc/w:p/w:pPr/w:pStyle", docx_ns)
  expect_identical(unique(xml2::xml_attr(cell_style, "w:val", docx_ns)),
                   "TableText")
  # A4 page
  page <- xml2::xml_find_first(doc, "//w:sectPr/w:pgSz", docx_ns)
  expect_identical(xml2::xml_attr(page, "w:w", docx_ns), "11906")
  expect_identical(xml2::xml_attr(page, "w:h", docx_ns), "16838")
})

test_that("a document given as template provides the styles", {
  template <- withr::local_tempfile(fileext = ".docx")
  .render_report(list(), template)
  blocks <- list(.report_block("heading", "Overview", level = 1L),
                 .report_block("paragraph", "Text"))
  path <- withr::local_tempfile(fileext = ".docx")
  .render_report(blocks, path, template = template)
  expect_identical(docx_paragraphs(docx_content(path)), c("Overview", "Text"))
})

test_that("invalid inputs are rejected", {
  x <- get_fits("report")
  path <- tempfile(fileext = ".docx")
  expect_error(write_survival_report(list(), path), "fit_survival_models")
  expect_error(write_survival_report(x, tempfile(fileext = ".doc")), "'file'")
  expect_error(write_survival_report(x, path, outcome = c("a", "b")),
               "'outcome'")
  expect_error(write_survival_report(x, path, conf_int = c(TRUE, TRUE)),
               "'conf_int'")
  expect_error(write_survival_report(x, path, template = "none.docx"),
               "'template'")
  expect_error(write_survival_report(x, path, dpi = 0), "'dpi'")
  expect_error(write_survival_report(x, path, n_sim = 1), "'n_sim'")
  expect_error(write_survival_report(x, path, seed = 0.5), "'seed'")
  expect_error(write_survival_report(x, path, interval = 0), "'interval'")
  expect_false(file.exists(path))
})
