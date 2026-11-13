# Excel workbook and Word report of inst/extdata/fit_example.csv, for the
# independent checks by python/check_workbook.py and python/check_report.py.
#
# All the models of fit_survival_models() with its defaults (the seven
# standard distributions and the spline models on the hazard, odds and
# normal scales with 1 to 3 knots), the default landmark times and the
# hazard ratio of the arm without adjustment. The checks are run with
# Python (numpy, scipy, openpyxl and python-docx; see README.md):
#   python3 -I inst/validation/python/check_workbook.py \
#     inst/validation/results/outputs/fit_example_results.xlsx \
#     inst/validation/python/expected inst/extdata/fit_example.csv \
#     --csv inst/validation/results/outputs/check_workbook.csv
# then the same with --mutate and --csv
# inst/validation/results/outputs/mutation/check_workbook.csv, and likewise
# check_report.py with fit_example_report.docx. The workbook and the report
# are not kept in git or in the package (.gitignore, .Rbuildignore): the
# report has more than 100 figures (about 16 MB). The verdicts are kept.
#
# Run from the root of the package after devtools::load_all():
#   source(file.path("inst", "validation", "validate_outputs.R"))
# It takes 2 to 3 minutes (the figures of the report).

dat <- processing_dataset(
  file.path("inst", "extdata", "fit_example.csv"),
  column_arm = "TRT01P", column_survtime = "AVAL", column_cnsr = "CNSR",
  control_arm = "Control"
)
fits <- fit_survival_models(dat)
folder <- file.path("inst", "validation", "results", "outputs")
dir.create(folder, showWarnings = FALSE, recursive = TRUE)
write_survival_workbook(fits, file.path(folder, "fit_example_results.xlsx"))
write_survival_report(fits, file.path(folder, "fit_example_report.docx"))
print(fits$models[fits$models$status != "ok",
                  c("fit_type", "arm", "model", "status", "message")])
