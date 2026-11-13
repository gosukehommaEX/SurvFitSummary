# Runs all the validation and reproduction scripts twice and summarizes the
# verdicts.
#
# 1. Normal run: the verdicts are written to inst/validation/results/ and
#    inst/reproduce/results/.
# 2. Mutation run (SURVFITSUMMARY_MUTATION = "1"): the values of the package
#    are shifted before each comparison (see verdict.R) and the verdicts are
#    written to results/mutation/. No comparison may PASS: otherwise the
#    comparison is vacuous.
# The summary (number of verdicts of each kind per script and run, and the
# check of the mutation run) is written to
# inst/validation/results/validation_summary.csv.
#
# Run from the root of the package after devtools::load_all():
#   source(file.path("inst", "validation", "run_validation.R"))
# It takes a few minutes.

scripts <- c(
  file.path("inst", "validation", "validate_survreg.R"),
  file.path("inst", "validation", "validate_closed_forms.R"),
  file.path("inst", "validation", "validate_splines_rstpm2.R"),
  file.path("inst", "reproduce", "reproduce_Jackson2016.R"),
  file.path("inst", "reproduce", "reproduce_RoystonParmar2002.R")
)
summary <- NULL
for (run in c("normal", "mutation")) {
  Sys.setenv(SURVFITSUMMARY_MUTATION = if (run == "mutation") "1" else "0")
  for (path in scripts) {
    started <- Sys.time()
    source(path, local = new.env())
    name <- sub("\\.R$", "", basename(path))
    folder <- file.path(dirname(path), "results")
    if (run == "mutation") folder <- file.path(folder, "mutation")
    verdicts <- utils::read.csv(file.path(folder, paste0(name, ".csv")),
                                stringsAsFactors = FALSE)
    counts <- table(factor(verdicts$verdict,
                           levels = c("PASS", "EXPLAINED", "FAIL", "INFO")))
    summary <- rbind(summary, data.frame(
      run = run, script = name, PASS = counts[["PASS"]],
      EXPLAINED = counts[["EXPLAINED"]], FAIL = counts[["FAIL"]],
      INFO = counts[["INFO"]],
      seconds = round(as.numeric(difftime(Sys.time(), started,
                                          units = "secs"))),
      check = if (run == "normal") {
        if (counts[["FAIL"]] == 0L) "OK" else "FAIL present"
      } else {
        if (counts[["PASS"]] == 0L) "OK" else "PASS in mutation run"
      },
      stringsAsFactors = FALSE
    ))
  }
}
Sys.unsetenv("SURVFITSUMMARY_MUTATION")
utils::write.csv(summary, file.path("inst", "validation", "results",
                                    "validation_summary.csv"),
                 row.names = FALSE)
print(summary)
