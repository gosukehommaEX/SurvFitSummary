# Small dataset shared by the tests. The same vectors are written in
# inst/validation/python/expected_km_cox.py, which computes the expected
# values independently with numpy and scipy.
make_test_raw <- function() {
  data.frame(
    USUBJID = sprintf("S%02d", 1:25),
    TRT01P = rep(c("Control", "Drug"), c(13, 12)),
    AVAL = c(2, 3, 3, 5, 6, 7, 7, 9, 11, 12, 14, 15, 18,
             3, 4, 6, 7, 8, 10, 12, 13, 16, 17, 20, 22),
    CNSR = 1 - c(1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0,
                 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0),
    AGE = c(61, 55, 70, 48, 66, 59, 73, 52, 64, 68, 57, 62, 71,
            58, 63, 49, 67, 54, 72, 60, 65, 51, 69, 56, 74),
    SEX = c(1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0,
            0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0),
    REGION = c("A", "B", "A", "A", "B", "B", "A", "B", "A", "B", "A", "B",
               "A", "B", "A", "B", "A", "B", "A", "A", "B", "A", "B", "A",
               "B"),
    stringsAsFactors = FALSE
  )
}

make_test_data <- function(control_arm = "Control") {
  processing_dataset(make_test_raw(), column_arm = "TRT01P",
                     column_survtime = "AVAL", column_cnsr = "CNSR",
                     column_covariates = c("AGE", "SEX"),
                     column_strata = "REGION", control_arm = control_arm)
}
