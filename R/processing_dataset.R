#' Process Dataset for Survival Analysis
#'
#' This function processes a dataset for survival analysis by standardizing
#' column names and creating necessary variables. It handles datasets with
#' either censoring indicators (CNSR) or event indicators (EVENT), or both.
#' The function supports both in-memory data objects (data.frame/tibble) and
#' file paths (CSV, XLSX, SAS formats).
#'
#' @param dataset A data frame/tibble containing survival analysis data, or a
#'   character string specifying the file path to a dataset. Supported file
#'   formats: CSV (.csv), Excel (.xlsx, .xls), SAS (.sas7bdat). The dataset
#'   must contain at least treatment arm, survival time, and event/censoring
#'   information
#' @param column_arm Character string specifying the column name for treatment arm
#'   in the input dataset. Will be renamed to "ARM" in output
#' @param column_survtime Character string specifying the column name for survival
#'   time in the input dataset. Will be renamed to "SURVTIME" in output
#' @param column_cnsr Character string specifying the column name for censoring
#'   indicator (1 = censored, 0 = event) in the input dataset. Default is NULL.
#'   If specified without column_event, EVENT will be calculated as 1 - CNSR
#' @param column_event Character string specifying the column name for event
#'   indicator (1 = event, 0 = censored) in the input dataset. Default is NULL.
#'   If both column_event and column_cnsr are specified, column_event takes priority
#'
#' @return A tibble containing standardized survival analysis data with columns:
#'   \item{ARM}{Treatment arm (factor)}
#'   \item{SURVTIME}{Survival time (numeric)}
#'   \item{EVENT}{Event indicator: 1 = event, 0 = censored (numeric)}
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item If dataset is a file path, reads the file based on extension
#'   \item Validates that at least one of column_cnsr or column_event is specified
#'   \item Renames ARM and SURVTIME columns to standardized names
#'   \item If only column_cnsr is specified, calculates EVENT = 1 - CNSR
#'   \item If only column_event is specified, uses it directly as EVENT
#'   \item If both are specified, column_event takes priority and column_cnsr is ignored
#' }
#'
#' File reading uses the following packages:
#' \itemize{
#'   \item CSV files: \code{readr::read_csv()}
#'   \item Excel files: \code{readxl::read_excel()}
#'   \item SAS files: \code{haven::read_sas()}
#' }
#'
#' @importFrom dplyr as_tibble rename mutate select
#' @importFrom rlang sym !!
#' @importFrom magrittr %>%
#' @importFrom readr read_csv
#' @importFrom readxl read_excel
#' @importFrom haven read_sas
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Using censoring indicator
#' dataset1 <- generate_dummy_survival_data(
#'   arm = c("Treatment", "Control"),
#'   n = c(100, 100),
#'   hazards = log(2) / c(20, 15),
#'   dropout_per_year = 0.1,
#'   seed = 123
#' )
#'
#' result1 <- processing_dataset(
#'   dataset = dataset1,
#'   column_arm = "ARM",
#'   column_survtime = "SURVTIME",
#'   column_cnsr = "CNSR",
#'   column_event = NULL
#' )
#'
#' head(result1)
#' # Output columns: ARM, SURVTIME, EVENT
#'
#' # Example 2: Using event indicator
#' dataset2 <- generate_dummy_survival_data(
#'   arm = c("Trt_A", "Trt_B"),
#'   n = c(75, 75),
#'   hazards = log(2) / c(24, 18),
#'   dropout_per_year = 0.15,
#'   seed = 456
#' )
#'
#' # Add EVENT column by converting CNSR
#' dataset2$EVENT <- 1 - dataset2$CNSR
#'
#' result2 <- processing_dataset(
#'   dataset = dataset2,
#'   column_arm = "ARM",
#'   column_survtime = "SURVTIME",
#'   column_cnsr = NULL,
#'   column_event = "EVENT"
#' )
#'
#' head(result2)
#' # Output columns: ARM, SURVTIME, EVENT
#'
#' # Example 3: Both censoring and event indicators (event takes priority)
#' dataset3 <- generate_dummy_survival_data(
#'   arm = c("Active", "Placebo"),
#'   n = c(80, 80),
#'   hazards = log(2) / c(22, 16),
#'   dropout_per_year = 0.12,
#'   seed = 789
#' )
#'
#' # Add EVENT column (will take priority over CNSR)
#' dataset3$EVENT <- 1 - dataset3$CNSR
#'
#' result3 <- processing_dataset(
#'   dataset = dataset3,
#'   column_arm = "ARM",
#'   column_survtime = "SURVTIME",
#'   column_cnsr = "CNSR",
#'   column_event = "EVENT"
#' )
#'
#' head(result3)
#' # Output columns: ARM, SURVTIME, EVENT
#'
#' # Verify that EVENT from column_event was used (not calculated from CNSR)
#' identical(result3$EVENT, dataset3$EVENT)
#' }
processing_dataset <- function(dataset,
                               column_arm,
                               column_survtime,
                               column_cnsr = NULL,
                               column_event = NULL) {

  # Validate inputs: at least one of column_cnsr or column_event must be specified
  if (is.null(column_cnsr) && is.null(column_event)) {
    stop("Either column_cnsr or column_event must be specified (not both NULL)")
  }

  # Check if dataset is a file path or data object
  if (is.character(dataset) && length(dataset) == 1) {
    # dataset is a file path
    if (!file.exists(dataset)) {
      stop("File not found: ", dataset)
    }

    # Get file extension
    file_ext <- tolower(tools::file_ext(dataset))

    # Read file based on extension
    dataset <- switch(
      file_ext,
      "csv" = {
        if (!requireNamespace("readr", quietly = TRUE)) {
          stop("Package 'readr' is required to read CSV files but not installed.")
        }
        readr::read_csv(dataset, show_col_types = FALSE)
      },
      "xlsx" = ,
      "xls" = {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          stop("Package 'readxl' is required to read Excel files but not installed.")
        }
        readxl::read_excel(dataset)
      },
      "sas7bdat" = {
        if (!requireNamespace("haven", quietly = TRUE)) {
          stop("Package 'haven' is required to read SAS files but not installed.")
        }
        haven::read_sas(dataset)
      },
      stop("Unsupported file format: ", file_ext, ". Supported formats: csv, xlsx, xls, sas7bdat")
    )
  } else if (!is.data.frame(dataset)) {
    # dataset is neither a file path nor a data.frame/tibble
    stop("'dataset' must be either a file path (character string) or a data.frame/tibble object")
  }

  # Start processing with renaming ARM and SURVTIME
  processed_data <- dataset %>%
    dplyr::as_tibble() %>%
    dplyr::rename(
      ARM = !!rlang::sym(column_arm),
      SURVTIME = !!rlang::sym(column_survtime)
    )

  # Handle EVENT column based on priority
  if (!is.null(column_event)) {
    # Priority 1: Use EVENT directly if specified
    processed_data <- processed_data %>%
      dplyr::rename(EVENT = !!rlang::sym(column_event))
  } else if (!is.null(column_cnsr)) {
    # Priority 2: Calculate EVENT from CNSR if EVENT is not specified
    processed_data <- processed_data %>%
      dplyr::rename(CNSR = !!rlang::sym(column_cnsr)) %>%
      dplyr::mutate(EVENT = 1 - CNSR) %>%
      dplyr::select(-CNSR)
  }

  # Select final columns
  processed_data <- processed_data %>%
    dplyr::select(ARM, SURVTIME, EVENT)

  return(processed_data)
}
