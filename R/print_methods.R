# ============================================================
# PRINT METHODS FOR MODEL CONTRAST TABLE OBJECTS (INTERNAL)
# ============================================================

# ============================================================
# 1. COMPLETE MODEL CONTRAST TABLE
# ============================================================

#' Print a complete model contrast summary
#'
#' @param x An object of class `"model_contrast_table"`.
#' @param ... Additional arguments passed to [print()].
#'
#' @return The input object `x`, invisibly.
#'
#' @importFrom cli cli_h1 cli_h2
#' @importFrom tibble as_tibble
#' @exportS3Method
#' @noRd
print.model_contrast_table <- function(x, ...) {

  cli::cli_h1("Model Contrast Summary")

  cli::cli_h2("Model-level information")

  model_summary <- x[["Model-Level Summary"]]

  # Capture extra arguments passed in '...'
  extra_args <- list(...)

  # Set default printing parameters if not provided in '...'
  if (!"n" %in% names(extra_args)) extra_args$n <- Inf
  if (!"width" %in% names(extra_args)) extra_args$width <- Inf

  do.call(print, c(list(x = tibble::as_tibble(model_summary)), extra_args))

  cli::cli_h2("Variable-level information")

  variable_summary <- x[["Variable-Level Summary"]]

  do.call(print, c(list(x = tibble::as_tibble(variable_summary)), extra_args))

  invisible(x)
}


# ============================================================
# 2. MODEL-LEVEL SUMMARY
# ============================================================

#' Print a model-level contrast summary
#'
#' @param x An object of class `"model_contrast_summary"`.
#' @param ... Additional arguments passed to [print()].
#'
#' @return The input object `x`, invisibly.
#'
#' @importFrom cli cli_h1
#' @importFrom tibble as_tibble
#' @exportS3Method
#' @noRd
print.model_contrast_summary <- function(x, ...) {

  cli::cli_h1("Model-level summary")

  printable_object <- as.data.frame(x, stringsAsFactors = FALSE)

  extra_args <- list(...)
  if (!"n" %in% names(extra_args)) extra_args$n <- Inf
  if (!"width" %in% names(extra_args)) extra_args$width <- Inf

  do.call(print, c(list(x = tibble::as_tibble(printable_object)), extra_args))

  invisible(x)
}


# ============================================================
# 3. VARIABLE-LEVEL SUMMARY
# ============================================================

#' Print a variable-level contrast summary
#'
#' @param x An object of class `"model_contrast_variables"`.
#' @param ... Additional arguments passed to [print()].
#'
#' @return The input object `x`, invisibly.
#'
#' @importFrom cli cli_h1
#' @importFrom tibble as_tibble
#' @exportS3Method
#' @noRd
print.model_contrast_variables <- function(x, ...) {

  cli::cli_h1("Variable-level summary")

  printable_object <- as.data.frame(x, stringsAsFactors = FALSE)

  extra_args <- list(...)
  if (!"n" %in% names(extra_args)) extra_args$n <- Inf
  if (!"width" %in% names(extra_args)) extra_args$width <- Inf

  do.call(print, c(list(x = tibble::as_tibble(printable_object)), extra_args))

  invisible(x)
}
