# ============================================================
# INTERNAL MODEL HELPERS
# ============================================================
#
# This file contains internal utility functions shared by the
# package's public functions.
#
# The helpers provide a common interface for extracting:
#
#   - model class information;
#   - fixed-effect coefficients;
#   - fixed-effect variance-covariance matrices;
#   - model families and link functions;
#   - fixed-effect formulas;
#   - random-effect grouping variables;
#   - coefficient tables and p-values;
#   - fixed-effect model matrices;
#   - population-level predictions;
#   - weighted averages of model-based predictions.
#
# These functions are intentionally not exported. Their names begin
# with a dot to indicate that they are internal package utilities.
#
# Currently supported model classes include:
#
#   - stats::lm
#   - stats::glm
#   - lme4::lmer
#   - lme4::glmer
#   - lmerTest::lmer
#   - nlme::lme
#   - glmmTMB::glmmTMB
#
# For glmmTMB models, only the conditional component ("cond") is
# considered. Zero-inflation and dispersion components are not
# currently analyzed.
# ============================================================


# ============================================================
# 1. PACKAGE AVAILABILITY
# ============================================================

.require_optional_package <- function(package, reason = NULL) {

  if (!requireNamespace(package, quietly = TRUE)) {

    message_text <- paste0(
      "Package '",
      package,
      "' is required"
    )

    if (!is.null(reason)) {
      message_text <- paste0(
        message_text,
        " ",
        reason
      )
    }

    stop(
      paste0(message_text, "."),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# ============================================================
# 2. MODEL-CLASS IDENTIFICATION
# ============================================================

.is_supported_model <- function(model) {

  inherits(model, "lm") ||
    inherits(model, "glm") ||
    inherits(model, "merMod") ||
    inherits(model, "lme") ||
    inherits(model, "glmmTMB")
}


.is_mixed_model <- function(model) {

  inherits(model, "merMod") ||
    inherits(model, "lme") ||
    inherits(model, "glmmTMB")
}


.is_generalized_model <- function(model) {

  family_object <- .extract_family(model)

  if (is.null(family_object)) {
    return(FALSE)
  }

  is_gaussian_identity <- identical(
    family_object$family,
    "gaussian"
  ) &&
    identical(
      family_object$link,
      "identity"
    )

  !is_gaussian_identity
}


.model_class_label <- function(model) {

  paste(
    class(model),
    collapse = ", "
  )
}


# ============================================================
# 3. MODEL FRAME
# ============================================================

.extract_model_frame <- function(model) {

  # The factors of the model frame carry only the contrasts stored in the
  # data; those supplied through the 'contrasts' argument of the model
  # call are added from the fitted model.
  .apply_fitted_contrasts(model, .model_frame_data(model))
}


.model_frame_data <- function(model) {

  # For nlme::lme, stats::model.frame() returns the internal model
  # structure rather than the data. getData() returns the data used
  # for fitting, after removal of incomplete rows.
  if (inherits(model, "lme")) {

    .require_optional_package(
      "nlme",
      "to extract the data of this model"
    )

    model_data <- tryCatch(
      nlme::getData(model),
      error = function(e) NULL
    )

    if (is.data.frame(model_data)) {

      # getData() returns every column of the original data; keep only
      # the variables of the fixed and random formulas, as a model
      # frame would.
      model_variables <- unique(
        c(
          all.vars(stats::formula(model)),
          all.vars(model$call$random),
          all.vars(nlme::getGroupsFormula(model))
        )
      )

      model_variables <- intersect(model_variables, names(model_data))

      return(model_data[, model_variables, drop = FALSE])
    }
  }

  model_frame <- tryCatch(
    stats::model.frame(model),
    error = function(e) NULL
  )

  if (!is.null(model_frame)) {
    return(model_frame)
  }

  # Some model classes store the original data directly.
  if (!is.null(model$data) &&
      is.data.frame(model$data)) {

    return(model$data)
  }

  stop(
    paste0(
      "The model frame could not be extracted from this model. ",
      "Provide the original data explicitly to the calling function."
    ),
    call. = FALSE
  )
}


# ============================================================
# 4. MODEL FAMILY AND LINK
# ============================================================

.extract_family <- function(model) {

  family_object <- tryCatch(
    stats::family(model),
    error = function(e) NULL
  )

  if (!is.null(family_object)) {
    return(family_object)
  }

  family_object
}


.family_link_label <- function(model) {

  family_object <- .extract_family(model)

  if (is.null(family_object)) {
    return(NA_character_)
  }

  paste(
    family_object$family,
    family_object$link,
    sep = " / "
  )
}


# ============================================================
# 5. FIXED-EFFECT COEFFICIENTS
# ============================================================

.extract_fixef <- function(model) {

  if (inherits(model, "merMod")) {

    .require_optional_package(
      "lme4",
      "to extract fixed effects from this model"
    )

    return(
      lme4::fixef(model)
    )
  }

  if (inherits(model, "glmmTMB")) {

    .require_optional_package(
      "glmmTMB",
      "to extract fixed effects from this model"
    )

    fixed_effects <- glmmTMB::fixef(model)$cond

    if (is.null(fixed_effects)) {

      stop(
        paste0(
          "The conditional fixed-effect coefficients could not ",
          "be extracted from the glmmTMB model."
        ),
        call. = FALSE
      )
    }

    return(fixed_effects)
  }

  if (inherits(model, "lme")) {

    .require_optional_package(
      "nlme",
      "to extract fixed effects from this model"
    )

    return(
      nlme::fixef(model)
    )
  }

  fixed_effects <- stats::coef(model)

  if (!is.numeric(fixed_effects)) {

    stop(
      paste0(
        "Fixed-effect coefficients could not be extracted as a ",
        "numeric vector from model class: ",
        .model_class_label(model),
        "."
      ),
      call. = FALSE
    )
  }

  fixed_effects
}


# ============================================================
# 6. FIXED-EFFECT VARIANCE-COVARIANCE MATRIX
# ============================================================

.extract_vcov <- function(model) {

  variance_object <- stats::vcov(model)

  # glmmTMB may return separate variance-covariance matrices for
  # the conditional, zero-inflation, and dispersion components.
  if (inherits(model, "glmmTMB") &&
      is.list(variance_object) &&
      !is.null(variance_object$cond)) {

    variance_object <- variance_object$cond
  }

  variance_matrix <- as.matrix(variance_object)

  fixed_effects <- .extract_fixef(model)
  coefficient_names <- names(fixed_effects)

  # Reorder the matrix so that rows and columns correspond exactly
  # to the fixed-effect coefficient vector.
  if (!is.null(rownames(variance_matrix)) &&
      !is.null(colnames(variance_matrix)) &&
      all(coefficient_names %in% rownames(variance_matrix)) &&
      all(coefficient_names %in% colnames(variance_matrix))) {

    variance_matrix <- variance_matrix[
      coefficient_names,
      coefficient_names,
      drop = FALSE
    ]
  }

  expected_dimension <- length(fixed_effects)

  if (!identical(
    dim(variance_matrix),
    c(expected_dimension, expected_dimension)
  )) {

    stop(
      paste0(
        "The fixed-effect variance-covariance matrix has dimensions ",
        "that are incompatible with the fixed-effect coefficient vector."
      ),
      call. = FALSE
    )
  }

  variance_matrix
}


# ============================================================
# 7. FIXED-EFFECT FORMULA
# ============================================================

.extract_fixed_formula <- function(model) {

  model_formula <- stats::formula(model)

  # lme4 and glmmTMB include random-effect expressions directly
  # in the model formula.
  if (inherits(model, "merMod") ||
      inherits(model, "glmmTMB")) {

    if (requireNamespace("reformulas", quietly = TRUE)) {

      return(
        reformulas::nobars(model_formula)
      )
    }

    # Fallback for simple random-effect terms if reformulas is not
    # installed. This regular expression removes parenthesized
    # expressions containing "|" or "||".
    formula_string <- paste(
      deparse(model_formula),
      collapse = " "
    )

    formula_string <- gsub(
      "\\([^\\)]*\\|\\|?[^\\)]*\\)",
      "",
      formula_string
    )

    formula_string <- gsub(
      "\\+\\s*\\+",
      "+",
      formula_string
    )

    formula_string <- gsub(
      "\\+\\s*$",
      "",
      formula_string
    )

    return(
      stats::as.formula(formula_string)
    )
  }

  # nlme::lme stores the random-effects specification separately.
  # Therefore, formula(model) already describes the fixed component.
  model_formula
}


# ============================================================
# 8. RANDOM-EFFECT GROUPING VARIABLES
# ============================================================

.extract_random_effects <- function(model) {

  if (inherits(model, "merMod") ||
      inherits(model, "glmmTMB")) {

    model_formula <- stats::formula(model)

    if (requireNamespace("reformulas", quietly = TRUE)) {

      random_terms <- tryCatch(
        reformulas::findbars(model_formula),
        error = function(e) list()
      )

      if (length(random_terms) == 0L) {
        return(character(0))
      }

      grouping_variables <- vapply(
        random_terms,
        function(term) {

          paste(
            deparse(term[[3L]]),
            collapse = ""
          )
        },
        character(1)
      )

      return(
        unique(grouping_variables)
      )
    }

    # Fallback based on the formula string.
    formula_string <- paste(
      deparse(model_formula),
      collapse = " "
    )

    matches <- gregexpr(
      "\\([^\\)]*\\|\\|?[^\\)]*\\)",
      formula_string
    )[[1L]]

    if (identical(matches[1L], -1L)) {
      return(character(0))
    }

    random_expressions <- regmatches(
      formula_string,
      list(matches)
    )[[1L]]

    grouping_variables <- sub(
      ".*\\|\\|?",
      "",
      random_expressions
    )

    grouping_variables <- gsub(
      "\\)",
      "",
      grouping_variables
    )

    return(
      unique(trimws(grouping_variables))
    )
  }

  if (inherits(model, "lme")) {

    .require_optional_package(
      "nlme",
      "to extract random-effect grouping variables"
    )

    grouping_formula <- tryCatch(
      nlme::getGroupsFormula(
        model,
        asList = TRUE
      ),
      error = function(e) NULL
    )

    if (is.null(grouping_formula)) {
      return(character(0))
    }

    if (is.list(grouping_formula)) {

      grouping_variables <- unique(
        unlist(
          lapply(
            grouping_formula,
            all.vars
          ),
          use.names = FALSE
        )
      )

    } else {

      grouping_variables <- all.vars(
        grouping_formula
      )
    }

    return(grouping_variables)
  }

  character(0)
}


.random_effects_label <- function(model) {

  random_effects <- .extract_random_effects(model)

  if (length(random_effects) == 0L) {
    return(NA_character_)
  }

  paste(
    random_effects,
    collapse = ", "
  )
}


# ============================================================
# 9. COEFFICIENT TABLE
# ============================================================

.extract_coef_table <- function(model) {

  model_summary <- summary(model)

  if (inherits(model, "glmmTMB")) {

    coefficient_table <- model_summary$coefficients$cond

    return(
      list(
        table = coefficient_table,
        source = "model"
      )
    )
  }

  if (inherits(model, "lme")) {

    return(
      list(
        table = model_summary$tTable,
        source = "model"
      )
    )
  }

  # Models fitted directly with lmerTest already contain p-values.
  if (inherits(model, "lmerModLmerTest")) {

    return(
      list(
        table = model_summary$coefficients,
        source = "lmerTest"
      )
    )
  }

  # Standard lme4::lmer models do not report p-values. If lmerTest
  # is installed, convert the model internally to obtain approximate
  # denominator degrees of freedom and p-values.
  if (inherits(model, "lmerMod") &&
      requireNamespace("lmerTest", quietly = TRUE)) {

    converted_model <- lmerTest::as_lmerModLmerTest(
      model
    )

    return(
      list(
        table = summary(converted_model)$coefficients,
        source = "lmerTest"
      )
    )
  }

  coefficient_table <- model_summary$coefficients

  list(
    table = coefficient_table,
    source = "model"
  )
}


# ============================================================
# 10. P-VALUE EXTRACTION
# ============================================================

.extract_p_values <- function(model,
                              coefficient_names,
                              alpha = 0.05) {

  if (!is.character(coefficient_names)) {

    stop(
      "'coefficient_names' must be a character vector.",
      call. = FALSE
    )
  }

  coefficient_information <- tryCatch(
    .extract_coef_table(model),
    error = function(e) NULL
  )

  unavailable_result <- function() {

    list(
      table = data.frame(
        coefficient = coefficient_names,
        p_value = NA_real_,
        significance = "NoPValue",
        stringsAsFactors = FALSE
      ),
      source = "not_available"
    )
  }

  if (is.null(coefficient_information) ||
      is.null(coefficient_information$table)) {

    return(
      unavailable_result()
    )
  }

  coefficient_table <- as.matrix(
    coefficient_information$table
  )

  possible_p_columns <- c(
    "Pr(>|t|)",
    "Pr(>|z|)",
    "p-value"
  )

  p_column <- intersect(
    possible_p_columns,
    colnames(coefficient_table)
  )

  if (length(p_column) == 0L) {

    return(
      unavailable_result()
    )
  }

  all_p_values <- coefficient_table[ ,p_column[1L]]

  names(all_p_values) <- rownames(coefficient_table)

  requested_p_values <- all_p_values[coefficient_names]

  significance <- ifelse(
    is.na(requested_p_values),
    "NoPValue",
    ifelse(
      requested_p_values < alpha,
      "Sig",
      "NotSig"
    )
  )

  source <- if (all(is.na(requested_p_values))) {
    "not_available"
  } else {
    coefficient_information$source
  }

  list(
    table = data.frame(
      coefficient = coefficient_names,
      p_value = as.numeric(requested_p_values),
      significance = significance,
      stringsAsFactors = FALSE
    ),
    source = source
  )
}


#' Coefficient names associated with the main effect of one predictor
#'
#' Internal helper. Identifies the fixed-effect coefficients of a main-effect
#' term through the `assign` attribute of the fixed-effect model matrix, in
#' the order of the contrast columns. Matching by name prefix is avoided
#' because another variable may share the predictor's prefix (for example
#' `group` and `group_size`).
#'
#' @param model A fitted model supported by the package.
#' @param predictor Name of the predictor.
#'
#' @return A character vector of coefficient names.
#'
#' @keywords internal
#' @noRd
.focal_coefficient_names <- function(model,
                                     predictor) {

  terms_object <- stats::delete.response(
    stats::terms(
      .extract_fixed_formula(model)
    )
  )

  term_labels <- gsub(
    "`",
    "",
    attr(terms_object, "term.labels"),
    fixed = TRUE
  )

  term_index <- which(term_labels == predictor)

  if (length(term_index) != 1L) {

    stop(
      paste0(
        "The predictor '",
        predictor,
        "' does not appear as a main-effect term in the model."
      ),
      call. = FALSE
    )
  }

  # The model frame carries its own terms, so transformed variables
  # such as log(x) do not need to be re-evaluated.
  model_matrix <- stats::model.matrix(
    terms_object,
    data = .extract_model_frame(model)
  )

  coefficient_names <- colnames(model_matrix)[attr(model_matrix, "assign") == term_index]

  unmatched_names <- setdiff(
    coefficient_names,
    names(.extract_fixef(model))
  )

  if (length(unmatched_names) > 0L) {

    stop(
      paste0(
        "The following model-matrix columns of '",
        predictor,
        "' were not found among the fixed-effect coefficients: ",
        paste(unmatched_names, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  coefficient_names
}


# ============================================================
# 11. ORIGINAL CONTRASTS AND FACTOR LEVELS
# ============================================================

.extract_original_contrasts <- function(model) {

  fixed_formula <- .extract_fixed_formula(model)

  model_frame <- tryCatch(
    .extract_model_frame(model),
    error = function(e) NULL
  )

  if (is.null(model_frame)) {
    return(NULL)
  }

  fixed_variables <- all.vars(fixed_formula)

  fixed_variables <- fixed_variables[fixed_variables %in% names(model_frame)]

  factor_variables <- fixed_variables[
    vapply(
      model_frame[fixed_variables],
      is.factor,
      logical(1)
    )
  ]

  if (length(factor_variables) == 0L) {
    return(NULL)
  }

  contrast_list <- lapply(
    factor_variables,
    function(variable_name) {
      stats::contrasts(
        model_frame[[variable_name]]
      )
    }
  )

  names(contrast_list) <- factor_variables

  contrast_list <- contrast_list[
    !vapply(
      contrast_list,
      is.null,
      logical(1)
    )
  ]

  if (length(contrast_list) == 0L) {
    return(NULL)
  }

  contrast_list
}

.extract_xlevels <- function(model,
                             terms_object = NULL) {

  if (is.null(terms_object)) {

    terms_object <- stats::delete.response(
      stats::terms(
        .extract_fixed_formula(model)
      )
    )
  }

  original_model_frame <- tryCatch(
    .extract_model_frame(model),
    error = function(e) NULL
  )

  if (is.null(original_model_frame)) {
    return(NULL)
  }

  tryCatch(
    stats::.getXlevels(
      terms_object,
      original_model_frame
    ),
    error = function(e) NULL
  )
}




# ============================================================
# 12. FIXED-EFFECT MODEL MATRIX FOR NEW DATA
# ============================================================

.build_fixed_model_matrix <- function(model,
                                      newdata) {

  if (!is.data.frame(newdata)) {

    stop(
      "'newdata' must be a data frame.",
      call. = FALSE
    )
  }

  fixed_formula <- .extract_fixed_formula(model)

  terms_object <- stats::delete.response(stats::terms(fixed_formula))

  fixed_effects <- .extract_fixef(model)
  coefficient_names <- names(fixed_effects)

  original_contrasts <- .extract_original_contrasts(model)

  xlevels <- .extract_xlevels(
    model,
    terms_object = terms_object
  )

  # ----------------------------------------------------------
  # REMOVE CONTRAST ATTRIBUTES BEFORE model.frame()
  # ----------------------------------------------------------
  #
  # Subsets of the original data retain contrast attributes on
  # their factors. When model.frame() restores the complete set of
  # factor levels through xlev, it may remove those attributes and
  # issue repeated "contrasts dropped from factor" warnings.
  #
  # The contrasts are therefore removed before model.frame() and
  # restored explicitly afterwards.
  # ----------------------------------------------------------

  clean_newdata <- newdata

  for (variable_name in names(clean_newdata)) {

    if (is.factor(clean_newdata[[variable_name]])) {

      attr(
        clean_newdata[[variable_name]],
        "contrasts"
      ) <- NULL
    }
  }

  # Construct a model frame containing the complete factor-level
  # structure used when fitting the original model.
  new_model_frame <- stats::model.frame(
    terms_object,
    data = clean_newdata,
    xlev = xlevels,
    drop.unused.levels = FALSE,
    na.action = stats::na.pass
  )

  # ----------------------------------------------------------
  # RESTORE THE ORIGINAL CONTRAST MATRICES
  # ----------------------------------------------------------

  if (!is.null(original_contrasts)) {

    contrast_variables <- intersect(
      names(original_contrasts),
      names(new_model_frame)
    )

    for (variable_name in contrast_variables) {

      variable <- new_model_frame[[variable_name]]
      contrast_matrix <- original_contrasts[[variable_name]]

      if (!is.factor(variable) ||
          is.null(contrast_matrix)) {

        next
      }

      # Ensure that the rows of the contrast matrix follow the
      # factor-level order used in the new model frame.
      if (!is.null(rownames(contrast_matrix))) {

        missing_levels <- setdiff(
          levels(variable),
          rownames(contrast_matrix)
        )

        if (length(missing_levels) > 0L) {

          stop(
            paste0(
              "The stored contrast matrix for variable '",
              variable_name,
              "' does not contain the following levels: ",
              paste(
                missing_levels,
                collapse = ", "
              ),
              "."
            ),
            call. = FALSE
          )
        }

        contrast_matrix <- contrast_matrix[
          levels(variable),
          ,
          drop = FALSE
        ]
      }

      stats::contrasts(
        new_model_frame[[variable_name]]
      ) <- contrast_matrix
    }
  }

  # Contrasts are already stored on the factors in the reconstructed
  # model frame, so contrasts.arg is intentionally not used.
  new_model_matrix <- stats::model.matrix(
    terms_object,
    data = new_model_frame
  )

  # ----------------------------------------------------------
  # ALIGN THE NEW MATRIX WITH THE FITTED COEFFICIENT VECTOR
  # ----------------------------------------------------------

  missing_columns <- setdiff(
    coefficient_names,
    colnames(new_model_matrix)
  )

  # Missing columns may occur when newdata represents only a subset
  # of the possible factor-level configurations. Such columns are
  # added as zeros to preserve conformity with the fitted model.
  if (length(missing_columns) > 0L) {

    additional_columns <- matrix(
      0,
      nrow = nrow(new_model_matrix),
      ncol = length(missing_columns)
    )

    colnames(additional_columns) <- missing_columns

    new_model_matrix <- cbind(
      new_model_matrix,
      additional_columns
    )
  }

  extra_columns <- setdiff(
    colnames(new_model_matrix),
    coefficient_names
  )

  if (length(extra_columns) > 0L) {

    new_model_matrix <- new_model_matrix[
      ,
      setdiff(
        colnames(new_model_matrix),
        extra_columns
      ),
      drop = FALSE
    ]
  }

  remaining_missing_columns <- setdiff(
    coefficient_names,
    colnames(new_model_matrix)
  )

  if (length(remaining_missing_columns) > 0L) {

    stop(
      paste0(
        "The following fixed-effect columns could not be constructed: ",
        paste(
          remaining_missing_columns,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  new_model_matrix[ ,coefficient_names,drop = FALSE]
}
# ============================================================
# 13. POPULATION-LEVEL PREDICTIONS
# ============================================================

.compute_fixed_predictions <- function(model,
                                       newdata,
                                       prediction_scale = c(
                                         "response",
                                         "link"
                                       ),
                                       conf_level = 0.95) {

  prediction_scale <- match.arg(prediction_scale)

  if (!is.numeric(conf_level) ||
      length(conf_level) != 1L ||
      is.na(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1) {

    stop(
      "'conf_level' must be a single numeric value between 0 and 1.",
      call. = FALSE
    )
  }

  fixed_effects <- .extract_fixef(model)
  variance_matrix <- .extract_vcov(model)

  model_matrix <- .build_fixed_model_matrix(
    model = model,
    newdata = newdata
  )

  linear_predictor <- as.numeric(model_matrix %*% fixed_effects)

  # Only the diagonal of X V X' is needed. Computing it row by row
  # avoids building an n x n matrix.
  linear_predictor_variance <- rowSums((model_matrix %*% variance_matrix) * model_matrix)

  linear_predictor_se <- sqrt(
    pmax(
      linear_predictor_variance,
      0
    )
  )

  critical_value <- stats::qnorm(1 - (1 - conf_level) / 2)

  ci_low_link <- linear_predictor - critical_value * linear_predictor_se

  ci_high_link <- linear_predictor + critical_value * linear_predictor_se

  family_object <- .extract_family(model)

  if (prediction_scale == "link" ||
      is.null(family_object)) {

    predictions <- data.frame(
      estimate = linear_predictor,
      se = linear_predictor_se,
      ci_low = ci_low_link,
      ci_high = ci_high_link,
      stringsAsFactors = FALSE
    )

  } else {

    inverse_link <- family_object$linkinv
    inverse_link_derivative <- family_object$mu.eta

    response_estimate <- inverse_link(linear_predictor)

    response_se <- abs(inverse_link_derivative(linear_predictor)) * linear_predictor_se

    # Decreasing links (for example the inverse link) reverse the order
    # of the transformed bounds.
    transformed_low <- inverse_link(ci_low_link)
    transformed_high <- inverse_link(ci_high_link)

    predictions <- data.frame(
      estimate = response_estimate,
      se = response_se,
      ci_low = pmin(transformed_low, transformed_high),
      ci_high = pmax(transformed_low, transformed_high),
      stringsAsFactors = FALSE
    )
  }

  # Store intermediate quantities needed for variance propagation
  # and delta-method calculations of averages or other combinations.
  attr(predictions, "model_matrix") <- model_matrix
  attr(predictions, "vcov_beta") <- variance_matrix
  attr(predictions, "eta") <- linear_predictor
  attr(predictions, "family") <- family_object

  predictions
}


# ============================================================
# 14. WEIGHTED AVERAGE OF MODEL-BASED PREDICTIONS
# ============================================================

.compute_average_prediction <- function(model,
                                        newdata,
                                        prediction_scale = c(
                                          "response",
                                          "link"
                                        ),
                                        conf_level = 0.95,
                                        weights = NULL) {

  prediction_scale <- match.arg(prediction_scale)

  predictions <- .compute_fixed_predictions(
    model = model,
    newdata = newdata,
    prediction_scale = prediction_scale,
    conf_level = conf_level
  )

  .combine_fixed_predictions(
    predictions = predictions,
    weights = weights,
    prediction_scale = prediction_scale,
    conf_level = conf_level
  )
}


# ============================================================
# 15. WEIGHTED COMBINATION OF EXISTING PREDICTIONS
# ============================================================

#' Combine population-level predictions with arbitrary weights
#'
#' Internal helper. Operates on the object returned by
#' `.compute_fixed_predictions()`, so that several weighted averages of the
#' same predictions can be obtained without rebuilding the model matrix.
#' The standard error uses the complete fixed-effect covariance matrix, so
#' correlations between the combined predictions are taken into account.
#'
#' @param predictions Object returned by `.compute_fixed_predictions()`.
#' @param weights Non-negative weights, one per prediction row. `NULL` gives
#'   equal weights. Weights are normalized to sum to one.
#' @param prediction_scale `"response"` or `"link"`.
#' @param conf_level Confidence level.
#'
#' @return A one-row data frame with `estimate`, `se`, `ci_low`, `ci_high`.
#'
#' @keywords internal
#' @noRd
.combine_fixed_predictions <- function(predictions,
                                       weights = NULL,
                                       prediction_scale = c(
                                         "response",
                                         "link"
                                       ),
                                       conf_level = 0.95) {

  prediction_scale <- match.arg(prediction_scale)

  model_matrix <- attr(
    predictions,
    "model_matrix"
  )

  variance_matrix <- attr(
    predictions,
    "vcov_beta"
  )

  linear_predictor <- attr(
    predictions,
    "eta"
  )

  family_object <- attr(
    predictions,
    "family"
  )

  if (is.null(model_matrix) ||
      is.null(variance_matrix) ||
      is.null(linear_predictor)) {

    stop(
      paste0(
        "'predictions' must be an object returned by ",
        ".compute_fixed_predictions()."
      ),
      call. = FALSE
    )
  }

  number_predictions <- nrow(
    model_matrix
  )

  if (number_predictions == 0L) {

    return(
      data.frame(
        estimate = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  if (is.null(weights)) {

    weights <- rep(1,number_predictions)
  }

  if (!is.numeric(weights) ||
      length(weights) != number_predictions) {

    stop(
      paste0(
        "'weights' must be a numeric vector containing one value ",
        "for each prediction."
      ),
      call. = FALSE
    )
  }

  if (any(is.na(weights)) ||
      any(weights < 0) ||
      sum(weights) <= 0) {

    stop(
      paste0(
        "'weights' must contain non-negative, non-missing values ",
        "and must have a positive sum."
      ),
      call. = FALSE
    )
  }

  normalized_weights <- weights / sum(weights)

  critical_value <- stats::qnorm(1 - (1 - conf_level) / 2)

  if (prediction_scale == "link" ||
      is.null(family_object)) {

    average_model_row <- matrix(colSums(model_matrix *normalized_weights),nrow = 1L)

    estimate <- sum(normalized_weights *linear_predictor)

    estimate_variance <- as.numeric(
      average_model_row %*%
        variance_matrix %*%
        t(average_model_row)
    )

  } else {

    inverse_link <- family_object$linkinv
    inverse_link_derivative <- family_object$mu.eta

    response_predictions <- inverse_link(linear_predictor)

    estimate <- sum(normalized_weights * response_predictions)

    # For
    #
    #   f(beta) = sum_i a_i g^{-1}(X_i beta),
    #
    # the gradient is:
    #
    #   sum_i a_i g^{-1}'(eta_i) X_i.
    gradient <- colSums(model_matrix * as.numeric(inverse_link_derivative(linear_predictor)) * normalized_weights)

    gradient <- matrix(gradient,nrow = 1L)

    estimate_variance <- as.numeric(gradient %*% variance_matrix %*% t(gradient))
  }

  estimate_se <- sqrt(max(estimate_variance,0))

  ci_low <- estimate - critical_value * estimate_se
  ci_high <- estimate + critical_value * estimate_se

  # On the response scale, a symmetric interval can leave the range of
  # the response (for example a negative probability). The interval is
  # therefore built on the link scale and transformed back:
  #
  #   g^{-1}( g(theta) +/- z * SE(theta) / g^{-1}'(g(theta)) ),
  #
  # where SE(theta) / g^{-1}'(g(theta)) is the delta-method standard
  # error of g(theta). For a single prediction this reproduces exactly
  # the interval computed on the linear predictor and transformed.
  if (prediction_scale == "response" &&
      !is.null(family_object)) {

    link_estimate <- suppressWarnings(family_object$linkfun(estimate))

    link_derivative <- if (is.finite(link_estimate)) {
      as.numeric(family_object$mu.eta(link_estimate))
    } else {
      NA_real_
    }

    if (is.finite(link_derivative) &&
        link_derivative != 0) {

      link_se <- estimate_se / abs(link_derivative)

      # Decreasing links (for example the inverse link) reverse the
      # order of the bounds.
      bounds <- sort(
        family_object$linkinv(
          link_estimate + c(-1, 1) * critical_value * link_se
        )
      )

      ci_low <- bounds[1L]
      ci_high <- bounds[2L]
    }
  }

  data.frame(
    estimate = estimate,
    se = estimate_se,
    ci_low = ci_low,
    ci_high = ci_high,
    stringsAsFactors = FALSE
  )
}

#' Generate contrast colours for MoPlot
#'
#' Internal helper used to assign a distinct colour to each displayed contrast.
#'
#' @param n Number of contrast colours required.
#'
#' @return A character vector containing `n` colours.
#'
#' @keywords internal
#' @noRd
.moplot_contrast_palette <- function(n) {

  if (
    length(n) != 1L ||
    is.na(n) ||
    !is.numeric(n) ||
    !is.finite(n) ||
    n < 1 ||
    n != as.integer(n)
  ) {
    stop(
      "`n` must be a single positive integer.",
      call. = FALSE
    )
  }

  n <- as.integer(n)

  wong <- c(
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#CC79A7",
    "#D55E00",
    "#56B4E9",
    "#F0E442"
  )

  if (n <= length(wong)) {
    return(
      wong[seq_len(n)]
    )
  }

  grDevices::hcl.colors(
    n = n,
    palette = "Dark 3"
  )
}

#' Extract analysis data from a fitted model
#'
#' Internal helper. The model frame contains transformed columns such as
#' `log(x)` rather than the original variables, which are needed to build
#' prediction data. The data supplied in the model call are therefore
#' recovered when available, restricted to the rows used for fitting, and
#' given the contrasts used for fitting. The model frame is returned when
#' the original data cannot be recovered.
#'
#' @param model A fitted model object.
#'
#' @return A data frame, or `NULL`.
#'
#' @keywords internal
#' @noRd
.extract_analysis_data <- function(model) {

  model_frame <- tryCatch(
    .extract_model_frame(model),
    error = function(e) NULL
  )

  # nlme::getData() already returns the original variables.
  if (inherits(model, "lme")) {
    return(model_frame)
  }

  call_data <- tryCatch(
    eval(
      stats::getCall(model)$data,
      environment(stats::formula(model))
    ),
    error = function(e) NULL
  )

  if (!is.data.frame(call_data) ||
      is.null(model_frame) ||
      !all(all.vars(stats::formula(model)) %in% names(call_data))) {

    return(if (is.data.frame(model_frame)) model_frame else NULL)
  }

  .align_data_with_model(model, call_data, model_frame)
}


#' Align data with a fitted model
#'
#' Internal helper. Restricts `data` to the rows used for fitting and gives
#' its factors the contrasts used for fitting. See
#' `.restrict_to_model_rows()` and `.apply_fitted_contrasts()`.
#'
#' @param model A fitted model object.
#' @param data A data frame with the original variables.
#' @param model_frame The model frame of `model`, if already extracted.
#'
#' @return A data frame.
#'
#' @keywords internal
#' @noRd
.align_data_with_model <- function(model,
                                   data,
                                   model_frame = NULL) {

  .apply_fitted_contrasts(
    model,
    .restrict_to_model_rows(model, data, model_frame)
  )
}


#' Contrasts used for fitting
#'
#' Internal helper. Contrasts may come from the factors in the data or from
#' the `contrasts` argument of the model call; in both cases the model
#' records the coding actually used. Models of class `lme` store it in
#' `model$contrasts`, the other supported classes as the `"contrasts"`
#' attribute of their fixed-effect model matrix.
#'
#' @param model A fitted model object.
#'
#' @return A named list with one element per factor (a matrix, or the name
#'   of a contrast function), or `NULL`.
#'
#' @keywords internal
#' @noRd
.fitted_contrasts <- function(model) {

  if (inherits(model, "lme")) {
    return(model$contrasts)
  }

  model_matrix <- tryCatch(
    stats::model.matrix(model),
    error = function(e) NULL
  )

  attr(model_matrix, "contrasts")
}


#' Give factors the contrasts used for fitting
#'
#' Internal helper. Assigns to each factor of `data` the contrast matrix
#' with which the model was fitted, so that the hypothesis matrices match
#' the coefficients also when the coding was given through the `contrasts`
#' argument of the model call, or when the data were modified after
#' fitting. Levels without observations are dropped, as the fitting
#' functions do, unless the fitted contrast matrix includes them.
#'
#' @param model A fitted model object.
#' @param data A data frame with the original variables.
#'
#' @return A data frame.
#'
#' @keywords internal
#' @noRd
.apply_fitted_contrasts <- function(model,
                                    data) {

  if (!is.data.frame(data)) {
    return(data)
  }

  fitted_contrasts <- .fitted_contrasts(model)

  for (variable_name in intersect(names(fitted_contrasts), names(data))) {

    if (!is.factor(data[[variable_name]])) {
      next
    }

    coding <- fitted_contrasts[[variable_name]]

    # Contrast functions are given by name, as in the model call.
    if (is.character(coding)) {
      coding <- tryCatch(
        get(coding, mode = "function"),
        error = function(e) NULL
      )
    }

    empty_levels <- table(data[[variable_name]]) == 0L

    if (any(empty_levels) &&
        !(is.matrix(coding) && nrow(coding) == length(empty_levels))) {

      data[[variable_name]] <- droplevels(data[[variable_name]])
    }

    factor_levels <- levels(data[[variable_name]])

    if (is.function(coding)) {
      coding <- coding(factor_levels)
    }

    if (is.matrix(coding) &&
        nrow(coding) == length(factor_levels)) {

      stats::contrasts(data[[variable_name]], how.many = ncol(coding)) <- coding
    }
  }

  data
}


#' Restrict data to the rows used for fitting
#'
#' Internal helper. Rows excluded from the fit, through `subset` or by the
#' removal of incomplete rows, are dropped from `data`, so that observed
#' estimates and reference grids describe the same observations as the
#' model. Rows are matched by row name; `data` is returned unchanged when it
#' has no more rows than the model frame or when the row names do not match.
#'
#' @param model A fitted model object.
#' @param data A data frame with the original variables.
#' @param model_frame The model frame of `model`, if already extracted.
#'
#' @return A data frame.
#'
#' @keywords internal
#' @noRd
.restrict_to_model_rows <- function(model,
                                    data,
                                    model_frame = NULL) {

  if (!is.data.frame(data)) {
    return(data)
  }

  if (is.null(model_frame)) {
    model_frame <- tryCatch(
      .extract_model_frame(model),
      error = function(e) NULL
    )
  }

  if (!is.data.frame(model_frame) ||
      nrow(data) <= nrow(model_frame)) {

    return(data)
  }

  model_rows <- rownames(model_frame)

  if (!all(model_rows %in% rownames(data))) {
    return(data)
  }

  data[model_rows, , drop = FALSE]
}


#' Numeric response values for observed estimates
#'
#' Internal helper. Evaluates the response of the model in `data`, so that a
#' transformed response such as `log(y)` yields the values on the scale of
#' the model, and converts the result to numbers: logical responses and
#' factors of binomial models become 0/1 with the convention of
#' [stats::glm()] (the first level is 0).
#'
#' @param model A fitted model.
#' @param data A data frame with the original variables.
#' @param outcome Name of the response variable, or the deparsed response
#'   expression.
#'
#' @return A numeric vector with one value per row of `data`.
#'
#' @keywords internal
#' @noRd
.response_values <- function(model,
                             data,
                             outcome) {

  fixed_formula <- .extract_fixed_formula(model)

  values <- NULL

  if (length(fixed_formula) == 3L) {

    response_expression <- fixed_formula[[2L]]

    if (outcome %in% all.vars(response_expression) ||
        identical(outcome, paste(deparse(response_expression), collapse = ""))) {

      values <- tryCatch(
        eval(response_expression, data, environment(fixed_formula)),
        error = function(e) NULL
      )
    }
  }

  if (is.null(values)) {
    values <- data[[outcome]]
  }

  if (is.null(values)) {

    stop(
      paste0("The outcome '", outcome, "' was not found in 'data'."),
      call. = FALSE
    )
  }

  if (is.matrix(values)) {

    stop(
      paste0(
        "Observed estimates are not available for a response with several ",
        "columns, such as cbind(successes, failures). Use ",
        "estimate_type = \"predicted\"."
      ),
      call. = FALSE
    )
  }

  if (is.logical(values)) {
    return(as.numeric(values))
  }

  if (is.factor(values)) {

    family_object <- .extract_family(model)

    if (is.null(family_object) ||
        !family_object$family %in% c("binomial", "quasibinomial")) {

      stop(
        "A factor outcome is supported only for binomial models.",
        call. = FALSE
      )
    }

    # As in glm(): the first level denotes failure, all others success.
    return(as.numeric(values != levels(values)[1L]))
  }

  if (!is.numeric(values)) {

    stop(
      "The dependent variable must be numeric, logical, or a factor of a binomial model.",
      call. = FALSE
    )
  }

  as.numeric(values)
}
