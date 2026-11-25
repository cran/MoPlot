# ============================================================
# REFERENCE GRID FOR MODEL-BASED PREDICTIONS
# ============================================================
#
# Model-based estimates are computed on a reference grid rather
# than on the observed rows. The grid contains every combination of
# the levels of the focal factors and of the other factors in the
# fixed-effect formula; numeric predictors are fixed at their mean.
#
# Within each focal cell, the rows of the grid are weighted:
#
#   - factors that do not interact with the focal factors receive
#     equal weights (a balanced average over their levels);
#   - factors that share an interaction term with a focal factor
#     receive the weights of the intercept row of their hypothesis
#     matrix. The focal coefficients are defined at exactly this
#     combination of their levels: the balanced average for sum,
#     Helmert, and sliding-difference coding, the reference level for
#     treatment coding.
#
# On the linear-predictor scale, the displayed comparisons then
# reconstruct the focal coefficients exactly. Averaging observed rows
# instead would mix in the effects of covariates whose distribution
# differs between cells.
# ============================================================


#' Variables on the right-hand side of the fixed-effect formula
#'
#' @keywords internal
#' @noRd
.fixed_predictor_variables <- function(model) {

  fixed_formula <- .extract_fixed_formula(model)

  right_hand_side <- fixed_formula[[length(fixed_formula)]]

  all.vars(right_hand_side)
}


#' Name of the dependent variable of a model
#'
#' The first variable of the left-hand side of the fixed-effect formula,
#' such as `y` for `log(y) ~ ...`.
#'
#' @param model A fitted model.
#'
#' @return A character string.
#'
#' @keywords internal
#' @noRd
.outcome_name <- function(model) {

  fixed_formula <- .extract_fixed_formula(model)

  if (length(fixed_formula) < 3L) {
    stop("The model formula has no dependent variable.", call. = FALSE)
  }

  all.vars(fixed_formula[[2L]])[1L]
}


#' Resolve `estimate_type = "auto"`
#'
#' Internal helper shared by `contrasts_onepred()`, `contrasts_multpred()`,
#' and the plotting functions. Observed means are used only when they
#' answer the same question as the model: the fixed part contains the focal
#' predictors alone, the model has no random effects (observed means would
#' ignore the clustering), the model has no unequal prior weights (observed
#' means are unweighted), and, for
#' generalized models, estimates are requested on the response scale
#' (observed means are not on the link scale). Otherwise model-based
#' predictions are used.
#'
#' @param model A fitted model.
#' @param predictor Names of the focal predictors (several for an interaction).
#' @param estimate_type `"auto"`, `"observed"`, or `"predicted"`.
#' @param prediction_scale `"response"` or `"link"`.
#'
#' @return `"observed"` or `"predicted"`.
#'
#' @keywords internal
#' @noRd
.resolve_estimate_type <- function(model,
                                   predictor,
                                   estimate_type,
                                   prediction_scale) {

  if (!identical(estimate_type, "auto")) {
    return(estimate_type)
  }

  # 'predictor' may contain several variables, such as those of an
  # interaction: the fixed part must contain exactly these variables.
  focal_predictor_only <- setequal(
    .fixed_predictor_variables(model),
    predictor
  )

  link_scale_requested <- .is_generalized_model(model) &&
    identical(prediction_scale, "link")

  if (focal_predictor_only &&
      !.is_mixed_model(model) &&
      !.has_unequal_weights(model) &&
      !link_scale_requested) {

    return("observed")
  }

  "predicted"
}


#' Detect unequal prior weights
#'
#' Internal helper. Prior weights change the estimates only when they are
#' not all equal; constant weights leave means and coefficients unchanged.
#'
#' @param model A fitted model.
#'
#' @return `TRUE` when the model has prior weights that are not all equal.
#'
#' @keywords internal
#' @noRd
.has_unequal_weights <- function(model) {

  prior_weights <- tryCatch(
    stats::weights(model, type = "prior"),
    error = function(e) NULL
  )

  is.numeric(prior_weights) &&
    length(unique(prior_weights[!is.na(prior_weights)])) > 1L
}


#' Build a reference grid for a set of focal factors
#'
#' Internal helper shared by `contrasts_onepred()` and `contrasts_multpred()`.
#'
#' @param model A fitted model supported by the package.
#' @param data Data frame containing the variables of the fixed-effect
#'   formula.
#' @param focal_variables Names of the focal factors.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{grid}{Data frame of prediction rows.}
#'     \item{weights}{Weight of each row within its focal cell; the weights
#'       of the rows of one focal cell sum to one.}
#'     \item{conditional_variables}{Factors weighted by the intercept row of
#'       their hypothesis matrix because they interact with a focal factor.}
#'     \item{interacting_numeric}{Numeric predictors that interact with a
#'       focal factor. The focal coefficients are defined at zero for these
#'       predictors, whereas the grid uses their mean.}
#'   }
#'
#' @keywords internal
#' @noRd
.reference_grid <- function(model,
                            data,
                            focal_variables) {

  predictor_variables <- .fixed_predictor_variables(model)

  missing_variables <- setdiff(predictor_variables, names(data))

  if (length(missing_variables) > 0L) {

    stop(
      paste0(
        "The following variables of the model formula were not found in ",
        "'data': ",
        paste(missing_variables, collapse = ", "),
        ". Supply the data used to fit the model; for transformed terms such ",
        "as log(x), the untransformed variable x is required."
      ),
      call. = FALSE
    )
  }

  for (variable_name in focal_variables) {

    if (!is.factor(data[[variable_name]])) {

      stop(
        paste0("The focal variable '", variable_name, "' must be a factor."),
        call. = FALSE
      )
    }
  }

  # ----------------------------------------------------------
  # VARIABLES THAT INTERACT WITH THE FOCAL FACTORS
  # ----------------------------------------------------------

  fixed_formula <- .extract_fixed_formula(model)

  term_labels <- attr(stats::terms(fixed_formula), "term.labels")

  term_variables <- lapply(
    term_labels,
    function(term_label) all.vars(str2lang(term_label))
  )

  interacting_variables <- unique(
    unlist(
      lapply(
        term_variables,
        function(variables) {
          if (length(variables) > 1L && any(focal_variables %in% variables)) variables
        }
      )
    )
  )

  interacting_variables <- setdiff(interacting_variables, focal_variables)

  model_levels <- .extract_xlevels(model)

  other_variables <- setdiff(predictor_variables, focal_variables)

  # ----------------------------------------------------------
  # LEVELS AND WEIGHTS OF EACH VARIABLE
  # ----------------------------------------------------------

  grid_levels <- list()
  grid_weights <- list()
  constant_values <- list()
  conditional_variables <- character()
  interacting_numeric <- character()

  for (variable_name in focal_variables) {
    grid_levels[[variable_name]] <- levels(data[[variable_name]])
  }

  for (variable_name in other_variables) {

    variable <- data[[variable_name]]

    if (is.numeric(variable) && is.null(dim(variable))) {

      constant_values[[variable_name]] <- mean(variable, na.rm = TRUE)

      if (variable_name %in% interacting_variables) {
        interacting_numeric <- c(interacting_numeric, variable_name)
      }

      next
    }

    factor_variable <- if (is.factor(variable)) variable else factor(variable)

    variable_levels <- model_levels[[variable_name]]

    if (is.null(variable_levels)) {
      variable_levels <- levels(factor_variable)
    }

    level_weights <- rep(1 / length(variable_levels), length(variable_levels))

    if (variable_name %in% interacting_variables) {

      coding_matrix <- stats::contrasts(factor_variable)

      if (!is.null(coding_matrix) &&
          nrow(coding_matrix) == length(variable_levels)) {

        intercept_row <- solve(cbind(1, unname(as.matrix(coding_matrix))))[1L, ]

        intercept_row[abs(intercept_row) < 1e-10] <- 0

        # A negative weight cannot define an average of predictions; the
        # balanced average is kept in that unusual case.
        if (all(intercept_row >= 0)) {

          level_weights <- intercept_row / sum(intercept_row)

          if (!isTRUE(all.equal(level_weights, rep(1 / length(level_weights), length(level_weights))))) {
            conditional_variables <- c(conditional_variables, variable_name)
          }
        }
      }
    }

    grid_levels[[variable_name]] <- variable_levels
    grid_weights[[variable_name]] <- stats::setNames(level_weights, variable_levels)
  }

  # ----------------------------------------------------------
  # GRID
  # ----------------------------------------------------------

  grid <- expand.grid(
    grid_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  row_weights <- rep(1, nrow(grid))

  for (variable_name in names(grid_weights)) {
    row_weights <- row_weights * grid_weights[[variable_name]][grid[[variable_name]]]
  }

  keep_rows <- row_weights > 0

  grid <- grid[keep_rows, , drop = FALSE]
  row_weights <- unname(row_weights[keep_rows])

  for (variable_name in names(grid_levels)) {

    original <- data[[variable_name]]

    grid[[variable_name]] <- if (is.logical(original)) {
      as.logical(grid[[variable_name]])
    } else {
      factor(
        grid[[variable_name]],
        levels = grid_levels[[variable_name]],
        ordered = is.ordered(original)
      )
    }
  }

  for (variable_name in names(constant_values)) {
    grid[[variable_name]] <- constant_values[[variable_name]]
  }

  rownames(grid) <- NULL

  list(
    grid = grid,
    weights = row_weights,
    focal_variables = focal_variables,
    conditional_variables = conditional_variables,
    interacting_numeric = interacting_numeric
  )
}


#' Weights of the grid rows belonging to a focal cell or group of cells
#'
#' @param grid_information Output of `.reference_grid()`.
#' @param conditions Named list mapping focal variables to one level each.
#'
#' @return Numeric weights over the grid rows, summing to one.
#'
#' @keywords internal
#' @noRd
.grid_cell_weights <- function(grid_information,
                               conditions) {

  in_cell <- rep(TRUE, nrow(grid_information$grid))

  for (variable_name in names(conditions)) {
    in_cell <- in_cell &
      as.character(grid_information$grid[[variable_name]]) == as.character(conditions[[variable_name]])
  }

  weights <- grid_information$weights * in_cell

  if (sum(weights) <= 0) {

    stop(
      "A combination of focal levels has no row in the reference grid.",
      call. = FALSE
    )
  }

  weights / sum(weights)
}


#' Summarize predictions with arbitrary weights over the grid rows
#'
#' When a single row is involved, its own estimate and confidence interval
#' are returned, so that intervals transformed from the link scale are
#' preserved. Otherwise the weighted average is computed with the complete
#' fixed-effect covariance matrix.
#'
#' @keywords internal
#' @noRd
.summarize_grid_predictions <- function(predictions,
                                        weights,
                                        prediction_scale,
                                        conf_level) {

  active_rows <- which(weights > 0)

  if (length(active_rows) == 1L) {

    return(
      data.frame(
        estimate = predictions$estimate[active_rows],
        se = predictions$se[active_rows],
        ci_low = predictions$ci_low[active_rows],
        ci_high = predictions$ci_high[active_rows],
        stringsAsFactors = FALSE
      )
    )
  }

  .combine_fixed_predictions(
    predictions = predictions,
    weights = weights,
    prediction_scale = prediction_scale,
    conf_level = conf_level
  )
}
