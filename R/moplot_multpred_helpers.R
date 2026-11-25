#' Extract interaction terms and associated coefficients
#'
#' @param model A fitted model.
#'
#' Only interactions among factors are decomposed; interactions involving a
#' numeric variable are left out of the map and listed, with their numeric
#' variables, in the attribute `"unsupported"`.
#'
#' @return A data frame containing interaction terms, their order, and the
#'   corresponding model coefficients.
#'
#' @keywords internal
#' @noRd
.extract_moplot_interaction_map <- function(model) {

  model_terms <- tryCatch(stats::terms(model), error = function(e) NULL)
  model_matrix <- tryCatch(stats::model.matrix(model), error = function(e) NULL)

  if (is.null(model_terms) || is.null(model_matrix)) {
    stop("The interaction structure could not be extracted from `model`.", call. = FALSE)
  }

  term_labels <- attr(model_terms, "term.labels")
  assign <- attr(model_matrix, "assign")

  if (is.null(term_labels) || is.null(assign)) {
    stop("The interaction structure could not be extracted from `model`.", call. = FALSE)
  }

  interaction_indices <- which(grepl(":", term_labels, fixed = TRUE))

  if (length(interaction_indices) == 0L) {
    stop("No interaction terms were found in `model`.", call. = FALSE)
  }

  # Interactions involving a numeric variable are not decomposed.
  model_frame <- tryCatch(.extract_model_frame(model), error = function(e) NULL)
  term_factors <- attr(model_terms, "factors")

  unsupported <- list()

  if (is.data.frame(model_frame) && is.matrix(term_factors)) {

    for (term_index in interaction_indices) {

      term_variables <- rownames(term_factors)[term_factors[, term_index] > 0]
      term_variables <- gsub("`", "", term_variables, fixed = TRUE)

      numeric_variables <- term_variables[
        !vapply(
          term_variables,
          function(variable) is.factor(model_frame[[variable]]),
          logical(1)
        )
      ]

      if (length(numeric_variables) > 0L) {
        unsupported[[term_labels[term_index]]] <- numeric_variables
      }
    }

    interaction_indices <- interaction_indices[!term_labels[interaction_indices] %in% names(unsupported)]

    if (length(interaction_indices) == 0L) {
      stop(
        paste0(
          "No interactions among factors were found in `model`. Interactions ",
          "involving numeric variables (",
          paste(unique(unlist(unsupported)), collapse = ", "),
          ") are not supported."
        ),
        call. = FALSE
      )
    }
  }

  interaction_map <- lapply(
    interaction_indices,
    function(term_index) {

      interaction <- term_labels[term_index]
      coefficients <- colnames(model_matrix)[assign == term_index]

      if (length(coefficients) == 0L) {
        return(NULL)
      }

      data.frame(
        interaction = interaction,
        interaction_order = lengths(strsplit(interaction, ":", fixed = TRUE)),
        coefficient = coefficients,
        stringsAsFactors = FALSE
      )
    }
  )

  interaction_map <- interaction_map[!vapply(interaction_map, is.null, logical(1))]

  if (length(interaction_map) == 0L) {
    stop("No interaction coefficients were found in `model`.", call. = FALSE)
  }

  interaction_map <- do.call(rbind, interaction_map)
  rownames(interaction_map) <- NULL

  attr(interaction_map, "unsupported") <- unsupported

  interaction_map
}

#' Print available interaction and coefficient selections
#'
#' @param interaction_map Output returned by `.extract_moplot_interaction_map()`.
#'
#' @return The information table invisibly.
#'
#' @keywords internal
#' @noRd
.print_moplot_interaction_info <- function(interaction_map) {

  interactions <- unique(interaction_map$interaction)

  info_table <- do.call(
    rbind,
    lapply(
      interactions,
      function(current_interaction) {

        current_coefficients <- interaction_map$coefficient[interaction_map$interaction == current_interaction]

        interaction_code <- paste0(
          'interaction = "',
          current_interaction,
          '"'
        )

        coefficient_code <- if (length(current_coefficients) == 1L) {
          paste0('coefficients = "', current_coefficients, '"')
        } else {
          paste0(
            'coefficients = c("',
            paste(current_coefficients, collapse = '", "'),
            '")'
          )
        }

        data.frame(
          `Plot full interaction` = interaction_code,
          `Select specific coefficients` = coefficient_code,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  cat("\nAvailable interaction selections:\n\n")

  print(
    info_table,
    row.names = FALSE,
    right = FALSE
  )

  cat("\n")

  invisible(info_table)
}



#' Resolve interaction and coefficient selection
#'
#' @param model A fitted model.
#' @param interaction Optional interaction term.
#' @param coefficients Optional coefficient vector.
#' @param interaction_info Logical. Should available selections be printed?
#'
#' @return A list containing the selected interaction, coefficients, and
#'   complete interaction map.
#'
#' @keywords internal
#' @noRd
.resolve_moplot_interaction_selection <- function(
    model,
    interaction = NULL,
    coefficients = NULL,
    interaction_info = TRUE
) {

  interaction_map <- .extract_moplot_interaction_map(model)

  if (!is.logical(interaction_info) || length(interaction_info) != 1L || is.na(interaction_info)) {
    stop("`interaction_info` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  available_interactions <- unique(interaction_map$interaction)
  available_coefficients <- interaction_map$coefficient

  if (!is.null(interaction)) {

    if (!is.character(interaction) || length(interaction) != 1L || is.na(interaction) || !nzchar(interaction)) {
      stop("`interaction` must be `NULL` or a single non-empty character string.", call. = FALSE)
    }

    unsupported <- attr(interaction_map, "unsupported")

    if (interaction %in% names(unsupported)) {
      stop(
        paste0(
          "Interaction `",
          interaction,
          "` involves numeric variables (",
          paste(unsupported[[interaction]], collapse = ", "),
          "); only interactions among factors are supported."
        ),
        call. = FALSE
      )
    }

    if (!interaction %in% available_interactions) {
      if (isTRUE(interaction_info)) {
        .print_moplot_interaction_info(interaction_map)
      }

      stop(
        paste0(
          "Interaction `",
          interaction,
          "` was not found in the fitted model."
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(coefficients)) {

    if (!is.character(coefficients) || length(coefficients) == 0L || anyNA(coefficients) || any(!nzchar(coefficients))) {
      stop("`coefficients` must be `NULL` or a non-empty character vector.", call. = FALSE)
    }

    unknown_coefficients <- setdiff(coefficients,available_coefficients)

    if (length(unknown_coefficients) > 0L) {
      if (isTRUE(interaction_info)) {
        .print_moplot_interaction_info(interaction_map)
      }

      stop(
        paste0(
          "The following coefficients were not found among the model interaction coefficients: ",
          paste(unknown_coefficients, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }
  }

  # --------------------------------------------------------------------------
  # Case 1: neither interaction nor coefficients supplied
  # --------------------------------------------------------------------------

  if (is.null(interaction) && is.null(coefficients)) {

    maximum_order <- max(interaction_map$interaction_order)

    highest_order_interactions <- unique(interaction_map$interaction[interaction_map$interaction_order == maximum_order])

    if (length(highest_order_interactions) > 1L) {

      .print_moplot_interaction_info(interaction_map)

      stop(
        paste0(
          "More than one interaction of the highest order was found. ",
          "Specify the interaction to display using `interaction = ...`."
        ),
        call. = FALSE
      )
    }

    interaction <- highest_order_interactions[1L]

    coefficients <- interaction_map$coefficient[interaction_map$interaction == interaction]
  }

  # --------------------------------------------------------------------------
  # Case 2: interaction supplied, coefficients omitted
  # --------------------------------------------------------------------------

  if (!is.null(interaction) && is.null(coefficients)) {

    coefficients <- interaction_map$coefficient[interaction_map$interaction == interaction]
  }

  # --------------------------------------------------------------------------
  # Case 3: coefficients supplied, interaction omitted
  # --------------------------------------------------------------------------

  if (is.null(interaction) && !is.null(coefficients)) {

    coefficient_interactions <- unique(
      interaction_map$interaction[
        match(
          coefficients,
          interaction_map$coefficient
        )
      ]
    )

    if (length(coefficient_interactions) != 1L) {

      if (isTRUE(interaction_info)) {
        .print_moplot_interaction_info(
          interaction_map
        )
      }

      stop(
        paste0(
          "The selected coefficients belong to more than one interaction. ",
          "Select coefficients from a single interaction or specify ",
          "`interaction` explicitly."
        ),
        call. = FALSE
      )
    }

    interaction <- coefficient_interactions[1L]
  }

  # --------------------------------------------------------------------------
  # Case 4: both supplied — verify consistency
  # --------------------------------------------------------------------------

  interaction_coefficients <- interaction_map$coefficient[interaction_map$interaction == interaction]

  invalid_for_interaction <- setdiff(coefficients,interaction_coefficients)

  if (length(invalid_for_interaction) > 0L) {

    if (isTRUE(interaction_info)) {
      .print_moplot_interaction_info(
        interaction_map
      )
    }

    stop(
      paste0(
        "The following coefficients do not belong to interaction `",
        interaction,
        "`: ",
        paste(
          invalid_for_interaction,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  if (isTRUE(interaction_info)) {

    .print_moplot_interaction_info(
      interaction_map
    )

    cat(
      "Current selection:\n",
      "  ",
      paste0(
        'interaction = "',
        interaction,
        '"'
      ),
      "\n",
      "  ",
      if (length(coefficients) == 1L) {
        paste0(
          'coefficients = "',
          coefficients,
          '"'
        )
      } else {
        paste0(
          'coefficients = c("',
          paste(
            coefficients,
            collapse = '", "'
          ),
          '")'
        )
      },
      "\n\n",
      sep = ""
    )
  }

  list(
    interaction = interaction,
    coefficients = coefficients,
    interaction_map = interaction_map
  )
}


#' Build a methodological note for multi-predictor contrasts
#'
#' @param model A fitted model.
#' @param estimate_type Character string specifying whether raw or predicted
#'   quantities are used.
#' @param prediction_scale Character string specifying the prediction scale.
#' @param split_rule Character string specifying the split-selection rule.
#' @param top_split Optional explicitly selected top-level split variable.
#' @param configuration Optional configuration-output mode. Use `NULL` when
#'   configuration information is not relevant to the calling function.
#'
#' @return A single character string containing the methodological note.
#'
#' @keywords internal
#' @noRd
.build_multpred_methodological_note <- function(
    model,
    estimate_type,
    prediction_scale,
    split_rule,
    top_split = NULL,
    configuration = NULL
) {

  methodological_points <- character(0)

  # --------------------------------------------------------------------------
  # Estimate type
  # --------------------------------------------------------------------------

  if (identical(estimate_type, "observed")) {

    methodological_points <- c(
      methodological_points,
      "estimate_type='observed' returns observed aggregated means."
    )

  } else if (identical(estimate_type, "predicted")) {

    methodological_points <- c(
      methodological_points,
      paste0(
        "estimate_type='predicted' returns model-based predictions computed ",
        "from the fixed effects and their variance-covariance matrix."
      ),
      paste0(
        "Predictions are computed on a reference grid: numeric predictors are ",
        "fixed at their mean, and the other factors are averaged with equal ",
        "weights, or held at the levels at which the coefficient is defined ",
        "when they interact with its variables."
      )
    )
  }

  # --------------------------------------------------------------------------
  # Mixed-effects predictions
  # --------------------------------------------------------------------------

  is_mixed_model <- inherits(model,c("merMod","glmmTMB"))

  if (
    identical(estimate_type, "predicted") &&
    is_mixed_model
  ) {

    methodological_points <- c(
      methodological_points,
      paste0(
        "Random effects are not included in the predictions; ",
        "the displayed predictions are therefore population-level estimates."
      )
    )
  }

  # --------------------------------------------------------------------------
  # Model family
  # --------------------------------------------------------------------------

  family_object <- tryCatch(
    .extract_family(model),
    error = function(e) NULL
  )

  is_generalized_model <- !is.null(family_object) &&
    !identical(
      tolower(
        family_object$family
      ),
      "gaussian"
    )

  # --------------------------------------------------------------------------
  # Response-scale generalized predictions
  # --------------------------------------------------------------------------

  if (
    identical(estimate_type, "predicted") &&
    is_generalized_model &&
    identical(prediction_scale, "response")
  ) {

    methodological_points <- c(
      methodological_points,
      paste0(
        "For response-scale generalized-model predictions, the function ",
        "averages inverse-link-transformed predictions and therefore computes ",
        "mean(linkinv(eta)), not linkinv(mean(eta))."
      )
    )
  }

  # --------------------------------------------------------------------------
  # Prediction uncertainty
  # --------------------------------------------------------------------------

  if (identical(estimate_type, "predicted")) {

    uncertainty_note <- paste0(
      "Cell-level and aggregated predicted standard errors are computed using ",
      "the complete fixed-effect variance-covariance matrix."
    )

    if (
      is_generalized_model &&
      identical(prediction_scale, "response")
    ) {

      uncertainty_note <- paste0(
        uncertainty_note,
        " For response-scale generalized models, uncertainty is propagated ",
        "through the multivariate delta method, and confidence intervals are ",
        "computed on the link scale and transformed with the inverse link."
      )
    }

    uncertainty_note <- paste0(
      uncertainty_note,
      " Predicted cell estimates are not treated as independent during aggregation."
    )

    methodological_points <- c(methodological_points,uncertainty_note)
  }

  # --------------------------------------------------------------------------
  # glmmTMB
  # --------------------------------------------------------------------------

  if (inherits(model, "glmmTMB")) {

    methodological_points <- c(
      methodological_points,
      "For glmmTMB objects, only the conditional model component is analyzed."
    )
  }

  # --------------------------------------------------------------------------
  # lmer p-values
  # --------------------------------------------------------------------------

  if (
    inherits(model, "lmerMod") &&
    requireNamespace(
      "lmerTest",
      quietly = TRUE
    )
  ) {

    methodological_points <- c(
      methodological_points,
      paste0(
        "For lmer models fitted with lme4, p-values are obtained through ",
        "lmerTest."
      )
    )
  }

  # --------------------------------------------------------------------------
  # Configuration representation
  # --------------------------------------------------------------------------

  if (!is.null(configuration)) {

    if (identical(configuration, "all")) {

      methodological_points <- c(
        methodological_points,
        paste0(
          "configuration_explicit shows the full hierarchical composition, ",
          "whereas configuration_semantic provides a compact description."
        )
      )

    } else if (identical(configuration, "explicit")) {

      methodological_points <- c(
        methodological_points,
        "configuration_explicit reports the full hierarchical composition."
      )

    } else if (identical(configuration, "semantic")) {

      methodological_points <- c(
        methodological_points,
        paste0(
          "configuration_semantic provides a compact description of the ",
          "hierarchical composition."
        )
      )
    }
  }

  # --------------------------------------------------------------------------
  # Split-variable selection
  # --------------------------------------------------------------------------

  if (!is.null(top_split)) {

    methodological_points <- c(
      methodological_points,
      paste0(
        "The top-level split variable was explicitly set to '",
        top_split,
        "' using top_split."
      )
    )

  } else {

    methodological_points <- c(
      methodological_points,
      paste0(
        "The top-level split variable was selected using split_rule='",
        split_rule,
        "'."
      )
    )
  }

  # --------------------------------------------------------------------------
  # Assemble note
  # --------------------------------------------------------------------------

  numbered_points <- paste0(
    seq_along(methodological_points),
    ") ",
    methodological_points
  )

  paste(
    c(
      "Methodological note:",
      numbered_points
    ),
    collapse = " "
  )
}

#' Prepare multi-predictor MoPlot data for interactions of arbitrary order
#'
#' Experimental internal helper used to construct the hierarchical graphical
#' representation of interaction contrasts returned by `contrasts_multpred()`.
#'
#' @param contrast_data A data frame returned by `contrasts_multpred()`.
#' @param within_split_gap Horizontal distance between the two quantities
#'   forming an elementary contrast.
#' @param between_split_gap Base horizontal distance between neighbouring
#'   contrast blocks. Higher-order blocks receive proportionally larger gaps.
#' @param between_coefficient_gap Horizontal distance between coefficient
#'   blocks.
#'
#' @return A list containing data frames used to construct an n-way
#'   interaction plot.
#'
#' @keywords internal
#' @noRd
.prepare_moplot_multpred_data <- function(
    contrast_data,
    within_split_gap = 0.75,
    between_split_gap = 0.35,
    between_coefficient_gap = 0.75
) {

  # ==========================================================================
  # 1. BASIC VALIDATION
  # ==========================================================================

  if (!is.data.frame(contrast_data)) {
    stop(
      "`contrast_data` must be a data frame returned by `contrasts_multpred()`.",
      call. = FALSE
    )
  }

  required_columns <- c(
    "coefficient",
    "contrast_id",
    "significance",
    "p_value",
    "split_var",
    "split_side",
    "estimate",
    "se",
    "ci_low",
    "ci_high",
    "configuration_explicit",
    "configuration_semantic"
  )

  missing_columns <- setdiff(
    required_columns,
    names(contrast_data)
  )

  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        "`contrast_data` is missing the required column",
        if (length(missing_columns) > 1L) "s: " else ": ",
        paste(missing_columns, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  if (nrow(contrast_data) == 0L) {
    stop(
      "`contrast_data` contains no contrasts to display.",
      call. = FALSE
    )
  }

  gap_values <- c(within_split_gap,between_split_gap,between_coefficient_gap)

  if (
    !is.numeric(gap_values) ||
    any(
      lengths(
        list(
          within_split_gap,
          between_split_gap,
          between_coefficient_gap
        )
      ) != 1L
    ) ||
    any(is.na(gap_values)) ||
    any(!is.finite(gap_values)) ||
    within_split_gap <= 0 ||
    between_split_gap < 0 ||
    between_coefficient_gap < 0
  ) {
    stop(
      paste0(
        "`within_split_gap` must be a positive finite number, and ",
        "`between_split_gap` and `between_coefficient_gap` must be ",
        "non-negative finite numbers."
      ),
      call. = FALSE
    )
  }

  decomposition <- attr(contrast_data,"decomposition")

  if (
    is.null(decomposition) ||
    !is.data.frame(decomposition) ||
    nrow(decomposition) == 0L
  ) {
    stop(
      paste0(
        "`contrast_data` does not contain the hierarchical decomposition ",
        "required by `.prepare_moplot_multpred_data()`."
      ),
      call. = FALSE
    )
  }

  required_decomposition_columns <- c(
    "coefficient",
    "contrast_id",
    "split_side",
    "inner_side",
    "decomposition_path",
    "decomposition_depth",
    "variable_path",
    "estimate",
    "se",
    "ci_low",
    "ci_high"
  )

  missing_decomposition_columns <- setdiff(required_decomposition_columns,names(decomposition))

  if (length(missing_decomposition_columns) > 0L) {
    stop(
      paste0(
        "The hierarchical decomposition is missing the required column",
        if (length(missing_decomposition_columns) > 1L) "s: " else ": ",
        paste(missing_decomposition_columns, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  # ==========================================================================
  # 2. IDENTIFY PREDICTOR COLUMNS
  # ==========================================================================

  structural_columns <- c(
    "coefficient",
    "contrast_id",
    "significance",
    "p_value",
    "split_var",
    "split_side",
    "n",
    "estimate",
    "se",
    "ci_low",
    "ci_high",
    "configuration_explicit",
    "configuration_semantic"
  )

  predictor_columns <- setdiff(names(contrast_data),structural_columns)

  if (length(predictor_columns) < 2L) {
    stop(
      paste0(
        "`.prepare_moplot_multpred_data()` requires at least two ",
        "categorical predictors."
      ),
      call. = FALSE
    )
  }

  split_variables <- unique(contrast_data$split_var)

  split_variables <- split_variables[!is.na(split_variables) & nzchar(split_variables)]

  if (length(split_variables) != 1L) {
    stop(
      "Exactly one top-level split variable must be identified.",
      call. = FALSE
    )
  }

  split_var <- split_variables[1L]

  # ==========================================================================
  # 3. INTERNAL PATH HELPERS
  # ==========================================================================

  split_path <- function(path) {
    if (
      is.na(path) ||
      !nzchar(path) ||
      identical(path, "root")
    ) {
      return(character(0))
    }

    strsplit(
      path,
      "_",
      fixed = TRUE
    )[[1L]]
  }

  collapse_path <- function(tokens) {
    if (length(tokens) == 0L) {
      return("root")
    }

    paste(
      tokens,
      collapse = "_"
    )
  }

  parent_path <- function(path) {
    tokens <- split_path(path)

    if (length(tokens) <= 1L) {
      return("root")
    }

    collapse_path(tokens[-length(tokens)])
  }

  common_prefix_length <- function(path_a, path_b) {
    a <- split_path(path_a)

    b <- split_path(path_b)

    n <- min(length(a),length(b))

    if (n == 0L) {
      return(0L)
    }

    same <- a[seq_len(n)] == b[seq_len(n)]

    mismatch <- which(!same)

    if (length(mismatch) == 0L) {
      return(n)
    }

    mismatch[1L] - 1L
  }

  path_order_key <- function(path) {
    tokens <- split_path(path)

    paste0(
      ifelse(
        tokens == "positive",
        "0",
        "1"
      ),
      collapse = ""
    )
  }

  # ==========================================================================
  # 4. PRESERVE COEFFICIENT ORDER
  # ==========================================================================

  coefficient_order <- unique(contrast_data$coefficient)

  point_parts <- vector("list",length(coefficient_order))

  segment_parts <- vector("list",length(coefficient_order))

  delta_parts <- vector("list",length(coefficient_order))

  coefficient_parts <- vector("list",length(coefficient_order))

  group_label_parts <- vector("list",length(coefficient_order))

  variable_paths <- vector("list",length(coefficient_order))

  hierarchy_depths <- integer(length(coefficient_order))

  current_x <- 1

  # ==========================================================================
  # 5. PROCESS EACH COEFFICIENT
  # ==========================================================================

  for (coefficient_index in seq_along(coefficient_order)) {

    coefficient_name <- coefficient_order[coefficient_index]

    coefficient_rows <- contrast_data[contrast_data$coefficient == coefficient_name, ,drop = FALSE]

    leaf_rows <- decomposition[decomposition$coefficient == coefficient_name, ,drop = FALSE]

    depths <- unique(leaf_rows$decomposition_depth)

    depths <- depths[!is.na(depths)]

    if (length(depths) != 1L) {
      stop(
        paste0(
          "Coefficient `",
          coefficient_name,
          "` does not have a unique decomposition depth."
        ),
        call. = FALSE
      )
    }

    depth <- as.integer(depths[1L])

    if (depth < 2L) {
      stop(
        paste0(
          "Coefficient `",
          coefficient_name,
          "` does not represent a multi-predictor interaction."
        ),
        call. = FALSE
      )
    }

    coefficient_variable_paths <- unique(leaf_rows$variable_path)

    coefficient_variable_paths <- coefficient_variable_paths[!is.na(coefficient_variable_paths) & nzchar(coefficient_variable_paths)]

    if (length(coefficient_variable_paths) != 1L) {
      stop(
        paste0(
          "Coefficient `",
          coefficient_name,
          "` does not have a unique variable hierarchy."
        ),
        call. = FALSE
      )
    }

    variable_path <- strsplit(
      coefficient_variable_paths[1L],
      .moplot_variable_separator,
      fixed = TRUE
    )[[1L]]

    if (length(variable_path) != depth) {
      stop(
        paste0(
          "Coefficient `",
          coefficient_name,
          "` has inconsistent variable and decomposition depths."
        ),
        call. = FALSE
      )
    }

    if (!all(variable_path %in% predictor_columns)) {
      stop(
        paste0(
          "Coefficient `",
          coefficient_name,
          "` contains variables that were not identified among the predictors."
        ),
        call. = FALSE
      )
    }

    hierarchy_depths[coefficient_index] <- depth

    variable_paths[[coefficient_index]] <- variable_path

    target_var <- variable_path[depth]

    top_split_var <- variable_path[1L]

    # ------------------------------------------------------------------------
    # 5A. ORDER LEAVES
    # ------------------------------------------------------------------------

    path_lengths <- lengths(lapply(leaf_rows$decomposition_path,split_path))

    if (any(path_lengths != depth)) {
      stop(
        paste0(
          "Coefficient `",
          coefficient_name,
          "` contains decomposition paths of inconsistent depth."
        ),
        call. = FALSE
      )
    }

    if (anyDuplicated(leaf_rows$decomposition_path)) {
      stop(
        paste0(
          "Coefficient `",
          coefficient_name,
          "` contains duplicated decomposition paths."
        ),
        call. = FALSE
      )
    }

    leaf_order <- order(
      vapply(
        leaf_rows$decomposition_path,
        path_order_key,
        character(1)
      )
    )

    leaf_rows <- leaf_rows[leaf_order, ,drop = FALSE]

    # ------------------------------------------------------------------------
    # 5B. ASSIGN HORIZONTAL POSITIONS
    # ------------------------------------------------------------------------

    leaf_x <- numeric(nrow(leaf_rows))

    leaf_x[1L] <- current_x

    if (nrow(leaf_rows) > 1L) {
      for (leaf_index in 2:nrow(leaf_rows)) {

        previous_path <- leaf_rows$decomposition_path[leaf_index - 1L]

        current_path <- leaf_rows$decomposition_path[leaf_index]

        common_depth <- common_prefix_length(previous_path,current_path)

        if (common_depth == depth - 1L) {
          current_gap <- within_split_gap
        } else {
          separation_level <- max(1L,depth - common_depth - 1L)

          current_gap <- between_split_gap * separation_level
        }

        leaf_x[leaf_index] <- leaf_x[leaf_index - 1L] + current_gap
      }
    }

    leaf_rows$x <- leaf_x
    leaf_rows$coefficient_index <- coefficient_index
    leaf_rows$target_group <- as.character(leaf_rows[[target_var]])
    leaf_rows$point_label <- as.character(leaf_rows$configuration_semantic)
    leaf_rows$split_group <- as.character(leaf_rows[[top_split_var]])

    top_split_paths <- vapply(
      leaf_rows$decomposition_path,
      function(path) {
        split_path(path)[1L]
      },
      character(1)
    )

    leaf_rows$split_index <- match(top_split_paths,unique(top_split_paths))

    point_parts[[coefficient_index]] <- leaf_rows

    # ------------------------------------------------------------------------
    # 5C. BUILD ELEMENTARY CONTRASTS
    # ------------------------------------------------------------------------

    elementary_parent_paths <- vapply(
      leaf_rows$decomposition_path,
      parent_path,
      character(1)
    )

    elementary_paths <- unique(unname(elementary_parent_paths))

    elementary_segments <- vector("list",length(elementary_paths))

    elementary_nodes <- vector("list",length(elementary_paths))

    for (elementary_index in seq_along(elementary_paths)) {

      node_path <- elementary_paths[elementary_index]

      child_rows <- leaf_rows[unname(elementary_parent_paths) == node_path, ,drop = FALSE]

      child_sides <- unname(
        vapply(
          child_rows$decomposition_path,
          function(path) {
            utils::tail(
              split_path(path),
              1L
            )
          },
          character(1)
        )
      )

      if (
        nrow(child_rows) != 2L ||
        length(child_sides) != 2L ||
        anyNA(child_sides) ||
        !setequal(
          child_sides,
          c(
            "positive",
            "negative"
          )
        )
      ) {
        stop(
          paste0(
            "Elementary node `",
            node_path,
            "` for coefficient `",
            coefficient_name,
            "` does not contain exactly one positive and one negative child."
          ),
          call. = FALSE
        )
      }

      child_order <- match(
        c(
          "positive",
          "negative"
        ),
        child_sides
      )

      child_rows <- child_rows[child_order, ,drop = FALSE]

      delta_number <- elementary_index

      delta_symbol <- paste0(
        "Delta[",
        delta_number,
        "]"
      )

      grouping_var <- variable_path[depth - 1L]

      grouping_value <- unique(
        as.character(
          child_rows[[grouping_var]]
        )
      )

      grouping_value <- if (length(grouping_value) == 1L) {
        grouping_value
      } else {
        collapse_path(
          split_path(node_path)
        )
      }

      node_tokens <- split_path(node_path)

      top_side <- if (length(node_tokens) > 0L) {
        node_tokens[1L]
      } else {
        NA_character_
      }

      elementary_segments[[
        elementary_index
      ]] <- data.frame(
        coefficient = coefficient_name,
        contrast_id = coefficient_rows$contrast_id[1L],
        coefficient_index = coefficient_index,
        split_index = elementary_index,
        split_group = grouping_value,
        split_side = top_side,
        node_path = node_path,
        node_level = 1L,
        x = child_rows$x[1L],
        xend = child_rows$x[2L],
        y = child_rows$estimate[1L],
        yend = child_rows$estimate[2L],
        delta_number = delta_number,
        delta_label = delta_symbol,
        delta_x = mean(
          child_rows$x
        ),
        delta_y = mean(
          child_rows$estimate
        ),
        positive_child = child_rows$decomposition_path[1L],
        negative_child = child_rows$decomposition_path[2L],
        stringsAsFactors = FALSE
      )

      elementary_nodes[[
        elementary_index
      ]] <- data.frame(
        coefficient = coefficient_name,
        contrast_id = coefficient_rows$contrast_id[1L],
        coefficient_index = coefficient_index,
        node_path = node_path,
        node_depth = depth - 1L,
        node_level = 1L,
        x = mean(
          child_rows$x
        ),
        y = mean(
          child_rows$estimate
        ),
        delta_number = delta_number,
        symbol_label = delta_symbol,
        expression = NA_character_,
        positive_child = child_rows$decomposition_path[1L],
        negative_child = child_rows$decomposition_path[2L],
        is_final = FALSE,
        stringsAsFactors = FALSE
      )
    }

    coefficient_segments <- do.call(rbind,elementary_segments)

    current_nodes <- do.call(rbind,elementary_nodes)

    next_delta_number <- nrow(current_nodes) + 1L

    higher_node_parts <- list()

    # ------------------------------------------------------------------------
    # 5D. BUILD HIGHER-ORDER CONTRAST NODES
    # ------------------------------------------------------------------------

    for (node_level in 2:depth) {

      current_parent_paths <- vapply(
        current_nodes$node_path,
        parent_path,
        character(1)
      )

      parent_paths <- unique(unname(current_parent_paths))

      next_nodes <- vector("list",length(parent_paths))

      for (parent_index in seq_along(parent_paths)) {

        node_path <- parent_paths[parent_index]

        child_nodes <- current_nodes[unname(current_parent_paths) == node_path, ,drop = FALSE]

        child_sides <- unname(
          vapply(
            child_nodes$node_path,
            function(path) {
              utils::tail(
                split_path(path),
                1L
              )
            },
            character(1)
          )
        )

        if (
          nrow(child_nodes) != 2L ||
          length(child_sides) != 2L ||
          anyNA(child_sides) ||
          !setequal(
            child_sides,
            c(
              "positive",
              "negative"
            )
          )
        ) {
          stop(
            paste0(
              "Hierarchical node `",
              node_path,
              "` for coefficient `",
              coefficient_name,
              "` does not contain exactly one positive and one negative child."
            ),
            call. = FALSE
          )
        }

        child_order <- match(
          c(
            "positive",
            "negative"
          ),
          child_sides
        )

        child_nodes <- child_nodes[child_order, ,drop = FALSE]

        is_final <- node_level == depth

        if (is_final) {
          delta_number <- NA_integer_
          delta_symbol <- "Delta"
        } else {
          delta_number <- next_delta_number

          delta_symbol <- paste0(
            "Delta[",
            delta_number,
            "]"
          )

          next_delta_number <- next_delta_number + 1L
        }

        expression <- paste0(
          delta_symbol,
          " == ",
          child_nodes$symbol_label[1L],
          " - ",
          child_nodes$symbol_label[2L]
        )

        next_nodes[[
          parent_index
        ]] <- data.frame(
          coefficient = coefficient_name,
          contrast_id = coefficient_rows$contrast_id[1L],
          coefficient_index = coefficient_index,
          node_path = node_path,
          node_depth = depth - node_level,
          node_level = node_level,
          x = mean(
            child_nodes$x
          ),
          y = mean(
            child_nodes$y
          ),
          child1_x = child_nodes$x[1L],
          child1_y = child_nodes$y[1L],
          child2_x = child_nodes$x[2L],
          child2_y = child_nodes$y[2L],
          delta_number = delta_number,
          symbol_label = delta_symbol,
          expression = expression,
          positive_child = child_nodes$node_path[1L],
          negative_child = child_nodes$node_path[2L],
          is_final = is_final,
          stringsAsFactors = FALSE
        )
      }

      current_nodes <- do.call(rbind,next_nodes)

      higher_node_parts[[length(higher_node_parts) + 1L]] <- current_nodes
    }

    coefficient_delta_data <- do.call(rbind,higher_node_parts)

    coefficient_delta_data$label <- coefficient_delta_data$expression

    delta_parts[[coefficient_index]] <- coefficient_delta_data

    # ------------------------------------------------------------------------
    # 5E. FINAL COEFFICIENT POSITION
    # ------------------------------------------------------------------------

    final_node <- coefficient_delta_data[coefficient_delta_data$is_final, ,drop = FALSE]

    if (nrow(final_node) != 1L) {
      stop(
        paste0(
          "Coefficient `",
          coefficient_name,
          "` does not define exactly one final contrast node."
        ),
        call. = FALSE
      )
    }

    coefficient_parts[[
      coefficient_index
    ]] <- data.frame(
      coefficient = coefficient_name,
      contrast_id = coefficient_rows$contrast_id[1L],
      significance = coefficient_rows$significance[1L],
      p_value = coefficient_rows$p_value[1L],
      coefficient_index = coefficient_index,
      x = final_node$x[1L],
      stringsAsFactors = FALSE
    )

    segment_parts[[coefficient_index]] <- coefficient_segments

    # ------------------------------------------------------------------------
    # 5F. BUILD HIERARCHICAL GROUP LABEL POSITIONS
    # ------------------------------------------------------------------------

    grouping_variables <- variable_path[seq_len(depth - 1L)]

    coefficient_group_labels <- list()

    group_counter <- 1L

    for (path_position in seq_along(grouping_variables)) {

      grouping_var <- grouping_variables[
        path_position
      ]

      group_prefix <- vapply(
        leaf_rows$decomposition_path,
        function(path) {

          tokens <- split_path(
            path
          )

          collapse_path(
            tokens[
              seq_len(
                path_position
              )
            ]
          )
        },
        character(1)
      )

      unique_prefixes <- unique(unname(group_prefix))

      for (prefix_index in seq_along(unique_prefixes)) {

        current_prefix <- unique_prefixes[prefix_index]

        current_rows <- leaf_rows[unname(group_prefix) == current_prefix, ,drop = FALSE]

        group_value <- unique(as.character(current_rows[[grouping_var]]))

        group_value <- if (length(group_value) == 1L) {
          group_value
        } else {
          current_prefix
        }

        coefficient_group_labels[[
          group_counter
        ]] <- data.frame(
          coefficient = coefficient_name,
          contrast_id = coefficient_rows$contrast_id[1L],
          coefficient_index = coefficient_index,
          variable = grouping_var,
          group_value = group_value,
          split_group = group_value,
          path_prefix = current_prefix,
          path_position = path_position,
          display_level = depth - path_position,
          x = mean(
            range(
              current_rows$x
            )
          ),
          stringsAsFactors = FALSE
        )

        group_counter <- group_counter + 1L
      }
    }

    group_label_parts[[coefficient_index]] <- do.call(rbind,coefficient_group_labels)

    current_x <- max(leaf_rows$x) + between_coefficient_gap
  }

  # ==========================================================================
  # 6. COMBINE OUTPUT
  # ==========================================================================

  point_data <- do.call(rbind,point_parts)

  segment_data <- do.call(rbind,segment_parts)

  delta_data <- do.call(rbind,delta_parts)

  coefficient_data <- do.call(rbind,coefficient_parts)

  group_label_data <- do.call(rbind,group_label_parts)

  rownames(point_data) <- NULL
  rownames(segment_data) <- NULL
  rownames(delta_data) <- NULL
  rownames(coefficient_data) <- NULL
  rownames(group_label_data) <- NULL

  # ==========================================================================
  # 7. BACKWARD-COMPATIBLE SPLIT LABEL DATA
  # ==========================================================================

  split_label_data <- group_label_data[group_label_data$variable == split_var, ,drop = FALSE]

  # ==========================================================================
  # 8. STORE STRUCTURAL INFORMATION
  # ==========================================================================

  unique_variable_paths <- unique(
    vapply(
      variable_paths,
      paste,
      collapse = .moplot_variable_separator,
      FUN.VALUE = character(1)
    )
  )

  if (length(unique_variable_paths) != 1L) {
    stop(
      paste0(
        "The selected coefficients do not share a common variable hierarchy. ",
        "This experimental implementation currently requires a common ",
        "hierarchy across displayed coefficients."
      ),
      call. = FALSE
    )
  }

  common_variable_path <- variable_paths[[1L]]

  common_depth <- hierarchy_depths[1L]

  if (!all(hierarchy_depths == common_depth)) {
    stop(
      "The selected coefficients do not share a common interaction order.",
      call. = FALSE
    )
  }

  target_var <- common_variable_path[
    common_depth
  ]

  attr(point_data,"split_var") <- split_var

  attr(point_data,"target_var") <- target_var

  attr(point_data,"variable_path") <- common_variable_path

  attr(point_data,"hierarchy_depth") <- common_depth

  x_limits <- c(min(point_data$x) - 0.5,max(point_data$x) + 0.5)

  list(
    point_data = point_data,
    segment_data = segment_data,
    split_label_data = split_label_data,
    group_label_data = group_label_data,
    coefficient_data = coefficient_data,
    delta_data = delta_data,
    split_var = split_var,
    target_var = target_var,
    variable_path = common_variable_path,
    hierarchy_depth = common_depth,
    x_limits = x_limits
  )
}

#' Shortest unique abbreviations of a set of labels
#'
#' Internal helper. A label whose first three characters are unique among the
#' labels is abbreviated as its first half (at least three characters)
#' followed by a period. Labels that share their first three characters, such
#' as `Treatment1` and `Treatment2`, keep a three-character prefix followed by
#' the part that distinguishes them (`Tre.1`, `Tre.2`). A label whose
#' abbreviation would still collide is kept in full.
#'
#' @param labels Character vector of unique labels to abbreviate.
#' @param reserved Labels that must not be produced as abbreviations.
#'
#' @return A character vector of abbreviations, one per label.
#'
#' @keywords internal
#' @noRd
.unique_abbreviations <- function(labels,
                                  reserved = character(0)) {

  # Labels with a unique three-character start: first half of the label.
  abbreviations <- paste0(
    substr(labels, 1L, pmax(3L, ceiling(nchar(labels) / 2))),
    "."
  )

  # Labels sharing their first three characters, such as Treatment1 and
  # Treatment2 or TreatmentHigh and TreatmentLow, would be hard to tell
  # apart by prefixes alone: they keep a three-character prefix followed
  # by the part that distinguishes them within the group.
  start <- substr(labels, 1L, 3L)

  for (shared_start in unique(start[duplicated(start)])) {

    group <- which(start == shared_start)

    group_labels <- labels[group]

    shortest <- min(nchar(group_labels))

    common_length <- 0L

    while (common_length < shortest &&
           length(unique(substr(group_labels, common_length + 1L, common_length + 1L))) == 1L) {

      common_length <- common_length + 1L
    }

    abbreviations[group] <- paste0(
      substr(group_labels, 1L, min(3L, common_length)),
      ".",
      substr(group_labels, common_length + 1L, nchar(group_labels))
    )
  }

  collides <- duplicated(abbreviations) |
    duplicated(abbreviations, fromLast = TRUE) |
    abbreviations %in% reserved

  # Any remaining collision keeps the full label.
  abbreviations[collides] <- labels[collides]

  abbreviations
}


#' Abbreviate long factor labels for graphical display
#'
#' Internal helper used to abbreviate long factor-level labels while
#' preserving short labels and Mean(...) expressions.
#'
#' @param labels Character vector of labels.
#' @param max_chars Maximum number of characters retained without abbreviation.
#'
#' @return A list containing the display labels and the abbreviation key.
#'
#' @keywords internal
#' @noRd
.abbreviate_moplot_labels <- function(labels, max_chars = 5L) {

  labels <- as.character(labels)

  extract_atomic_labels <- function(x) {

    if (is.na(x) || !nzchar(x)) return(character(0))

    if (grepl("^Mean\\(.*\\)$", x)) {

      inside <- sub("^Mean\\((.*)\\)$", "\\1", x)

      return(
        trimws(
          strsplit(
            inside,
            ",",
            fixed = TRUE
          )[[1L]]
        )
      )
    }

    # Weighted groups: "W(level weight, level weight)".
    if (grepl("^W\\(.*\\)$", x)) {

      entries <- trimws(
        strsplit(
          sub("^W\\((.*)\\)$", "\\1", x),
          ",",
          fixed = TRUE
        )[[1L]]
      )

      return(sub("\\s+\\S+$", "", entries))
    }

    x
  }

  atomic_labels <- unique(
    unlist(
      lapply(
        labels,
        extract_atomic_labels
      ),
      use.names = FALSE
    )
  )

  long_labels <- atomic_labels[
    nchar(atomic_labels) > max_chars
  ]

  if (length(long_labels) == 0L) {

    return(
      list(
        labels = labels,
        key = NULL
      )
    )
  }

  abbreviations <- .unique_abbreviations(
    labels = long_labels,
    reserved = setdiff(atomic_labels, long_labels)
  )

  # An abbreviation is useful only if it is shorter than the label.
  shortened <- nchar(abbreviations) < nchar(long_labels)

  long_labels <- long_labels[shortened]
  abbreviations <- abbreviations[shortened]

  if (length(long_labels) == 0L) {

    return(
      list(
        labels = labels,
        key = NULL
      )
    )
  }

  lookup <- stats::setNames(abbreviations,long_labels)

  replace_atomic_labels <- function(x) {

    if (is.na(x) || !nzchar(x)) return(x)

    if (grepl("^Mean\\(.*\\)$", x)) {

      inside <- sub(
        "^Mean\\((.*)\\)$",
        "\\1",
        x
      )

      values <- trimws(
        strsplit(
          inside,
          ",",
          fixed = TRUE
        )[[1L]]
      )

      replace <- values %in% names(lookup)

      values[replace] <- unname(lookup[values[replace]])

      return(
        paste0(
          "Mean(",
          paste(
            values,
            collapse = ","
          ),
          ")"
        )
      )
    }

    # Weighted groups: "W(level weight, level weight)". Level names may
    # contain spaces, so the weight is taken as the last token.
    if (grepl("^W\\(.*\\)$", x)) {

      entries <- trimws(
        strsplit(
          sub("^W\\((.*)\\)$", "\\1", x),
          ",",
          fixed = TRUE
        )[[1L]]
      )

      level_names <- sub("\\s+\\S+$", "", entries)
      weights <- sub("^.*\\s+(\\S+)$", "\\1", entries)

      replace <- level_names %in% names(lookup)

      level_names[replace] <- unname(lookup[level_names[replace]])

      return(
        paste0(
          "W(",
          paste(level_names, weights, collapse = ", "),
          ")"
        )
      )
    }

    if (x %in% names(lookup)) {
      return(
        unname(
          lookup[x]
        )
      )
    }

    x
  }

  display_labels <- vapply(
    labels,
    replace_atomic_labels,
    character(1)
  )

  key <- data.frame(
    abbreviated = unname(abbreviations),
    full = long_labels,
    stringsAsFactors = FALSE
  )

  list(
    labels = display_labels,
    key = key
  )
}

#' Collect abbreviated factor labels used in a multi-predictor MoPlot
#'
#' @param prep Output returned by `.prepare_moplot_multpred_data()`.
#' @param max_chars Maximum number of characters retained without abbreviation.
#'
#' @return A data frame containing the abbreviated and full labels, or `NULL`
#'   when no abbreviations are required.
#'
#' @keywords internal
#' @noRd
.collect_moplot_label_key <- function(
    prep,
    max_chars = 5L
) {

  all_labels <- as.character(
    prep$point_data$target_group
  )

  if (
    !is.null(prep$group_label_data) &&
    nrow(prep$group_label_data) > 0L
  ) {

    all_labels <- c(
      all_labels,
      as.character(
        prep$group_label_data$group_value
      )
    )
  }

  abbreviation_result <-
    .abbreviate_moplot_labels(
      labels = all_labels,
      max_chars = max_chars
    )

  if (
    is.null(abbreviation_result$key) ||
    nrow(abbreviation_result$key) == 0L
  ) {
    return(NULL)
  }

  unique(abbreviation_result$key)
}

#' Apply a MoPlot label abbreviation key
#'
#' @param labels Character vector of labels.
#' @param label_key Abbreviation key returned by `.collect_moplot_label_key()`.
#'
#' @return Character vector of display labels.
#'
#' @keywords internal
#' @noRd
.apply_moplot_label_key <- function(
    labels,
    label_key = NULL
) {

  labels <- as.character(labels)

  if (
    is.null(label_key) ||
    nrow(label_key) == 0L
  ) {
    return(labels)
  }

  lookup <- stats::setNames(
    label_key$abbreviated,
    label_key$full
  )

  replace_one <- function(x) {

    if (
      is.na(x) ||
      !nzchar(x)
    ) {
      return(x)
    }

    if (grepl("^Mean\\(.*\\)$", x)) {

      inside <- sub(
        "^Mean\\((.*)\\)$",
        "\\1",
        x
      )

      values <- trimws(
        strsplit(
          inside,
          ",",
          fixed = TRUE
        )[[1L]]
      )

      replace <- values %in% names(lookup)

      values[replace] <-unname(lookup[values[replace]])

      return(
        paste0(
          "Mean(",
          paste(
            values,
            collapse = ","
          ),
          ")"
        )
      )
    }

    # Weighted groups: "W(level weight, level weight)". Level names may
    # contain spaces, so the weight is taken as the last token.
    if (grepl("^W\\(.*\\)$", x)) {

      entries <- trimws(
        strsplit(
          sub("^W\\((.*)\\)$", "\\1", x),
          ",",
          fixed = TRUE
        )[[1L]]
      )

      level_names <- sub("\\s+\\S+$", "", entries)
      weights <- sub("^.*\\s+(\\S+)$", "\\1", entries)

      replace <- level_names %in% names(lookup)

      level_names[replace] <- unname(lookup[level_names[replace]])

      return(
        paste0(
          "W(",
          paste(level_names, weights, collapse = ", "),
          ")"
        )
      )
    }

    if (x %in% names(lookup)) {
      return(
        unname(
          lookup[x]
        )
      )
    }

    x
  }

  vapply(
    labels,
    replace_one,
    character(1)
  )
}
#' Plot the upper panel of an n-way multi-predictor MoPlot
#'
#' Experimental internal helper used with
#' `.prepare_moplot_multpred_data()`.
#'
#' @param prep Output returned by `.prepare_moplot_multpred_data()`.
#'
#' @param style Character string specifying the graphical style. Available
#'   options are `"minimal"` and `"tree"`.
#'
#' @param base_size Base font size.
#'
#' @param label_key Abbreviation key returned by `.collect_moplot_label_key()`.
#'
#' @param y_title Title of the y axis, usually the name of the response.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd
.plot_moplot_multpred_upper <- function(
    prep,
    style = c("minimal", "tree"),
    base_size = 13,
    label_key=NULL,
    y_title = NULL
) {

  # ==========================================================================
  # 1. ARGUMENT MATCHING
  # ==========================================================================

  style <- match.arg(style)


  # ==========================================================================
  # 2. DATA
  # ==========================================================================

  point_data <- prep$point_data
  segment_data <- prep$segment_data
  delta_data <- prep$delta_data
  group_label_data <- prep$group_label_data
  hierarchy_depth <- prep$hierarchy_depth

  # ==========================================================================
  # 2.1. X-AXIS LABEL MANAGEMENT
  # ==========================================================================

  point_data$target_group_display <-
    .apply_moplot_label_key(
      labels = point_data$target_group,
      label_key = label_key
    )

  x_label_angle <- 45
  x_label_hjust <- 1
  x_label_vjust <- 1

  # ==========================================================================
  # 3. BASIC VALIDATION
  # ==========================================================================

  if (
    is.null(hierarchy_depth) ||
    length(hierarchy_depth) != 1L ||
    is.na(hierarchy_depth) ||
    hierarchy_depth < 2L
  ) {
    stop(
      "`prep$hierarchy_depth` must be a single integer of at least 2.",
      call. = FALSE
    )
  }

  required_delta_columns <- c(
    "coefficient",
    "contrast_id",
    "node_path",
    "node_level",
    "x",
    "y",
    "label",
    "positive_child",
    "negative_child",
    "is_final"
  )

  missing_delta_columns <- setdiff(required_delta_columns,names(delta_data))

  if (length(missing_delta_columns) > 0L) {
    stop(
      paste0(
        "`prep$delta_data` is missing the required column",
        if (length(missing_delta_columns) > 1L) "s: " else ": ",
        paste(
          missing_delta_columns,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  # ==========================================================================
  # 4. COLOUR HELPER
  # ==========================================================================

  lighten_colour <- function(
    colour,
    amount
  ) {

    rgb_value <- grDevices::col2rgb(colour) / 255

    rgb_value <- rgb_value + (1 - rgb_value) * amount

    grDevices::rgb(
      rgb_value[1L, ],
      rgb_value[2L, ],
      rgb_value[3L, ]
    )
  }


  # ==========================================================================
  # 5. STYLE SETTINGS
  # ==========================================================================

  grayscale_palette <- c(
    estimate = "#737373",
    segment = "#525252",
    midpoint = "#252525",
    guide = "#BDBDBD",
    text = "#252525"
  )

  coloured_palette <- c(
    estimate = "#737373",
    guide = "#BDBDBD",
    text = "#252525"
  )

  estimate_point_size <- if (
    identical(
      style,
      "minimal"
    )
  ) {
    2.8
  } else {
    3.5
  }

  estimate_linewidth <- if (
    identical(
      style,
      "minimal"
    )
  ) {
    0.60
  } else {
    0.80
  }

  segment_linewidth <- if (
    identical(
      style,
      "minimal"
    )
  ) {
    0.75
  } else {
    1.00
  }

  delta_label_size <- if (
    identical(
      style,
      "minimal"
    )
  ) {
    base_size / 3.8
  } else {
    base_size / 3.5
  }

  # Formula labels scale with base_size, so that they grow with the rest
  # of the figure; the default base_size of 13 gives 3.5 and 3.8.
  contrast_label_size <- if (
    identical(
      style,
      "minimal"
    )
  ) {
    base_size * 3.5 / 13
  } else {
    base_size * 3.8 / 13
  }

  if (identical(style, "minimal")) {

    palette <- grayscale_palette

    segment_data$Contrast_colour <- palette[["segment"]]

    delta_data$Base_colour <- palette[["midpoint"]]

  } else {

    palette <- coloured_palette

    contrast_ids <- sort(unique(segment_data$contrast_id))

    contrast_colours <- stats::setNames(
      .moplot_contrast_palette(
        length(
          contrast_ids
        )
      ),
      as.character(
        contrast_ids
      )
    )

    segment_data$Contrast_colour <- unname(
      contrast_colours[
        as.character(
          segment_data$contrast_id
        )
      ]
    )

    delta_data$Base_colour <- unname(
      contrast_colours[
        as.character(
          delta_data$contrast_id
        )
      ]
    )
  }

  # ==========================================================================
  # 6. HIERARCHICAL COLOUR AND SIZE PROGRESSION
  # ==========================================================================

  delta_data$Level_fraction <- pmin(
    pmax(
      (
        delta_data$node_level -
          2
      ) /
        max(
          hierarchy_depth -
            2,
          1
        ),
      0
    ),
    1
  )

  # Lower levels of the hierarchy are drawn in lighter shades of the
  # contrast colour. The lightening is limited so that their formulas
  # remain legible in print.
  delta_data$Lighten_amount <-
    0.25 -
    0.15 *
    delta_data$Level_fraction

  delta_data$Contrast_colour <- vapply(
    seq_len(
      nrow(
        delta_data
      )
    ),
    function(i) {

      if (
        isTRUE(
          delta_data$is_final[i]
        )
      ) {
        return(
          delta_data$Base_colour[i]
        )
      }

      lighten_colour(
        colour =
          delta_data$Base_colour[i],
        amount =
          delta_data$Lighten_amount[i]
      )
    },
    character(1)
  )

  delta_data$Point_size <-
    4.0 +
    1.5 *
    delta_data$Level_fraction

  delta_data$Point_size[
    delta_data$is_final
  ] <- 5.5

  # ==========================================================================
  # 7. Y-AXIS LIMITS AND BREAKS
  # ==========================================================================

  ci_min <- min(point_data$ci_low,na.rm = TRUE)

  ci_max <- max(point_data$ci_high,na.rm = TRUE)

  y_range <- c(ci_min,ci_max)

  y_span <- diff(y_range)

  if (
    !is.finite(y_span) ||
    y_span == 0
  ) {
    y_span <- max(
      abs(
        y_range
      ),
      1
    )
  }

  candidate_breaks <- pretty(y_range,n = 5)

  break_step <- stats::median(diff(candidate_breaks))

  if (
    !is.finite(break_step) ||
    break_step <= 0
  ) {
    break_step <- 1
  }

  lower_candidates <- candidate_breaks[candidate_breaks < ci_min]

  y_lower <- if (
    length(
      lower_candidates
    ) > 0L
  ) {
    max(
      lower_candidates
    )
  } else {
    min(
      candidate_breaks
    ) -
      break_step
  }

  upper_candidates <- candidate_breaks[candidate_breaks > ci_max]

  y_upper <- if (
    length(
      upper_candidates
    ) > 0L
  ) {
    min(
      upper_candidates
    )
  } else {
    max(
      candidate_breaks
    ) +
      break_step
  }

  y_breaks <- seq(
    from = y_lower,
    to = y_upper,
    by = break_step
  )

  # seq() can leave a residue such as -5.55e-17 in place of zero, which
  # makes the default labels switch the whole axis to scientific notation.
  y_breaks[abs(y_breaks) < 1e-8 * break_step] <- 0

  y_plot_lower <- y_lower - 0.10 * break_step

  # ==========================================================================
  # 8. HIERARCHICAL GROUP LABEL POSITIONS
  # ==========================================================================

  # Stacked labels (three-way and higher interactions) need a full line of
  # text between rows; a single row keeps the closer position.
  group_label_step <- if (
    nrow(group_label_data) > 0L &&
    max(group_label_data$display_level, na.rm = TRUE) > 1L
  ) {
    0.10 * y_span
  } else {
    0.075 * y_span
  }

  if (
    nrow(
      group_label_data
    ) > 0L
  ) {

    group_label_data$label_y <-
      ci_max +
      group_label_data$display_level *
      group_label_step

    group_label_data$display_label <-
      .apply_moplot_label_key(
        labels =
          group_label_data$group_value,
        label_key =
          label_key
      )

    highest_group_label_y <- max(
      group_label_data$label_y,
      na.rm = TRUE
    )

  } else {

    highest_group_label_y <- ci_max
  }

  y_plot_upper <- max(
    y_upper +
      0.10 *
      break_step,
    highest_group_label_y +
      0.06 *
      y_span
  )
   # ==========================================================================
  # 10. DATA FOR TREE STYLE
  # ==========================================================================

  tree_node_data <- NULL
  tree_segment_data <- NULL

  if (
    identical(
      style,
      "tree"
    )
  ) {

    tree_node_data <- delta_data

    # ------------------------------------------------------------------------
    # 10.1. VERTICAL SPACING
    # ------------------------------------------------------------------------
    #
    # Each hierarchical level occupies its own vertical band:
    #
    # incoming branch
    #       |
    #   Entry_y
    #
    #   formula
    #
    #    Tree_y
    #       |
    # outgoing branch
    #
    # This prevents formulas from overlapping the horizontal connectors and
    # ensures that all subsequent levels are shifted downward automatically.
    # ------------------------------------------------------------------------

    # The spacing below is computed for a given height of one formula line,
    # in data units; see 10.4b for how that height is chosen.
    layout_tree <- function(formula_line) {

      tree_top_gap <- max(
        0.075 *
          y_span,
        0.22 *
          break_step
      )

      tree_level_step <- max(
        0.060 *
          y_span,
        0.22 *
          break_step
      )

      tree_branch_drop <-
        0.22 *
        tree_level_step

      tree_central_drop <- max(
        0.010 *
          y_span,
        0.10 *
          tree_level_step
      )

      tree_formula_gap <- max(
        0.018 * y_span,
        0.12 * tree_level_step
      )

      tree_outgoing_gap <- max(
        0.035 * y_span,
        0.25 * tree_level_step
      )

      tree_top_formula_offset <-
        0.025 * y_span

      tree_final_x_gap <-
        0.065 * y_span

      # Every formula needs the height of one formula line: above the
      # outgoing branch (first level) and between the incoming and outgoing
      # branches (higher levels). A final formula is centred beside its X,
      # whose arms reach beyond half a line, so the X keeps a whole line
      # below the bridge of its incoming branch.
      tree_top_gap <- max(tree_top_gap, tree_top_formula_offset + formula_line)

      tree_outgoing_gap <- max(tree_outgoing_gap, formula_line)

      tree_final_x_gap <- max(tree_final_x_gap, formula_line - tree_central_drop)

      tree_top_y <-
        ci_min -
        tree_top_gap

      # ------------------------------------------------------------------------
      # 10.2. INITIALISE TREE COORDINATES
      # ------------------------------------------------------------------------

      tree_node_data$Entry_y <- NA_real_

      tree_node_data$Tree_y <- NA_real_

      tree_node_data$Formula_y <- NA_real_

      tree_node_data$Final_x_y <- NA_real_

      # Vertical band free for each formula, between the tree elements above
      # and below it. The formulas shrink to fit it when the panel is short
      # (see fit_text.R). The lower band of final formulas is the axis,
      # set once the lower limit is known (section 10.6).
      tree_node_data$Band_low <- NA_real_

      tree_node_data$Band_high <- NA_real_

      tree_levels <- sort(unique(tree_node_data$node_level))

      # ------------------------------------------------------------------------
      # 10.3. FIRST TREE LEVEL
      # ------------------------------------------------------------------------
      #
      # The first hierarchical formulas (e.g., Δ9-Δ12 in a four-way
      # interaction) retain the layout already used successfully: the formula
      # sits just above the point from which its branch descends.
      # ------------------------------------------------------------------------

      first_tree_level <- min(tree_levels)

      first_level_rows <- tree_node_data$node_level == first_tree_level

      first_level_final_rows <- first_level_rows & tree_node_data$is_final

      first_level_nonfinal_rows <- first_level_rows & !tree_node_data$is_final

      # For interactions of order >= 3, the first tree level contains
      # intermediate contrasts. Their formulas remain above the outgoing branches.
      tree_node_data$Tree_y[first_level_nonfinal_rows] <- tree_top_y

      tree_node_data$Formula_y[first_level_nonfinal_rows] <- tree_top_y + tree_top_formula_offset

      tree_node_data$Entry_y[first_level_nonfinal_rows] <- tree_node_data$Formula_y[first_level_nonfinal_rows]

      # For a two-way interaction, the first tree level is already the final
      # contrast. There is no intermediate branching structure, so the final
      # X is placed below the estimates, with its formula on the same line.
      tree_node_data$Entry_y[first_level_final_rows] <- tree_top_y

      tree_node_data$Formula_y[first_level_final_rows] <- tree_top_y - tree_final_x_gap

      tree_node_data$Tree_y[first_level_final_rows] <- tree_node_data$Formula_y[first_level_final_rows]

      tree_node_data$Final_x_y[first_level_final_rows] <- tree_node_data$Formula_y[first_level_final_rows]

      # First-level formulas have only the estimates above them: an
      # intermediate formula sits on its outgoing branch, a final formula
      # beside its X.
      tree_node_data$Band_low[first_level_nonfinal_rows] <- tree_node_data$Formula_y[first_level_nonfinal_rows]

      tree_node_data$Band_high[first_level_rows] <- ci_min

      # ------------------------------------------------------------------------
      # 10.4. POSITION SUCCESSIVE TREE LEVELS
      # ------------------------------------------------------------------------
      #
      # Every higher-order level is positioned recursively from the outgoing
      # branch position of the level immediately above it.
      # ------------------------------------------------------------------------

      if (
        length(
          tree_levels
        ) > 1L
      ) {

        for (
          current_level in tree_levels[-1L]
        ) {

          current_rows <- tree_node_data$node_level == current_level

          previous_level <-
            max(
              tree_levels[
                tree_levels <
                  current_level
              ]
            )

          previous_rows <- tree_node_data$node_level == previous_level

          previous_branch_y <- min(tree_node_data$Tree_y[previous_rows],na.rm = TRUE)

          # Horizontal bridge connecting the two child nodes.
          current_joint_y <- previous_branch_y - tree_branch_drop

          # End of the short central vertical segment.
          current_entry_y <- current_joint_y - tree_central_drop

          # Formula is always below the incoming central branch.
          current_formula_y <- current_entry_y - tree_formula_gap

          tree_node_data$Entry_y[current_rows] <- current_entry_y

          tree_node_data$Formula_y[current_rows] <- current_formula_y

          # For non-final nodes, subsequent branches begin below the formula.
          # This automatically pushes all later hierarchical levels downward.
          current_final_rows <- current_rows & tree_node_data$is_final

          current_nonfinal_rows <- current_rows & !tree_node_data$is_final

          tree_node_data$Tree_y[current_nonfinal_rows] <- current_formula_y - tree_outgoing_gap

          # The final node does not need an outgoing branch. Its X and its
          # formula share one line below the incoming branch, so that their
          # separation is horizontal and does not depend on the height of
          # the plot. Tree_y is retained as a valid internal coordinate.
          final_line_y <- current_entry_y - tree_final_x_gap

          tree_node_data$Formula_y[current_final_rows] <- final_line_y

          tree_node_data$Tree_y[current_final_rows] <- final_line_y

          tree_node_data$Final_x_y[current_final_rows] <- final_line_y

          # An intermediate formula hangs between its incoming branch and its
          # outgoing branch; a final formula, beside its X, must stay below
          # the horizontal bridge of its incoming branch.
          tree_node_data$Band_low[current_nonfinal_rows] <- tree_node_data$Tree_y[current_nonfinal_rows]

          tree_node_data$Band_high[current_nonfinal_rows] <- current_formula_y

          tree_node_data$Band_high[current_final_rows] <- current_joint_y
        }
      }


      list(
        nodes = tree_node_data,
        branch_drop = tree_branch_drop
      )
    }

    # ------------------------------------------------------------------------
    # 10.4b. HEIGHT OF ONE FORMULA LINE
    # ------------------------------------------------------------------------
    #
    # Formulas have a fixed size in points, so the room they need is a
    # fraction of the panel height, i.e. of the plotted range, which itself
    # grows with the tree. Iterating line = fraction * range converges
    # quickly because the tree takes up a small share of the range. At the
    # recommended export size a formula then fits at full size; on smaller
    # plots it shrinks to its band (see fit_text.R).
    # ------------------------------------------------------------------------

    formula_line_fraction <- 0.05

    tree_bottom_fraction <- 0.07

    # The final X and its formula keep clear of the x axis: the lowest tree
    # element lies at least tree_bottom_fraction of the plotted range above
    # it, and a final formula keeps half a line below its centre.
    tree_plot_lower <- function(nodes, formula_line) {

      lowest_tree_y <- min(
        c(
          nodes$Tree_y,
          nodes$Entry_y,
          nodes$Formula_y,
          nodes$Final_x_y,
          nodes$Final_x_y[nodes$is_final] - formula_line / 2
        ),
        na.rm = TRUE
      )

      min(
        y_plot_lower,
        (lowest_tree_y - tree_bottom_fraction * y_plot_upper) / (1 - tree_bottom_fraction)
      )
    }

    formula_line <- 0

    for (iteration in seq_len(50L)) {

      tree_layout <- layout_tree(formula_line)

      plotted_range <- y_plot_upper - tree_plot_lower(tree_layout$nodes, formula_line)

      next_formula_line <- formula_line_fraction * plotted_range

      if (abs(next_formula_line - formula_line) <= 1e-6 * plotted_range) {
        break
      }

      formula_line <- next_formula_line
    }

    tree_node_data <- tree_layout$nodes

    tree_branch_drop <- tree_layout$branch_drop

    first_tree_level <- min(tree_node_data$node_level)

    # ------------------------------------------------------------------------
    # 10.5. BUILD TREE CONNECTORS
    # ------------------------------------------------------------------------

    connector_parts <- list()
    connector_counter <- 1L

    parent_rows <- tree_node_data[tree_node_data$node_level > first_tree_level, ,drop = FALSE]

    if (
      nrow(
        parent_rows
      ) > 0L
    ) {

      for (
        i in seq_len(
          nrow(
            parent_rows
          )
        )
      ) {

        parent <- parent_rows[
          i,
          ,
          drop = FALSE
        ]

        children <- tree_node_data[
          tree_node_data$coefficient ==
            parent$coefficient &
            tree_node_data$node_path %in%
            c(
              parent$positive_child,
              parent$negative_child
            ),
          ,
          drop = FALSE
        ]

        if (
          nrow(
            children
          ) != 2L
        ) {
          next
        }

        children <- children[
          match(
            c(
              parent$positive_child,
              parent$negative_child
            ),
            children$node_path
          ),
          ,
          drop = FALSE
        ]

        child_start_y <- children$Tree_y

        # Both children belong to the same hierarchical level and therefore
        # normally share the same y coordinate. mean() keeps the construction
        # robust if tiny numerical differences are ever introduced.
        joint_y <- mean(
          child_start_y
        ) -
          tree_branch_drop

        connector_colour <- parent$Base_colour[1L]

        # --------------------------------------------------------------------
        # Left child vertical branch
        # --------------------------------------------------------------------

        connector_parts[[
          connector_counter
        ]] <- data.frame(
          contrast_id =
            parent$contrast_id,
          x =
            children$x[1L],
          y =
            child_start_y[1L],
          xend =
            children$x[1L],
          yend =
            joint_y,
          Contrast_colour =
            connector_colour,
          stringsAsFactors =
            FALSE
        )

        connector_counter <- connector_counter + 1L

        # --------------------------------------------------------------------
        # Right child vertical branch
        # --------------------------------------------------------------------

        connector_parts[[
          connector_counter
        ]] <- data.frame(
          contrast_id =
            parent$contrast_id,
          x =
            children$x[2L],
          y =
            child_start_y[2L],
          xend =
            children$x[2L],
          yend =
            joint_y,
          Contrast_colour =
            connector_colour,
          stringsAsFactors =
            FALSE
        )

        connector_counter <-
          connector_counter +
          1L

        # --------------------------------------------------------------------
        # Horizontal bridge between children
        # --------------------------------------------------------------------

        connector_parts[[
          connector_counter
        ]] <- data.frame(
          contrast_id =
            parent$contrast_id,
          x =
            children$x[1L],
          y =
            joint_y,
          xend =
            children$x[2L],
          yend =
            joint_y,
          Contrast_colour =
            connector_colour,
          stringsAsFactors =
            FALSE
        )

        connector_counter <- connector_counter + 1L

        # --------------------------------------------------------------------
        # Short central descent
        #
        # Crucially, this stops at Entry_y, ABOVE the parent formula.
        # It no longer crosses the formula itself.
        # --------------------------------------------------------------------

        connector_parts[[
          connector_counter
        ]] <- data.frame(
          contrast_id =
            parent$contrast_id,
          x =
            parent$x[1L],
          y =
            joint_y,
          xend =
            parent$x[1L],
          yend =
            parent$Entry_y[1L],
          Contrast_colour =
            connector_colour,
          stringsAsFactors =
            FALSE
        )

        connector_counter <-
          connector_counter +
          1L
      }
    }

    if (
      length(
        connector_parts
      ) > 0L
    ) {

      tree_segment_data <- do.call(rbind,connector_parts)
    }

    # ------------------------------------------------------------------------
    # 10.6. LOWER PLOT LIMIT
    # ------------------------------------------------------------------------

    # The lower limit keeps the final X and its formula clear of the x axis
    # (see tree_plot_lower() in 10.4b). Their size is fixed in points, so
    # the space below them is a fraction of the plotted range, i.e. of the
    # panel height, rather than of the estimates.
    y_plot_lower <- tree_plot_lower(tree_node_data, formula_line)

    tree_node_data$Band_low[tree_node_data$is_final] <- y_plot_lower
  }


  # ==========================================================================
  # 12. CREATE PLOT
  # ==========================================================================

  plot <- ggplot2::ggplot()

  # ==========================================================================
  # 13. ELEMENTARY CONTRAST SEGMENTS
  # ==========================================================================

  plot <- plot +
    ggplot2::geom_segment(
      data =
        segment_data,
      ggplot2::aes(
        x = x,
        xend = xend,
        y = y,
        yend = yend,
        colour =
          Contrast_colour
      ),
      linetype =
        "dashed",
      linewidth =
        segment_linewidth
    )

  # ==========================================================================
  # 14. ELEMENTARY DELTA LABELS
  # ==========================================================================

  if (
    !identical(
      style,
      "minimal"
    )
  ) {

    plot <- plot +
      ggplot2::geom_label(
        data =
          segment_data,
        ggplot2::aes(
          x =
            delta_x,
          y =
            delta_y,
          label =
            delta_label,
          colour =
            Contrast_colour
        ),
        inherit.aes =
          FALSE,
        parse=TRUE,
        size =
          delta_label_size,
        fill =
          "white",
        linewidth =
          0.2,
        label.r =
          grid::unit(
            0.10,
            "lines"
          )
      )
  }

  # ==========================================================================
  # 15. ESTIMATES AND CONFIDENCE INTERVALS
  # ==========================================================================

  plot <- plot +
    ggplot2::geom_errorbar(
      data =
        point_data,
      ggplot2::aes(
        x = x,
        ymin =
          ci_low,
        ymax =
          ci_high
      ),
      width =
        0.10,
      linewidth =
        estimate_linewidth,
      colour =
        palette[["estimate"]]
    )

  plot <- plot +
    ggplot2::geom_point(
      data =
        point_data,
      ggplot2::aes(
        x = x,
        y =
          estimate
      ),
      size =
        estimate_point_size,
      colour =
        palette[["estimate"]]
    )

  # ==========================================================================
  # 16. HIERARCHICAL GROUP LABELS
  # ==========================================================================

  if (
    nrow(
      group_label_data
    ) > 0L
  ) {

    # Neighbouring labels shrink when the plot is too narrow for them;
    # their size is fitted whenever the plot is drawn (see fit_text.R).
    plot <- plot +
      ggplot2::layer(
        geom =
          GeomMoplotFitText,
        stat =
          "identity",
        position =
          "identity",
        data =
          group_label_data,
        mapping =
          ggplot2::aes(
            x = x,
            y =
              label_y,
            label =
              display_label
          ),
        inherit.aes =
          FALSE,
        show.legend =
          FALSE,
        params =
          list(
            fontface =
              "bold",
            colour =
              palette[["text"]],
            vjust =
              0
          )
      )
  }



  # ==========================================================================
  # 19. TREE STYLE
  # ==========================================================================

  if (
    identical(
      style,
      "tree"
    )
  ) {

    # ------------------------------------------------------------------------
    # 19.1. TREE CONNECTORS
    # ------------------------------------------------------------------------

    if (
      !is.null(
        tree_segment_data
      ) &&
      nrow(
        tree_segment_data
      ) > 0L
    ) {

      plot <- plot +
        ggplot2::geom_segment(
          data =
            tree_segment_data,
          ggplot2::aes(
            x = x,
            xend = xend,
            y = y,
            yend = yend,
            colour =
              Contrast_colour
          ),
          inherit.aes =
            FALSE,
          linewidth =
            0.50
        )
    }

    # ------------------------------------------------------------------------
    # 19.2. TREE FORMULAS
    #
    # All formulas are drawn by one layer, so that they share one text size
    # when they have to shrink to fit their band on a short or narrow panel.
    #
    # - First hierarchical level: the formula sits above its outgoing branch.
    # - Higher intermediate levels: incoming branch -> formula -> outgoing
    #   branch.
    # - Final contrast: the X marks it at the centre of its block, aligned
    #   with the coefficient in the lower panel, and the formula is written
    #   on the same line, to its right. A negative hjust shifts the text by a
    #   fraction of its own width, so the gap is a physical distance:
    #   stacking the formula and the X vertically required a gap in data
    #   units, which overlapped on short panels or with a small range of the
    #   data.
    # ------------------------------------------------------------------------

    first_tree_level <- min(tree_node_data$node_level)

    tree_node_data$Formula_hjust <- ifelse(tree_node_data$is_final, -0.25, 0.5)

    tree_node_data$Formula_vjust <- ifelse(
      tree_node_data$is_final,
      0.5,
      ifelse(tree_node_data$node_level == first_tree_level, 0, 1)
    )

    plot <- plot +
      ggplot2::layer(
        geom =
          GeomMoplotFitText,
        stat =
          "identity",
        position =
          "identity",
        data =
          tree_node_data,
        mapping =
          ggplot2::aes(
            x = x,
            y = Formula_y,
            ymin = Band_low,
            ymax = Band_high,
            label = label,
            colour =
              Contrast_colour,
            hjust =
              Formula_hjust,
            vjust =
              Formula_vjust
          ),
        inherit.aes =
          FALSE,
        show.legend =
          FALSE,
        params =
          list(
            parse =
              TRUE,
            size =
              contrast_label_size
          )
      )

    # ------------------------------------------------------------------------
    # 19.3. FINAL TREE X
    # ------------------------------------------------------------------------

    tree_final_data <- tree_node_data[
      tree_node_data$is_final,
      ,
      drop = FALSE
    ]

    if (
      nrow(
        tree_final_data
      ) > 0L
    ) {

      plot <- plot +
        ggplot2::geom_point(
          data =
            tree_final_data,
          ggplot2::aes(
            x = x,
            y = Final_x_y,
            colour =
              Contrast_colour
          ),
          inherit.aes =
            FALSE,
          shape =
            4,
          size =
            5.5,
          stroke =
            1.2,
          show.legend =
            FALSE
        )
    }
  }

  # ==========================================================================
  # 20. IDENTITY SCALES
  # ==========================================================================

  plot <- plot + ggplot2::scale_colour_identity()


  # ==========================================================================
  # 21. AXIS SCALES
  # ==========================================================================

  plot <- plot +
    ggplot2::scale_x_continuous(
      breaks =
        point_data$x,
      labels =
        point_data$target_group_display,
      expand =
        ggplot2::expansion(
          mult = 0,
          add = 0
        )
    )

  plot <- plot +
    ggplot2::scale_y_continuous(
      breaks =
        y_breaks,
      expand =
        ggplot2::expansion(
          mult = 0,
          add = 0
        )
    )

  # ==========================================================================
  # 22. AXIS LABELS
  # ==========================================================================

  plot <- plot +
    ggplot2::labs(
      x =
        prep$target_var,
      y =
        y_title
    )

  # ==========================================================================
  # 23. COORDINATES
  # ==========================================================================

  plot <- plot +
    ggplot2::coord_cartesian(
      xlim =
        prep$x_limits,
      ylim =
        c(
          y_plot_lower,
          y_plot_upper
        ),
      clip =
        "off",
      expand =
        FALSE
    )

  # ==========================================================================
  # 24. THEME
  # ==========================================================================

  plot <- plot +
    ggplot2::theme_minimal(
      base_size =
        base_size
    )

  plot <- plot +
    ggplot2::theme(
      panel.grid.major.x =
        ggplot2::element_blank(),

      panel.grid.minor =
        ggplot2::element_blank(),

      axis.line.x =
        ggplot2::element_line(
          linewidth =
            0.4,
          colour =
            palette[["text"]]
        ),

      axis.line.y =
        ggplot2::element_line(
          linewidth =
            0.4,
          colour =
            palette[["text"]]
        ),

      axis.text.x =
        ggplot2::element_text(
          angle =
            x_label_angle,
          hjust =
            x_label_hjust,
          vjust =
            x_label_vjust,
          colour =
            palette[["text"]]
        ),

      axis.text.y =
        ggplot2::element_text(
          colour =
            palette[["text"]]
        ),

      axis.title =
        ggplot2::element_text(
          face = if (
            style %in%
            c(
              "tree"
            )
          ) {
            "bold"
          } else {
            "plain"
          }
        ),

      plot.margin =
        ggplot2::margin(
          t = 18,
          r = 10,
          b = 10,
          l = 10
        )
    )

  plot
}


#' Prepare model parameters for a multi-predictor MoPlot
#'
#' Internal helper that extracts the model parameters corresponding to the
#' interaction coefficients displayed by `moplot_multpred()`.
#'
#' @param model A fitted model supported by `extract_param()`.
#' @param coefficient_data Coefficient-layout data returned by
#'   `.prepare_moplot_multpred_data()`.
#' @param conf_level Confidence level used for parameter confidence intervals.
#' @param coefficient_scale Scale used for the lower-panel estimates.
#'   Available options are `"model"` and `"interpreted"`.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `parameter_data`: plot-ready coefficient information;
#'     \item `parameter_table`: complete `extract_param()` output for the
#'       interaction coefficients;
#'     \item `null_value`: null reference value for the displayed coefficient
#'       scale;
#'     \item `coefficient_scale`: coefficient scale used.
#'   }
#'
#' @keywords internal
#' @noRd
.prepare_moplot_multpred_parameters <- function(
    model,
    coefficient_data,
    conf_level = 0.95,
    coefficient_scale = c("model", "interpreted"),
    scale_factors = NULL,
    per_difference = FALSE
) {

  # ============================================================
  # 1. VALIDATION
  # ============================================================

  coefficient_scale <- match.arg(
    coefficient_scale
  )

  required_columns <- c(
    "coefficient",
    "contrast_id",
    "significance",
    "p_value",
    "coefficient_index",
    "x"
  )

  missing_columns <- setdiff(
    required_columns,
    names(coefficient_data)
  )

  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        "`coefficient_data` is missing the required column",
        if (length(missing_columns) > 1L) "s: " else ": ",
        paste(missing_columns, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  # ============================================================
  # 2. EXTRACT MODEL PARAMETERS
  # ============================================================

  parameter_table <- extract_param(
    model = model,
    conf_level = conf_level,
    include_intercept = FALSE
  )

  parameter_indices <- match(coefficient_data$coefficient,parameter_table$Parameter)

  if (anyNA(parameter_indices)) {

    unmatched_parameters <- coefficient_data$coefficient[is.na(parameter_indices)]

    stop(
      paste0(
        "The following interaction coefficients could not be matched to ",
        "`extract_param()` output: ",
        paste(
          unique(unmatched_parameters),
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  interaction_parameters <- parameter_table[parameter_indices, ,drop = FALSE]

  ratio_types <- c(
    "Odds ratio",
    "Risk ratio",
    "Incidence rate ratio",
    "Hazard ratio",
    "Multiplicative ratio",
    "Exponentiated coefficient"
  )

  # Cohen's d converted from the log odds ratio refers to the
  # coefficient as coded; dividing by the scale factor makes it refer to
  # the displayed interaction contrast. Conversions from t or z are
  # already independent of the scaling of the contrasts.
  d_from_log_odds <- interaction_parameters$d_method %in% "log odds ratio"

  if (any(d_from_log_odds) &&
      !is.null(scale_factors)) {

    d_scales <- unname(scale_factors[coefficient_data$coefficient])[d_from_log_odds]

    interaction_parameters$d[d_from_log_odds] <- interaction_parameters$d[d_from_log_odds] / d_scales
    interaction_parameters$d_CI_low[d_from_log_odds] <- interaction_parameters$d_CI_low[d_from_log_odds] / d_scales
    interaction_parameters$d_CI_high[d_from_log_odds] <- interaction_parameters$d_CI_high[d_from_log_odds] / d_scales
  }

  # Each interaction coefficient equals its scale factor multiplied by
  # the displayed difference of differences on the linear-predictor
  # scale. Dividing by the scale factor expresses the coefficient per
  # unit of that difference; test statistics and p-values are unchanged.
  if (isTRUE(per_difference)) {

    if (is.null(scale_factors) ||
        anyNA(scale_factors[coefficient_data$coefficient])) {

      stop(
        "Scale factors are required when `per_difference = TRUE`.",
        call. = FALSE
      )
    }

    coefficient_scales <- unname(scale_factors[coefficient_data$coefficient])
    is_ratio <- interaction_parameters$Estimate_type %in% ratio_types

    interaction_parameters$Estimate <- interaction_parameters$Estimate / coefficient_scales
    interaction_parameters$CI_low <- interaction_parameters$CI_low / coefficient_scales
    interaction_parameters$CI_high <- interaction_parameters$CI_high / coefficient_scales

    interaction_parameters$Interpreted_Estimate <- ifelse(
      is_ratio,
      exp(interaction_parameters$Estimate),
      interaction_parameters$Interpreted_Estimate / coefficient_scales
    )

    interaction_parameters$Interpreted_CI_low <- ifelse(
      is_ratio,
      exp(interaction_parameters$CI_low),
      interaction_parameters$Interpreted_CI_low / coefficient_scales
    )

    interaction_parameters$Interpreted_CI_high <- ifelse(
      is_ratio,
      exp(interaction_parameters$CI_high),
      interaction_parameters$Interpreted_CI_high / coefficient_scales
    )
  }

  # ============================================================
  # 3. DISPLAY SCALE
  # ============================================================

  if (identical(coefficient_scale, "model")) {

    display_estimate <- interaction_parameters$Estimate

    display_ci_low <- interaction_parameters$CI_low

    display_ci_high <- interaction_parameters$CI_high

    display_type <- rep("Model estimate",nrow(interaction_parameters))

    null_value <- 0

  } else {

    display_estimate <- interaction_parameters$Interpreted_Estimate

    display_ci_low <- interaction_parameters$Interpreted_CI_low

    display_ci_high <- interaction_parameters$Interpreted_CI_high

    display_type <- interaction_parameters$Estimate_type

    null_value <- if (
      all(interaction_parameters$Estimate_type %in% ratio_types)
    ) {
      1
    } else {
      0
    }
  }

  # ============================================================
  # 4. BUILD PLOT-READY PARAMETER DATA
  # ============================================================

  parameter_data <- data.frame(
    Contrast = coefficient_data$contrast_id,
    Contrast_number = coefficient_data$coefficient_index,
    Parameter = interaction_parameters$Parameter,
    x = coefficient_data$x,
    Display_estimate = display_estimate,
    Display_CI_low = display_ci_low,
    Display_CI_high = display_ci_high,
    Display_type = display_type,
    d = interaction_parameters$d,
    d_CI_low = interaction_parameters$d_CI_low,
    d_CI_high = interaction_parameters$d_CI_high,
    significance = coefficient_data$significance,
    p_value = coefficient_data$p_value,
    stringsAsFactors = FALSE
  )

  # ============================================================
  # 5. RETURN
  # ============================================================

  list(
    parameter_data = parameter_data,
    parameter_table = interaction_parameters,
    null_value = null_value,
    coefficient_scale = coefficient_scale,
    per_difference = isTRUE(per_difference)
  )
}



#' Build the lower panel of an n-way multi-predictor contrast plot
#'
#' Experimental internal helper used with
#' `.prepare_moplot_multpred_data()`.
#'
#' @param prep An object returned by `.prepare_moplot_multpred_data()`.
#'
#' @param parameters An object returned by
#'   `.prepare_moplot_multpred_parameters()`.
#'
#' @param style Character string specifying the graphical style. Available
#'   options are `"minimal"` and `"tree"`.
#'
#' @param digits Number of decimal places used in estimate and effect-size
#'   labels. The default is `2`.
#'
#' @param base_size Base font size passed to `ggplot2::theme_light()`.
#'   The default is `13`.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd
.plot_moplot_multpred_lower <- function(
    prep,
    parameters,
    style = c("minimal", "tree"),
    digits = 2L,
    base_size = 13
) {

  # ==========================================================================
  # 1. ARGUMENT MATCHING AND VALIDATION
  # ==========================================================================

  style <- match.arg(
    style
  )

  if (
    !is.list(prep) ||
    is.null(prep$point_data) ||
    is.null(prep$coefficient_data)
  ) {
    stop(
      paste0(
        "`prep` must be an object returned by ",
        "`.prepare_moplot_multpred_data()`."
      ),
      call. = FALSE
    )
  }

  if (
    !is.list(parameters) ||
    is.null(parameters$parameter_data)
  ) {
    stop(
      paste0(
        "`parameters` must be an object returned by ",
        "`.prepare_moplot_multpred_parameters()`."
      ),
      call. = FALSE
    )
  }

  if (
    !is.numeric(digits) ||
    length(digits) != 1L ||
    is.na(digits) ||
    !is.finite(digits) ||
    digits < 0 ||
    digits != as.integer(digits)
  ) {
    stop(
      "`digits` must be a single non-negative integer.",
      call. = FALSE
    )
  }

  digits <- as.integer(digits)

  if (
    !is.numeric(base_size) ||
    length(base_size) != 1L ||
    is.na(base_size) ||
    !is.finite(base_size) ||
    base_size <= 0
  ) {
    stop(
      "`base_size` must be a single positive numeric value.",
      call. = FALSE
    )
  }

  parameter_data <- parameters$parameter_data
  coefficient_data <- prep$coefficient_data

  required_parameter_columns <- c(
    "Contrast",
    "Contrast_number",
    "Parameter",
    "x",
    "Display_estimate",
    "Display_CI_low",
    "Display_CI_high",
    "Display_type",
    "d",
    "significance"
  )

  missing_parameter_columns <- setdiff(required_parameter_columns,names(parameter_data))

  if (length(missing_parameter_columns) > 0L) {
    stop(
      paste0(
        "`parameters$parameter_data` is missing the following columns: ",
        paste(
          missing_parameter_columns,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  parameter_data$x <- as.numeric(parameter_data$x)

  parameter_data$Display_estimate <- as.numeric(parameter_data$Display_estimate)

  parameter_data$Display_CI_low <- as.numeric(parameter_data$Display_CI_low)

  parameter_data$Display_CI_high <- as.numeric(parameter_data$Display_CI_high)

  parameter_data$d <- as.numeric(parameter_data$d)

  if (
    anyNA(parameter_data$x) ||
    any(!is.finite(parameter_data$x))
  ) {
    stop(
      "The horizontal positions of the model parameters are invalid.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # 2. VERIFY ALIGNMENT WITH UPPER-PANEL COEFFICIENT BLOCKS
  # ==========================================================================

  coefficient_midpoints <- as.numeric(coefficient_data$x)

  if (
    length(coefficient_midpoints) != nrow(parameter_data) ||
    !isTRUE(
      all.equal(
        parameter_data$x,
        coefficient_midpoints,
        tolerance = 1e-10
      )
    )
  ) {
    stop(
      paste0(
        "The parameter positions are not aligned with the coefficient ",
        "midpoints defined by the upper panel."
      ),
      call. = FALSE
    )
  }

  if (
    nrow(parameter_data) != nrow(coefficient_data) ||
    !identical(
      parameter_data$Parameter,
      coefficient_data$coefficient
    )
  ) {
    stop(
      paste0(
        "The lower-panel parameters could not be aligned with the ",
        "interaction coefficients represented in the upper panel."
      ),
      call. = FALSE
    )
  }

  # ==========================================================================
  # 3. ESTIMATE ABBREVIATIONS AND LABELS
  # ==========================================================================

  # Plotmath labels, shared with moplot_onepred(); see
  # .parameter_label_expressions().
  labels <- .parameter_label_expressions(
    display_type = parameter_data$Display_type,
    estimate = parameter_data$Display_estimate,
    d = parameter_data$d,
    per_difference = parameters$per_difference,
    digits = digits
  )

  parameter_data$Estimate_label <- labels$estimate_label
  parameter_data$Text_label <- labels$text_label

  significant <- !is.na(
    parameter_data$significance
  ) &
    parameter_data$significance == "Sig"

  parameter_data$Contrast_label <- paste0(
    "Coef. ",
    parameter_data$Contrast_number,
    ifelse(
      significant,
      "*",
      ""
    )
  )

  # ==========================================================================
  # 4. STYLE SETTINGS
  # ==========================================================================

  grayscale_palette <- c(
    midpoint = "#252525",
    interval = "#969696",
    guide = "#BDBDBD",
    text = "#252525",
    label_border = "#BDBDBD",
    label_fill = "white"
  )

  coloured_palette <- c(
    interval = "#969696",
    guide = "#BDBDBD",
    text = "#252525",
    label_border = "#BDBDBD",
    label_fill = "white"
  )

  if (identical(style, "minimal")) {

    palette <- grayscale_palette

    point_size <- 4.5
    point_stroke <- 1.0

    interval_linewidth <- 0.70
    guide_linewidth <- 0.40

    label_size <- 3.1
    label_border_width <- 0.20

    label_x_offset <- 0.16

    parameter_data$Point_colour <-
      palette[["midpoint"]]

  } else {

    palette <- coloured_palette

    point_size <- 5.5
    point_stroke <- 1.2

    interval_linewidth <- 0.85
    guide_linewidth <- 0.45

    label_size <- 3.3
    label_border_width <- 0.22

    label_x_offset <- 0.18

    contrast_ids <- sort(
      unique(
        parameter_data$Contrast_number
      )
    )

    contrast_colours <- stats::setNames(
      .moplot_contrast_palette(
        length(contrast_ids)
      ),
      as.character(contrast_ids)
    )

    parameter_data$Point_colour <- unname(
      contrast_colours[
        as.character(
          parameter_data$Contrast_number
        )
      ]
    )
  }

  # ==========================================================================
  # 5. AXIS LIMITS
  # ==========================================================================

  displayed_y <- c(
    parameter_data$Display_estimate,
    parameter_data$Display_CI_low,
    parameter_data$Display_CI_high
  )

  displayed_y <- displayed_y[
    is.finite(
      displayed_y
    )
  ]

  if (length(displayed_y) == 0L) {
    stop(
      "No finite values were available for the lower-panel y-axis.",
      call. = FALSE
    )
  }

  y_range <- range(displayed_y)

  y_span <- diff(y_range)

  if (
    !is.finite(y_span) ||
    y_span == 0
  ) {
    y_span <- max(
      abs(y_range),
      1
    )
  }

  y_lower <- y_range[1L] - 0.12 * y_span

  y_upper <- y_range[2L] + 0.12 * y_span

  # Use the full upper-panel horizontal range so that the upper and lower
  # panels remain perfectly aligned regardless of interaction order.
  x_limits <- prep$x_limits

  # ==========================================================================
  # 6. LABEL POSITIONS
  # ==========================================================================

  parameter_data$Label_x <- parameter_data$x + label_x_offset

  parameter_data$Label_y <- parameter_data$Display_estimate

  # ==========================================================================
  # 7. CREATE PLOT
  # ==========================================================================

  plot <- ggplot2::ggplot()

  # Vertical guides identify the centre of each coefficient block.
  plot <- plot +
    ggplot2::geom_vline(
      xintercept = coefficient_midpoints,
      linewidth = guide_linewidth,
      linetype = "dotted",
      colour = palette[["guide"]]
    )

  # Confidence intervals are drawn as plain vertical lines without
  # horizontal end caps.
  plot <- plot +
    ggplot2::geom_linerange(
      data = parameter_data,
      mapping = ggplot2::aes(
        x = x,
        ymin = Display_CI_low,
        ymax = Display_CI_high
      ),
      linewidth = interval_linewidth,
      colour = palette[["interval"]],
      inherit.aes = FALSE
    )

  # Interaction coefficients are represented by X markers.
  plot <- plot +
    ggplot2::geom_point(
      data = parameter_data,
      mapping = ggplot2::aes(
        x = x,
        y = Display_estimate,
        colour = Point_colour
      ),
      shape = 4,
      size = point_size,
      stroke = point_stroke,
      show.legend = FALSE,
      inherit.aes = FALSE
    )

  plot <- plot +
    ggplot2::scale_colour_identity()

  plot <- plot +
    ggplot2::geom_label(
      data = parameter_data,
      mapping = ggplot2::aes(
        x = Label_x,
        y = Label_y,
        label = Text_label
      ),
      parse = TRUE,
      hjust = 0,
      vjust = 0.5,
      size = label_size,
      fill = palette[["label_fill"]],
      colour = palette[["text"]],
      linewidth = label_border_width,
      label.r = grid::unit(
        0.12,
        "lines"
      ),
      show.legend = FALSE,
      inherit.aes = FALSE
    )

  # ==========================================================================
  # 8. SCALES AND THEME
  # ==========================================================================

  plot <- plot +
    ggplot2::scale_x_continuous(
      breaks = parameter_data$x,
      labels = parameter_data$Contrast_label,
      expand = ggplot2::expansion(
        mult = 0,
        add = 0
      )
    )

  plot <- plot +
    ggplot2::coord_cartesian(
      xlim = x_limits,
      ylim = c(
        y_lower,
        y_upper
      ),
      clip = "off",
      expand = FALSE
    )

  plot <- plot +
    ggplot2::labs(
      x = "Contrast",
      y = if (
        identical(
          parameters$coefficient_scale,
          "model"
        )
      ) {
        if (isTRUE(parameters$per_difference)) "Model estimate / s" else "Model estimate"
      } else {
        unique(
          parameter_data$Display_type
        )[1L]
      }
    )

  plot <- plot +
    ggplot2::theme_light(
      base_size = base_size
    )

  plot <- plot +
    ggplot2::theme(
      panel.grid.major.x =
        ggplot2::element_blank(),

      panel.grid.minor.x =
        ggplot2::element_blank(),

      panel.grid.major.y =
        ggplot2::element_line(
          colour = "grey92"
        ),

      panel.grid.minor.y =
        ggplot2::element_blank(),

      axis.title =
        ggplot2::element_text(
          face = if (
            style %in%
            c(
              "tree"
            )
          ) {
            "bold"
          } else {
            "plain"
          }
        ),

      axis.text.x =
        ggplot2::element_text(
          colour = palette[["text"]]
        ),

      axis.text.y =
        ggplot2::element_text(
          colour = palette[["text"]]
        ),

      plot.margin =
        ggplot2::margin(
          t = 8,
          r = 18,
          b = 10,
          l = 10
        )
    )

  plot
}

#' Assemble a multi-predictor MoPlot
#'
#' Internal helper that combines the upper contrast panel and lower
#' coefficient panel.
#'
#' @param upper_plot Upper-panel `ggplot` object.
#' @param lower_plot Lower-panel `ggplot` object.
#' @param panel_heights Numeric vector of length two specifying the relative
#'   heights of the upper and lower panels.
#'
#' @return A `patchwork` object.
#'
#' @keywords internal
#' @noRd
.assemble_moplot_multpred <- function(
    upper_plot,
    lower_plot,
    panel_heights = c(2, 1)
) {

  if (
    !is.numeric(panel_heights) ||
    length(panel_heights) != 2L ||
    anyNA(panel_heights) ||
    any(!is.finite(panel_heights)) ||
    any(panel_heights <= 0)
  ) {
    stop(
      "`panel_heights` must contain two positive finite numeric values.",
      call. = FALSE
    )
  }

  patchwork::wrap_plots(
    upper_plot,
    lower_plot,
    ncol = 1,
    heights = panel_heights
  )
}

#' Generate a textual description for a multi-predictor MoPlot
#'
#' Internal helper used by `moplot_multpred()` to generate a natural-language
#' description of an interaction-contrast plot of arbitrary order.
#'
#' @param model A fitted model supported by the package.
#' @param contrast_output A data frame returned by `contrasts_multpred()`.
#' @param prep An object returned by `.prepare_moplot_multpred_data()`.
#' @param interaction Character string identifying the interaction displayed
#'   in the plot.
#' @param style Character string specifying the graphical style. Available
#'   options are `"minimal"`, `"delta"`, and `"tree"`.
#' @param sig_level Significance threshold used for the contrast labels.
#' @param show_coef Logical. Whether the lower coefficient panel is displayed.
#' @param coefficient_scale Scale used for the lower-panel estimates.
#' @param estimate_type Character string identifying how the upper-panel
#'   quantities were computed.
#' @param prediction_scale Character string identifying the prediction scale.
#'
#' @return A single character string containing the plot description.
#'
#' @keywords internal
#' @noRd
.moplot_multpred_caption <- function(
    model,
    contrast_output,
    prep,
    interaction,
    style = c("minimal", "tree"),
    sig_level = 0.05,
    show_coef = TRUE,
    coefficient_scale = c("model", "interpreted"),
    per_difference = FALSE,
    estimate_type,
    prediction_scale,
    label_key= NULL
) {

  # ==========================================================================
  # 1. ARGUMENTS
  # ==========================================================================

  style <- match.arg(style)

  coefficient_scale <- match.arg(coefficient_scale)

  estimate_type <- tolower(as.character(estimate_type)[1L])

  prediction_scale <- tolower(as.character(prediction_scale)[1L])

  if (
    !is.character(interaction) ||
    length(interaction) != 1L ||
    is.na(interaction) ||
    !nzchar(interaction)
  ) {
    stop(
      "`interaction` must be a single non-empty character string.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # 2. MODEL AND INTERACTION INFORMATION
  # ==========================================================================

  model_information <- model_contrast_table(model,type = "all")

  variable_information <- model_information[["Variable-Level Summary"]]

  dependent_variables <- variable_information$Variable[variable_information$Role == "Dependent"]

  if (
    length(
      dependent_variables
    ) != 1L
  ) {
    stop(
      "`.moplot_multpred_caption()` requires exactly one dependent variable.",
      call. = FALSE
    )
  }

  dependent_variable <- dependent_variables[1L]

  interaction_variables <- strsplit(
    interaction,
    ":",
    fixed = TRUE
  )[[1L]]

  hierarchy_variables <- prep$variable_path

  hierarchy_depth <- prep$hierarchy_depth

  if (
    is.null(hierarchy_variables) ||
    length(hierarchy_variables) < 2L
  ) {
    stop(
      "The graphical interaction hierarchy could not be identified.",
      call. = FALSE
    )
  }

  if (
    is.null(hierarchy_depth) ||
    length(hierarchy_depth) != 1L ||
    is.na(hierarchy_depth)
  ) {
    hierarchy_depth <-
      length(
        hierarchy_variables
      )
  }

  split_variable <- hierarchy_variables[1L]

  target_variable <- hierarchy_variables[length(hierarchy_variables)]

  # ==========================================================================
  # 3. DISPLAYED ESTIMATE DESCRIPTION
  # ==========================================================================

  family_object <- tryCatch(
    stats::family(
      model
    ),
    error = function(e) NULL
  )

  family_name <- if (
    is.null(
      family_object
    )
  ) {
    "gaussian"
  } else {
    tolower(
      family_object$family
    )
  }

  link_name <- if (
    is.null(
      family_object
    )
  ) {
    "identity"
  } else {
    tolower(
      family_object$link
    )
  }

  mean_description <- if (
    identical(
      estimate_type,
      "observed"
    )
  ) {

    if (
      family_name %in%
      c(
        "binomial",
        "quasibinomial"
      )
    ) {

      "observed proportions"

    } else if (
      family_name %in%
      c(
        "poisson",
        "quasipoisson",
        "negative binomial",
        "negative binomial 1",
        "negative binomial 2"
      )
    ) {

      "observed mean counts"

    } else {

      "observed means"
    }

  } else if (
    identical(
      estimate_type,
      "predicted"
    )
  ) {

    if (
      identical(
        prediction_scale,
        "link"
      )
    ) {

      paste0(
        "model-based estimates on the ",
        link_name,
        " link scale"
      )

    } else if (
      family_name %in%
      c(
        "binomial",
        "quasibinomial"
      )
    ) {

      "model-based predicted probabilities"

    } else if (
      family_name %in%
      c(
        "poisson",
        "quasipoisson",
        "negative binomial",
        "negative binomial 1",
        "negative binomial 2"
      )
    ) {

      "model-based predicted mean counts"

    } else {

      "model-based estimates on the response scale"
    }

  } else {

    "displayed estimates"
  }

  # ==========================================================================
  # 4. INTRODUCTION
  # ==========================================================================

  interaction_display <- paste(
    interaction_variables,
    collapse = " \u00d7 "
  )

  introduction <- paste0(
    "This plot summarizes the contrasts associated with the ",
    interaction_display,
    " interaction in the fitted model for ",
    dependent_variable,
    "."
  )

  # ==========================================================================
  # 5. GRAPHICAL HIERARCHY
  # ==========================================================================

  hierarchy_display <- paste(
    hierarchy_variables,
    collapse = " \u2192 "
  )

  hierarchy_description <- paste0(
    "The graphical decomposition follows the hierarchy ",
    hierarchy_display,
    ", with ",
    split_variable,
    " defining the top-level split and ",
    target_variable,
    " defining the quantities represented on the x-axis."
  )

  # ==========================================================================
  # 6. CONTRAST CODING
  # ==========================================================================

  contrast_coding_description <- vapply(
    interaction_variables,
    function(current_variable) {

      current_row <- variable_information[
        variable_information$Variable ==
          current_variable &
          variable_information$Role ==
          "Predictor",
        ,
        drop = FALSE
      ]

      if (
        nrow(
          current_row
        ) == 0L
      ) {
        return(
          paste0(
            "The contrast coding for variable ",
            current_variable,
            " could not be identified."
          )
        )
      }

      current_contrast <- current_row$Contrast_Matrix_Type[1L]

      contrast_explanation <- switch(
        current_contrast,

        "Treatment" = {

          current_factor <- .extract_model_frame(
            model
          )[[current_variable]]

          current_matrix <- contrasts(
            current_factor
          )

          reference_index <- which(
            rowSums(
              abs(
                current_matrix
              )
            ) == 0
          )

          reference_level <- if (
            length(
              reference_index
            ) == 1L
          ) {
            rownames(
              current_matrix
            )[reference_index]
          } else {
            NA_character_
          }

          if (
            !is.na(
              reference_level
            )
          ) {
            paste0(
              "each non-reference level is compared with the reference category ",
              reference_level
            )
          } else {
            "each non-reference level is compared with the reference category"
          }
        },

        "Sum" =
          "each represented level is compared with the unweighted Grand Mean",

        "Helmert" =
          "each level is compared with the average of the preceding levels",

        "Reverse Helmert" =
          "each level is compared with the average of the subsequent levels",

        "Sliding Differences" =
          "each level is compared with the immediately preceding level",

        "Custom" =
          paste0(
            "each contrast represents a comparison specified by the ",
            "user-defined contrast matrix"
          ),

        paste0(
          "the comparisons follow the corresponding ",
          current_contrast,
          " contrast matrix"
        )
      )

      paste0(
        "Variable ",
        current_variable,
        " was coded using ",
        current_contrast,
        " contrasts, whereby ",
        contrast_explanation,
        "."
      )
    },
    character(1)
  )

  # ==========================================================================
  # 7. UPPER PANEL
  # ==========================================================================

  upper_description <- paste0(
    "Points represent the ",
    mean_description,
    " for the quantities involved in the displayed interaction contrasts, ",
    "and error bars represent their confidence intervals."
  )

  # ==========================================================================
  # 8. STYLE-SPECIFIC INTERACTION INTERPRETATION
  # ==========================================================================

  interaction_description <- switch(
    style,

    minimal = paste0(
      "Dashed segments connect the two quantities compared in each elementary contrast."
    ),

    tree = {
      if (prep$hierarchy_depth == 2L) {
        paste0(
          "Dashed segments connect the two quantities compared in each elementary contrast. ",
          "The resulting elementary contrasts are labelled \u03941, \u03942, and so forth. ",
          "The tree shows how these elementary contrasts are combined according to the contrast hierarchy ",
          "to obtain the final interaction contrast, which is marked by \u00d7."
        )
      } else {
        paste0(
          "Dashed segments connect the two quantities compared in each elementary contrast. ",
          "The resulting elementary contrasts are labelled \u03941, \u03942, and so forth. ",
          "The tree then shows how these elementary contrasts are progressively combined according to the contrast hierarchy. ",
          "Reading the tree from top to bottom, formulas at each successive level show how pairs of contrasts from the preceding level ",
          "are combined to form higher-order contrasts. ",
          "This process continues until the final interaction contrast is obtained, which is marked by \u00d7."
        )
      }
    }
  )


  # ==========================================================================
  # 8.1. ABBREVIATED FACTOR LABELS
  # ==========================================================================

  abbreviation_description <- NULL

  if (
    !is.null(label_key) &&
    is.data.frame(label_key) &&
    nrow(label_key) > 0L
  ) {

    abbreviation_mapping <- paste0(
      label_key$abbreviated,
      " = ",
      label_key$full
    )

    abbreviation_description <- paste0(
      "Abbreviated factor labels are used where needed for readability: ",
      paste(
        abbreviation_mapping,
        collapse = "; "
      ),
      "."
    )
  }

  # ==========================================================================
  # 9. COEFFICIENT LABELS
  # ==========================================================================

  label_description <- sprintf(
    paste0(
      "Labels Coef. 1, Coef. 2, and so forth identify the displayed interaction ",
      "coefficients; an asterisk indicates p < %.3f."
    ),
    sig_level
  )

  # ==========================================================================
  # 10. LOWER PANEL
  # ==========================================================================

  coefficient_description <- NULL

  if (
    isTRUE(
      show_coef
    )
  ) {

    lower_quantity <- if (
      identical(
        coefficient_scale,
        "model"
      )
    ) {
      "model coefficient"
    } else {
      "interpreted effect estimate"
    }

    coefficient_description <- paste0(
      "The lower panel displays the ",
      lower_quantity,
      " for each displayed interaction contrast, together with its confidence ",
      "interval and an approximate Cohen's d. Each coefficient is horizontally ",
      "aligned with the center of the corresponding hierarchical contrast block ",
      "in the upper panel."
    )
  }

  # Each interaction coefficient equals the product of the scale
  # factors of its components multiplied by the displayed difference
  # of differences on the linear-predictor scale.
  scale_description <- NULL

  scale_factors <- attr(contrast_output, "scale_factors")

  if (
    length(scale_factors) > 0L &&
    any(abs(scale_factors - 1) > 1e-8)
  ) {

    scale_list <- paste0(
      "s = ",
      .format_scale_factor(scale_factors),
      " for ",
      names(scale_factors),
      collapse = ", "
    )

    scale_description <- if (isTRUE(per_difference) && isTRUE(show_coef)) {
      paste0(
        "Under this coding, each interaction coefficient equals a scale factor ",
        "s multiplied by the displayed interaction contrast (",
        scale_list,
        "); the lower panel reports the coefficients divided by s, whereas ",
        "test statistics, p-values, and Cohen's d are unaffected by this ",
        "rescaling."
      )
    } else {
      paste0(
        "Under this coding, each interaction coefficient equals a scale factor ",
        "s multiplied by the displayed interaction contrast on the ",
        "linear-predictor scale (",
        scale_list,
        ")."
      )
    }
  }

  weighted_description <- if (
    any(grepl("W(", as.matrix(contrast_output), fixed = TRUE))
  ) {
    paste0(
      "Labels of the form W(a1 .20, a3 .80) denote weighted means of the ",
      "listed levels with the indicated weights, as implied by the ",
      "hypothesis matrix of the coding."
    )
  } else {
    NULL
  }

  # ==========================================================================
  # 11. CAPTION
  # ==========================================================================

  # Relation between the declared and the estimated comparisons. The
  # figure always shows the comparisons that the coefficients estimate.
  hypothesis_description <- NULL

  hypothesis_check <- attr(contrast_output, "hypothesis_check")

  if (!is.null(hypothesis_check) && nrow(hypothesis_check) > 0L) {

    failing <- hypothesis_check$Coefficient[hypothesis_check$Status %in% c("differs", "reversed")]
    partial <- hypothesis_check$Coefficient[hypothesis_check$Status == "partially declared"]

    hypothesis_description <- if (length(failing) > 0L) {
      paste0(
        "The coefficients ", paste(failing, collapse = ", "),
        " do not estimate the declared hypotheses of their components; ",
        "the figure shows the comparisons that the coefficients estimate."
      )
    } else if (length(partial) > 0L) {
      paste0(
        "The declared hypotheses are estimated by the components that were ",
        "declared; the other components were not checked."
      )
    } else {
      paste0(
        "Each displayed coefficient estimates the product of the declared ",
        "hypotheses of its components."
      )
    }
  }

  caption_parts <- c(
    introduction,
    hierarchy_description,
    contrast_coding_description,
    hypothesis_description,
    upper_description,
    interaction_description,
    weighted_description,
    abbreviation_description,
    label_description,
    coefficient_description,
    scale_description
  )

  caption_parts <- caption_parts[
    !is.na(
      caption_parts
    ) &
      nzchar(
        caption_parts
      )
  ]

  paste(
    caption_parts,
    collapse = " "
  )
}


#' Compute recommended export dimensions for a multi-predictor MoPlot
#'
#' @param prep Prepared plotting data returned by
#'   `.prepare_moplot_multpred_data()`.
#' @param style Graphical style used for the upper panel.
#' @param show_coef Logical. Whether the lower coefficient panel is present.
#'
#' @return A named numeric vector containing recommended width and height
#'   in inches.
#'
#' @keywords internal
#' @noRd
.moplot_multpred_export_size <- function(
    prep,
    style = c("minimal", "tree"),
    show_coef = TRUE
) {

  style <- match.arg(style)

  number_points <- nrow(prep$point_data)

  hierarchy_depth <- prep$hierarchy_depth

  if (
    is.null(hierarchy_depth) ||
    length(hierarchy_depth) != 1L ||
    is.na(hierarchy_depth)
  ) {
    hierarchy_depth <- 2L
  }

  number_coefficients <- nrow(prep$coefficient_data)

  # Horizontal space is driven primarily by the number of displayed
  # quantities. Delta and tree representations receive slightly more room.
  width_per_point <- switch(
    style,
    minimal = 0.42,
    tree = 0.50
  )

  width <- 6 +
    width_per_point *
    number_points

  # Multiple coefficient blocks need additional separation.
  if (number_coefficients > 1L) {

    width <- width +
      0.35 *
      (
        number_coefficients -
          1L
      )
  }

  width <- max(
    7,
    width
  )

  # Vertical space depends on the depth of the interaction hierarchy.
  height <- switch(
    style,
    minimal = 4.5 +
      0.35 *
      max(
        hierarchy_depth -
          2L,
        0L
      ),
    tree = 5.0 +
      1.25 *
      max(
        hierarchy_depth -
          2L,
        0L
      )
  )

  # Reserve space for the lower coefficient panel.
  if (isTRUE(show_coef)) {

    height <- height +
      2.5
  }

  c(
    width = round(
      width,
      1
    ),
    height = round(
      height,
      1
    )
  )
}
