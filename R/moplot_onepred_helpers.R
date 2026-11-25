#' Assign horizontal positions to displayed contrast quantities
#'
#' Internal helper used by `.prepare_moplot_onepred()`.
#'
#' Treatment and Sum contrasts are displayed as independent pairs.
#' For all other contrast types, repeated quantities reuse their previous
#' horizontal position.
#'
#' Structural matching is based on `Levels`, whereas axis labels are obtained
#' from `Display_Levels`.
#'
#' @param contrast_data A data frame returned by `contrasts_onepred()`.
#' @param contrast_type Character string identifying the contrast type.
#'
#' @return A list containing the updated contrast data and horizontal-axis
#'   information.
#'
#' @keywords internal
#' @noRd
.assign_moplot_positions <- function(
    contrast_data,
    contrast_type
) {

  required_columns <- c(
    "Contrast",
    "Levels",
    "Display_Levels"
  )

  missing_columns <- setdiff(required_columns,names(contrast_data))

  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        "`contrast_data` is missing the following required columns: ",
        paste(
          missing_columns,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  if (
    !is.character(contrast_type) ||
    length(contrast_type) != 1L ||
    is.na(contrast_type) ||
    !nzchar(contrast_type)
  ) {
    stop(
      "`contrast_type` must be a single non-empty character string.",
      call. = FALSE
    )
  }

  contrast_ids <- unique(contrast_data$Contrast)

  output <- contrast_data
  output$x <- NA_real_

  # Reusing the position of a quantity shared by two contrasts draws a
  # chain when the quantity lies on opposite sides of them (B - A and
  # C - B under sliding differences). When it lies on the same side, as
  # a common reference level, the segments would start from the same
  # point and overlap; each contrast then receives its own pair of
  # positions, as for treatment and sum coding. Within each contrast,
  # the first row is the positive side and the second the negative one.
  side <- stats::ave(
    seq_len(nrow(contrast_data)),
    contrast_data$Contrast,
    FUN = seq_along
  )

  shared_on_same_side <- anyDuplicated(
    paste(contrast_data$Levels, side, sep = "\r")
  ) > 0L

  paired_layout <- contrast_type %in% c("Treatment", "Sum") ||
    shared_on_same_side

  # ==========================================================================
  # PAIRED LAYOUT: TREATMENT AND SUM
  # ==========================================================================

  if (paired_layout) {

    for (contrast_index in seq_along(contrast_ids)) {

      contrast_id <- contrast_ids[contrast_index]

      row_indices <- which(output$Contrast == contrast_id)

      if (length(row_indices) != 2L) {
        stop(
          paste0(
            "Contrast `",
            contrast_id,
            "` must contain exactly two displayed quantities."
          ),
          call. = FALSE
        )
      }

      left_position <- 2 * (contrast_index - 1L) + 1L
      right_position <- left_position + 1L

      output$x[row_indices] <- c(right_position,left_position)
    }

    axis_information <- output[
      order(output$x),
      c(
        "x",
        "Display_Levels"
      ),
      drop = FALSE
    ]

    rownames(axis_information) <- NULL

    return(
      list(
        contrast_data = output,
        axis_breaks = axis_information$x,
        axis_labels = axis_information$Display_Levels,
        position_layout = "paired"
      )
    )
  }

  # ==========================================================================
  # REUSE LAYOUT: OTHER CONTRAST TYPES
  # ==========================================================================

  position_map <- numeric(0L)
  display_map <- character(0L)

  next_position <- 1L

  for (contrast_id in contrast_ids) {

    row_indices <- which(output$Contrast == contrast_id)

    structural_labels <- output$Levels[row_indices]

    display_labels <- output$Display_Levels[row_indices]

    previously_seen <- structural_labels %in% names(position_map)

    block_positions <- rep(NA_real_, length(row_indices))

    if (any(previously_seen)) {

      block_positions[previously_seen] <- unname(position_map[structural_labels[previously_seen]])

      new_indices <- which(!previously_seen)

      if (length(new_indices) > 0L) {

        block_positions[new_indices] <- seq(
          from = next_position,
          by = 1,
          length.out = length(new_indices)
        )
      }

    } else {

      block_positions <- rev(
        seq(
          from = next_position,
          by = 1,
          length.out = length(row_indices)
        )
      )
    }

    output$x[row_indices] <- block_positions

    for (index in seq_along(structural_labels)) {

      structural_label <- structural_labels[index]

      if (!structural_label %in% names(position_map)) {

        position_map[structural_label] <- block_positions[index]

        display_map[structural_label] <- display_labels[index]
      }
    }

    next_position <- max(unname(position_map)) + 1L
  }

  position_order <- order(unname(position_map))

  ordered_structural_labels <- names(position_map)[position_order]

  list(
    contrast_data = output,
    axis_breaks = unname(
      position_map[position_order]
    ),
    axis_labels = unname(
      display_map[ordered_structural_labels]
    ),
    position_layout = "reuse"
  )
}

#' Prepare data for a one-predictor contrast plot
#'
#' Internal helper used by `moplot_onepred()`.
#'
#' @param model A fitted model supported by the package.
#' @param predictor Character string identifying the focal categorical
#'   predictor.
#' @param sig_level Significance threshold.
#' @param conf_level Confidence level.
#' @param estimate_type Type of upper-panel estimates.
#' @param prediction_scale Prediction scale for model-based estimates.
#' @param coefficient_scale Scale used for the lower-panel estimate.
#' @param per_difference Logical. Should lower-panel coefficients be divided
#'   by their scale factor?
#' @param data Optional data frame with the original variables; recovered
#'   from the model when `NULL`.
#' @param bounds Optional natural range of the response, passed to
#'   `contrasts_onepred()`.
#' @param hypotheses Optional declared hypotheses, passed to
#'   `contrasts_onepred()`.
#'
#' @return An object of class `"moplot_onepred_data"`.
#'
#' @keywords internal
#' @noRd
.prepare_moplot_onepred <- function(
    model,
    predictor,
    sig_level = 0.05,
    conf_level = 0.95,
    estimate_type = c("auto", "observed", "predicted"),
    prediction_scale = c("response", "link"),
    coefficient_scale = c("model", "interpreted"),
    per_difference = FALSE,
    data = NULL,
    bounds = NULL,
    hypotheses = NULL
) {

  # ==========================================================================
  # 1. ARGUMENT MATCHING AND VALIDATION
  # ==========================================================================

  estimate_type <- match.arg(estimate_type)
  prediction_scale <- match.arg(prediction_scale)
  coefficient_scale <- match.arg(coefficient_scale)

  if (!.is_supported_model(model)) {
    stop(
      paste0(
        "Unsupported model class. Supported model classes include ",
        "lm, glm, lme4::lmer/glmer, lmerTest::lmer, nlme::lme, ",
        "and glmmTMB."
      ),
      call. = FALSE
    )
  }

  if (
    !is.character(predictor) ||
    length(predictor) != 1L ||
    is.na(predictor) ||
    !nzchar(predictor)
  ) {
    stop(
      "`predictor` must be a single non-empty character string.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(sig_level) ||
    length(sig_level) != 1L ||
    is.na(sig_level) ||
    !is.finite(sig_level) ||
    sig_level <= 0 ||
    sig_level >= 1
  ) {
    stop(
      "`sig_level` must be a single numeric value strictly between 0 and 1.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(conf_level) ||
    length(conf_level) != 1L ||
    is.na(conf_level) ||
    !is.finite(conf_level) ||
    conf_level <= 0 ||
    conf_level >= 1
  ) {
    stop(
      "`conf_level` must be a single numeric value strictly between 0 and 1.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # 2. MODEL INFORMATION
  # ==========================================================================

  model_frame <- .extract_model_frame(model)

  # Prediction data need the original variables (x for a term such as
  # log(x)), which the model frame does not contain.
  analysis_data <- if (is.null(data)) {
    .extract_analysis_data(model)
  } else {
    .align_data_with_model(model, data, model_frame)
  }

  if (!is.data.frame(analysis_data)) {
    stop(
      paste0(
        "The data used to fit the model could not be recovered. ",
        "Supply them through `data`."
      ),
      call. = FALSE
    )
  }

  if (!predictor %in% names(model_frame)) {
    stop(
      paste0(
        "The predictor `",
        predictor,
        "` was not found in the model frame."
      ),
      call. = FALSE
    )
  }

  if (!is.factor(model_frame[[predictor]])) {
    stop(
      paste0(
        "The focal predictor `",
        predictor,
        "` must be represented by a factor variable."
      ),
      call. = FALSE
    )
  }

  model_information <- model_contrast_table(model,type = "all")

  variable_summary <- model_information[["Variable-Level Summary"]]

  outcome_names <- variable_summary$Variable[variable_summary$Role == "Dependent"]

  if (length(outcome_names) != 1L) {
    stop(
      "`moplot_onepred()` requires exactly one dependent variable.",
      call. = FALSE
    )
  }

  outcome <- outcome_names[1L]

  predictor_row <- variable_summary[
    variable_summary$Variable == predictor &
      variable_summary$Role == "Predictor",
    ,
    drop = FALSE
  ]

  if (nrow(predictor_row) != 1L) {
    stop(
      paste0(
        "`",
        predictor,
        "` could not be uniquely identified as a model predictor."
      ),
      call. = FALSE
    )
  }

  contrast_type <- predictor_row$Contrast_Matrix_Type[1L]

  if (
    is.na(contrast_type) ||
    !nzchar(contrast_type)
  ) {
    stop(
      paste0(
        "The contrast matrix assigned to `",
        predictor,
        "` could not be identified."
      ),
      call. = FALSE
    )
  }

  is_polynomial <- identical(contrast_type, "Polynomial")

  # Polynomial components do not compare two quantities, so there is no
  # scale factor relating a coefficient to a displayed difference.
  if (is_polynomial && isTRUE(per_difference)) {
    warning(
      paste0(
        "`per_difference` is ignored for polynomial contrasts, whose ",
        "coefficients do not refer to a difference between two ",
        "displayed quantities."
      ),
      call. = FALSE
    )

    per_difference <- FALSE
  }

  # ==========================================================================
  # 3. CHECK FOCAL INTERACTIONS
  # ==========================================================================

  fixed_formula <- .extract_fixed_formula(model)

  terms_object <- stats::terms(fixed_formula)

  term_labels <- attr(terms_object,"term.labels")

  interaction_terms <- term_labels[
    grepl(
      ":",
      term_labels,
      fixed = TRUE
    )
  ]

  if (length(interaction_terms) > 0L) {

    focal_interaction <- vapply(
      interaction_terms,
      function(term_label) {

        interaction_variables <- all.vars(
          stats::as.formula(
            paste("~", term_label)
          )
        )

        predictor %in% interaction_variables
      },
      logical(1)
    )

    if (any(focal_interaction)) {
      stop(
        paste0(
          "`moplot_onepred()` does not currently support models in which ",
          "the focal predictor is involved in an interaction."
        ),
        call. = FALSE
      )
    }
  }

  # ==========================================================================
  # 4. RESOLVE AUTOMATIC ESTIMATE TYPE
  # ==========================================================================

  # The same rule as contrasts_onepred(); see .resolve_estimate_type().
  resolved_estimate_type <- .resolve_estimate_type(
    model = model,
    predictor = predictor,
    estimate_type = estimate_type,
    prediction_scale = prediction_scale
  )

  # ==========================================================================
  # 5. CONTRAST ESTIMATES
  # ==========================================================================

  contrast_data <- contrasts_onepred(
    model = model,
    predictor = predictor,
    outcome = outcome,
    data = analysis_data,
    sig_level = sig_level,
    conf_level = conf_level,
    estimate_type = resolved_estimate_type,
    prediction_scale = prediction_scale,
    bounds = bounds,
    hypotheses = hypotheses
  )

  contrast_family <- attr(contrast_data,"family")

  contrast_link <- attr(contrast_data,"link")

  contrast_estimate_type <- attr(contrast_data,"estimate_type")

  contrast_prediction_scale <- attr(contrast_data,"prediction_scale")



  if (
    !is.data.frame(contrast_data) ||
    nrow(contrast_data) == 0L
  ) {
    stop(
      "`contrasts_onepred()` did not return any interpretable contrasts.",
      call. = FALSE
    )
  }

  required_contrast_columns <- c(
    "Contrast",
    "Levels",
    "Estimate",
    "SE",
    "CI_low",
    "CI_high",
    "Significance"
  )

  missing_contrast_columns <- setdiff(required_contrast_columns,names(contrast_data))

  if (length(missing_contrast_columns) > 0L) {
    stop(
      paste0(
        "`contrasts_onepred()` did not return the expected columns: ",
        paste(missing_contrast_columns, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  contrast_ids <- unique(contrast_data$Contrast)

  rows_per_contrast <- table(contrast_data$Contrast)

  if (!is_polynomial && any(rows_per_contrast != 2L)) {
    stop(
      paste0(
        "Each non-polynomial contrast must contain exactly two displayed ",
        "quantities."
      ),
      call. = FALSE
    )
  }

  contrast_disclaimer <- attr(contrast_data,"disclaimer")

  contrast_hypothesis_check <- attr(contrast_data, "hypothesis_check")

  contrast_coding_notes <- attr(contrast_data, "coding_notes")

  # ==========================================================================
  # 6. HORIZONTAL POSITIONS
  # ==========================================================================

  if (is_polynomial) {

    # Every component involves all levels, which keep a single position
    # spaced according to the scores of the coding.
    level_rows <- contrast_data[contrast_data$Contrast == contrast_ids[1L], , drop = FALSE]

    level_positions <- .polynomial_level_positions(level_rows$Weight)

    contrast_data$x <- level_positions[match(contrast_data$Levels, level_rows$Levels)]

    axis_breaks <- level_positions

    axis_labels <- level_rows$Display_Levels

    position_layout <- "polynomial"

  } else {

    position_information <- .assign_moplot_positions(
      contrast_data = contrast_data,
      contrast_type = contrast_type
    )

    contrast_data <- position_information$contrast_data

    axis_breaks <- position_information$axis_breaks

    axis_labels <- position_information$axis_labels

    position_layout <- position_information$position_layout
  }

  # ==========================================================================
  # 7. SEGMENT TABLE
  # ==========================================================================

  # Polynomial components are drawn as curves (section 9) instead.
  segment_list <- if (is_polynomial) list() else lapply(
    seq_along(contrast_ids),
    function(contrast_index) {

      contrast_id <- contrast_ids[contrast_index]

      current_contrast <- contrast_data[
        contrast_data$Contrast == contrast_id,
        ,
        drop = FALSE
      ]

      current_contrast <- current_contrast[
        order(current_contrast$x),
        ,
        drop = FALSE
      ]

      significance_value <- current_contrast$Significance[1L]

      significant <- !is.na(significance_value) &&
        identical(significance_value, "sig")

      data.frame(
        Contrast = contrast_id,
        Contrast_number = contrast_index,
        x = current_contrast$x[1L],
        xend = current_contrast$x[2L],
        y = current_contrast$Estimate[1L],
        yend = current_contrast$Estimate[2L],
        x_mid = mean(current_contrast$x),
        y_mid = mean(current_contrast$Estimate),
        Significant = significant,
        Label = paste0(
          "C",
          contrast_index,
          if (significant) "*" else ""
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  segment_data <- if (is_polynomial) NULL else do.call(rbind,segment_list)

  if (!is_polynomial) {
    rownames(segment_data) <- NULL
  }

  # ==========================================================================
  # 8. MODEL PARAMETERS
  # ==========================================================================

  parameter_table <- extract_param(
    model = model,
    conf_level = conf_level,
    include_intercept = FALSE
  )

  focal_matrix_columns <- .focal_coefficient_names(
    model = model,
    predictor = predictor
  )

  parameter_indices <- match(focal_matrix_columns,parameter_table$Parameter)

  if (anyNA(parameter_indices)) {

    unmatched_parameters <- focal_matrix_columns[
      is.na(parameter_indices)
    ]

    stop(
      paste0(
        "The following model-matrix columns could not be matched to ",
        "`extract_param()` output: ",
        paste(unmatched_parameters, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  focal_parameters <- parameter_table[
    parameter_indices,
    ,
    drop = FALSE
  ]

  if (nrow(focal_parameters) != length(contrast_ids)) {
    stop(
      paste0(
        "The number of model parameters associated with `",
        predictor,
        "` does not match the number of displayed contrasts."
      ),
      call. = FALSE
    )
  }

  parameter_data <- data.frame(
    Contrast = contrast_ids,
    Contrast_number = seq_along(contrast_ids),
    Parameter = focal_parameters$Parameter,
    # Polynomial coefficients sit in their own panel, one per component.
    x = if (is_polynomial) seq_along(contrast_ids) else segment_data$x_mid,
    Estimate = focal_parameters$Estimate,
    CI_low = focal_parameters$CI_low,
    CI_high = focal_parameters$CI_high,
    Interpreted_Estimate =
      focal_parameters$Interpreted_Estimate,
    Interpreted_CI_low =
      focal_parameters$Interpreted_CI_low,
    Interpreted_CI_high =
      focal_parameters$Interpreted_CI_high,
    Estimate_type =
      focal_parameters$Estimate_type,
    d = focal_parameters$d,
    d_CI_low = focal_parameters$d_CI_low,
    d_CI_high = focal_parameters$d_CI_high,
    stringsAsFactors = FALSE
  )

  ratio_types <- c(
    "Odds ratio",
    "Risk ratio",
    "Incidence rate ratio",
    "Hazard ratio",
    "Multiplicative ratio",
    "Exponentiated coefficient"
  )

  # Each coefficient equals its scale factor multiplied by the displayed
  # difference on the linear-predictor scale (for example 1/3 for the
  # second Helmert contrast).
  # Polynomial coefficients refer to the components as coded.
  parameter_data$Scale_factor <- if (is_polynomial) {
    rep(1, length(contrast_ids))
  } else {
    contrast_data$Scale_factor[
      match(contrast_ids, contrast_data$Contrast)
    ]
  }

  # Cohen's d converted from the log odds ratio refers to the
  # coefficient as coded; dividing by the scale factor makes it refer to
  # the displayed difference. Conversions from t or z are already
  # independent of the scaling of the contrast.
  d_from_log_odds <- focal_parameters$d_method %in% "log odds ratio"

  parameter_data$d[d_from_log_odds] <- parameter_data$d[d_from_log_odds] /
    parameter_data$Scale_factor[d_from_log_odds]
  parameter_data$d_CI_low[d_from_log_odds] <- parameter_data$d_CI_low[d_from_log_odds] /
    parameter_data$Scale_factor[d_from_log_odds]
  parameter_data$d_CI_high[d_from_log_odds] <- parameter_data$d_CI_high[d_from_log_odds] /
    parameter_data$Scale_factor[d_from_log_odds]

  if (isTRUE(per_difference)) {

    # Dividing by the scale factor expresses each coefficient per unit
    # of the displayed difference. Test statistics, p-values, and
    # Cohen's d are unchanged by this rescaling.
    scale_factors <- parameter_data$Scale_factor
    is_ratio <- parameter_data$Estimate_type %in% ratio_types

    parameter_data$Estimate <- parameter_data$Estimate / scale_factors
    parameter_data$CI_low <- parameter_data$CI_low / scale_factors
    parameter_data$CI_high <- parameter_data$CI_high / scale_factors

    parameter_data$Interpreted_Estimate <- ifelse(
      is_ratio,
      exp(parameter_data$Estimate),
      parameter_data$Interpreted_Estimate / scale_factors
    )

    parameter_data$Interpreted_CI_low <- ifelse(
      is_ratio,
      exp(parameter_data$CI_low),
      parameter_data$Interpreted_CI_low / scale_factors
    )

    parameter_data$Interpreted_CI_high <- ifelse(
      is_ratio,
      exp(parameter_data$CI_high),
      parameter_data$Interpreted_CI_high / scale_factors
    )
  }

  if (identical(coefficient_scale, "model")) {

    parameter_data$Display_estimate <- parameter_data$Estimate

    parameter_data$Display_CI_low <- parameter_data$CI_low

    parameter_data$Display_CI_high <- parameter_data$CI_high

    parameter_data$Display_type <- "Model estimate"

    null_value <- 0

  } else {

    parameter_data$Display_estimate <- parameter_data$Interpreted_Estimate

    parameter_data$Display_CI_low <- parameter_data$Interpreted_CI_low

    parameter_data$Display_CI_high <- parameter_data$Interpreted_CI_high

    parameter_data$Display_type <- parameter_data$Estimate_type

    null_value <- if (
      all(parameter_data$Estimate_type %in% ratio_types)
    ) {
      1
    } else {
      0
    }
  }

  # ==========================================================================
  # 9. RAW DATA
  # ==========================================================================

  raw_level_values <- as.character(analysis_data[[predictor]])

  # Raw observations are shown on the scale of the model response
  # (log(y) for a log-transformed response, 0/1 for binary factors).
  raw_outcome_values <- .response_values(
    model = model,
    data = analysis_data,
    outcome = outcome
  )

  raw_data_list <- lapply(
    levels(analysis_data[[predictor]]),
    function(current_level) {

      # A level can intentionally appear at multiple horizontal positions,
      # as in the paired layout used for Treatment contrasts.
      current_positions <- unique(contrast_data$x[contrast_data$Levels == current_level])

      current_positions <- current_positions[is.finite(current_positions)]

      # Aggregated quantities such as "Grand Mean" or "A, B" have no
      # corresponding raw observations. Likewise, a factor level that is
      # represented only inside an aggregate is not plotted separately.
      if (length(current_positions) == 0L) {
        return(NULL)
      }

      current_rows <- which(raw_level_values == current_level)

      if (length(current_rows) == 0L) {
        return(NULL)
      }

      do.call(
        rbind,
        lapply(
          current_positions,
          function(current_position) {

            data.frame(
              x = rep(
                current_position,
                length(current_rows)
              ),
              y = raw_outcome_values[current_rows],
              Level = rep(
                current_level,
                length(current_rows)
              ),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )

  raw_data_list <- Filter(
    Negate(is.null),
    raw_data_list
  )

  if (length(raw_data_list) == 0L) {

    raw_data <- data.frame(
      x = numeric(),
      y = numeric(),
      Level = character(),
      stringsAsFactors = FALSE
    )

  } else {

    raw_data <- do.call(rbind,raw_data_list)

    rownames(raw_data) <- NULL
  }

  # ==========================================================================
  # 10. POLYNOMIAL COMPONENTS
  # ==========================================================================

  polynomial_data <- NULL

  if (is_polynomial) {

    family_object <- .extract_family(model)

    # Same rule as contrasts_onepred(): estimates of generalized models are
    # on the response scale unless link-scale predictions were requested.
    estimates_on_response_scale <- .is_generalized_model(model) &&
      !is.null(family_object) &&
      (
        identical(contrast_estimate_type, "observed") ||
          identical(contrast_prediction_scale, "response")
      )

    polynomial_data <- .polynomial_component_data(
      contrast_data = contrast_data,
      coefficients = focal_parameters$Estimate,
      significant = contrast_data$Significance[
        match(contrast_ids, contrast_data$Contrast)
      ] %in% "sig",
      linkfun = if (estimates_on_response_scale) family_object$linkfun else identity,
      linkinv = if (estimates_on_response_scale) family_object$linkinv else identity
    )
  }

  # ==========================================================================
  # 11. RETURN
  # ==========================================================================

  output <- list(
    model = model,
    outcome = outcome,
    predictor = predictor,
    contrast_type = contrast_type,

    family = contrast_family,
    link = contrast_link,
    estimate_type = contrast_estimate_type,
    prediction_scale = contrast_prediction_scale,

    coefficient_scale = coefficient_scale,
    per_difference = isTRUE(per_difference),
    position_layout = position_layout,
    axis_breaks = axis_breaks,
    axis_labels = axis_labels,
    contrast_data = contrast_data,
    segment_data = segment_data,
    parameter_data = parameter_data,
    raw_data = raw_data,
    polynomial_data = polynomial_data,
    null_value = null_value,
    disclaimer = contrast_disclaimer,
    hypothesis_check = contrast_hypothesis_check,
    coding_notes = contrast_coding_notes
  )

  class(output) <- c("moplot_onepred_data","list")

  output
}


#' Build the upper panel of a one-predictor contrast plot
#'
#' Internal helper used by `moplot_onepred()`.
#'
#' @param prep An object returned by `.prepare_moplot_onepred()`.
#' @param show_data Logical. Should raw observations be displayed using
#'   horizontal jittering?
#' @param base_size Base font size passed to `ggplot2::theme_light()`.
#'
#' @return A `ggplot` object.
#'
#' @keywords intrnal build for the upper panel
#'
#' @param prep An object returned by `.prepare_moplot_onepred()`.
#'
#' @param show_data Logical. Should raw observations be displayed using
#'   horizontal jittering? The default is `FALSE`.
#'
#' @param style Character string specifying the graphical style. Available
#'   options are:
#'
#'   \describe{
#'     \item{`"minimal"`}{
#'       Uses a restrained grayscale palette and relatively compact graphical
#'       elements.
#'     }
#'     \item{`"publication"`}{
#'       Uses a colorblind-safe palette, stronger graphical emphasis, and
#'       larger contrast-midpoint markers.
#'     }
#'   }
#'
#' @param base_size Base font size passed to `ggplot2::theme_light()`.
#'   The default is `13`.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd
.plot_moplot_onepred_upper <- function(
    prep,
    show_data = FALSE,
    style = c("minimal", "publication"),
    base_size = 13,
    label_key=NULL
) {

  # ==========================================================================
  # 1. ARGUMENT MATCHING AND VALIDATION
  # ==========================================================================

  style <- match.arg(style)

  if (
    !inherits(prep, "moplot_onepred_data") ||
    !is.list(prep)
  ) {
    stop(
      paste0(
        "`prep` must be an object returned by ",
        "`.prepare_moplot_onepred()`."
      ),
      call. = FALSE
    )
  }

  if (
    !is.logical(show_data) ||
    length(show_data) != 1L ||
    is.na(show_data)
  ) {
    stop(
      "`show_data` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

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

  contrast_data <- prep$contrast_data
  segment_data <- prep$segment_data
  raw_data <- prep$raw_data

  axis_breaks <- as.numeric(prep$axis_breaks)

  axis_labels <- as.character(prep$axis_labels)

  axis_labels <- .apply_moplot_label_key(labels = axis_labels,label_key = label_key)

  required_contrast_columns <- c(
    "Estimate",
    "CI_low",
    "CI_high",
    "x"
  )

  missing_contrast_columns <- setdiff(required_contrast_columns,names(contrast_data))

  if (length(missing_contrast_columns) > 0L) {
    stop(
      paste0(
        "`prep$contrast_data` is missing the following columns: ",
        paste(
          missing_contrast_columns,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  required_segment_columns <- c(
    "x",
    "xend",
    "y",
    "yend",
    "x_mid",
    "y_mid",
    "Significant",
    "Label"
  )

  missing_segment_columns <- setdiff(required_segment_columns,names(segment_data))

  if (length(missing_segment_columns) > 0L) {
    stop(
      paste0(
        "`prep$segment_data` is missing the following columns: ",
        paste(
          missing_segment_columns,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  if (
    length(axis_breaks) == 0L ||
    length(axis_breaks) != length(axis_labels) ||
    anyNA(axis_breaks) ||
    any(!is.finite(axis_breaks))
  ) {
    stop(
      "The horizontal-axis information in `prep` is invalid.",
      call. = FALSE
    )
  }

  contrast_data$x <- as.numeric(contrast_data$x)

  segment_data$x <- as.numeric(segment_data$x)

  segment_data$xend <- as.numeric(segment_data$xend)

  segment_data$x_mid <- as.numeric(segment_data$x_mid)

  segment_data$y_mid <- as.numeric(segment_data$y_mid)

  raw_data$x <- as.numeric(raw_data$x)

  raw_data$y <- as.numeric(raw_data$y)

  # ==========================================================================
  # 2. STYLE SETTINGS
  # ==========================================================================

  grayscale_palette <- c(
    raw_data = "#D9D9D9",
    estimate = "#737373",
    segment = "#525252",
    midpoint = "#252525",
    guide = "#BDBDBD",
    text = "#252525"
  )

  publication_palette <- c(
    raw_data = "#D9D9D9",
    estimate = "#737373",
    guide = "#BDBDBD",
    text = "#252525"
  )

  if (identical(style, "minimal")) {

    palette <- grayscale_palette

    estimate_point_size <- 2.8
    estimate_linewidth <- 0.60
    segment_linewidth <- 0.75

    midpoint_size <- 4.5
    midpoint_stroke <- 1.0

    contrast_label_size <- 4.0

    raw_point_size <- 1.5
    raw_alpha <- 0.25

    guide_linewidth <- 0.40

    segment_data$Segment_colour <- palette[["segment"]]
    segment_data$Midpoint_colour <- palette[["midpoint"]]

  } else {

    palette <- publication_palette

    estimate_point_size <- 3.5
    estimate_linewidth <- 0.80
    segment_linewidth <- 1.00

    midpoint_size <- 5.5
    midpoint_stroke <- 1.2

    contrast_label_size <- 4.6

    raw_point_size <- 1.8
    raw_alpha <- 0.30

    guide_linewidth <- 0.45

    contrast_ids <- sort(unique(segment_data$Contrast_number))

    contrast_colours <- stats::setNames(
      .moplot_contrast_palette(
        length(contrast_ids)
      ),
      as.character(contrast_ids)
    )

    segment_data$Segment_colour <- unname(
      contrast_colours[
        as.character(
          segment_data$Contrast_number
        )
      ]
    )

    segment_data$Midpoint_colour <- segment_data$Segment_colour
  }

  # ==========================================================================
  # 3. Y-AXIS LIMITS
  # ==========================================================================

  displayed_y <- c(
    contrast_data$Estimate,
    contrast_data$CI_low,
    contrast_data$CI_high,
    segment_data$y,
    segment_data$yend,
    segment_data$y_mid
  )

  if (isTRUE(show_data)) {
    displayed_y <- c(
      displayed_y,
      raw_data$y
    )
  }

  displayed_y <- displayed_y[is.finite(displayed_y)]

  if (length(displayed_y) == 0L) {
    stop(
      "No finite values were available for the upper-panel y-axis.",
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

  y_lower <- y_range[1L] - 0.10 * y_span

  y_upper <- y_range[2L] + 0.24 * y_span

  contrast_label_y <- y_range[2L] + 0.15 * y_span

  # ==========================================================================
  # 4. X-AXIS LIMITS
  # ==========================================================================

  displayed_x <- c(
    axis_breaks,
    contrast_data$x,
    segment_data$x,
    segment_data$xend,
    segment_data$x_mid
  )

  if (isTRUE(show_data)) {
    displayed_x <- c(
      displayed_x,
      raw_data$x
    )
  }

  displayed_x <- displayed_x[is.finite(displayed_x)]

  if (length(displayed_x) == 0L) {
    stop(
      "No finite values were available for the upper-panel x-axis.",
      call. = FALSE
    )
  }

  x_range <- range(displayed_x)

  x_limits <- c(x_range[1L] - 0.5,x_range[2L] + 0.5)

  # ==========================================================================
  # 5. BASE PLOT
  # ==========================================================================

  plot <- ggplot2::ggplot()

  # ==========================================================================
  # 6. OPTIONAL RAW DATA
  # ==========================================================================

  if (isTRUE(show_data)) {

    valid_raw_data <- raw_data[
      is.finite(raw_data$x) &
        is.finite(raw_data$y),
      ,
      drop = FALSE
    ]

    plot <- plot +
      ggplot2::geom_jitter(
        data = valid_raw_data,
        mapping = ggplot2::aes(
          x = x,
          y = y
        ),
        width = 0.06,
        height = 0,
        alpha = raw_alpha,
        size = raw_point_size,
        colour = palette[["raw_data"]],
        inherit.aes = FALSE
      )
  }

  # ==========================================================================
  # 7. VERTICAL MIDPOINT GUIDES
  # ==========================================================================

  plot <- plot +
    ggplot2::geom_vline(
      xintercept = segment_data$x_mid,
      linewidth = guide_linewidth,
      linetype = "dotted",
      colour = palette[["guide"]]
    )

  # ==========================================================================
  # 8. CONTRAST SEGMENTS
  # ==========================================================================

  plot <- plot +
    ggplot2::geom_segment(
      data = segment_data,
      mapping = ggplot2::aes(
        x = x,
        xend = xend,
        y = y,
        yend = yend,
        colour = Segment_colour
      ),
      linewidth = segment_linewidth,
      linetype = "dashed",
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_identity()

  # ==========================================================================
  # 9. ESTIMATES AND CONFIDENCE INTERVALS
  # ==========================================================================

  plot <- plot +
    ggplot2::geom_errorbar(
      data = contrast_data,
      mapping = ggplot2::aes(
        x = x,
        ymin = CI_low,
        ymax = CI_high
      ),
      width = 0.06,
      linewidth = estimate_linewidth,
      colour = palette[["estimate"]],
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = contrast_data,
      mapping = ggplot2::aes(
        x = x,
        y = Estimate
      ),
      size = estimate_point_size,
      colour = palette[["estimate"]],
      inherit.aes = FALSE
    )

  # ==========================================================================
  # 10. CONTRAST MIDPOINT MARKERS
  # ==========================================================================

  plot <- plot +
    ggplot2::geom_point(
      data = segment_data,
      mapping = ggplot2::aes(
        x = x_mid,
        y = y_mid,
        colour = Midpoint_colour
      ),
      shape = 4,
      size = midpoint_size,
      stroke = midpoint_stroke,
      show.legend = FALSE,
      inherit.aes = FALSE
    )

  # ==========================================================================
  # 11. CONTRAST LABELS
  # ==========================================================================

  plot <- plot +
    ggplot2::geom_text(
      data = segment_data,
      mapping = ggplot2::aes(
        x = x_mid,
        y = contrast_label_y,
        label = Label,
        colour = Midpoint_colour
      ),
      size = contrast_label_size,
      fontface = if (
        identical(
          style,
          "publication"
        )
      ) {
        "bold"
      } else {
        "plain"
      },
      vjust = 0,
      show.legend = FALSE,
      inherit.aes = FALSE
    )

  # ==========================================================================
  # 12. SCALES AND THEME
  # ==========================================================================

  plot <- plot +
    ggplot2::scale_x_continuous(
      breaks = axis_breaks,
      labels = axis_labels
    ) +
    ggplot2::coord_cartesian(
      xlim = x_limits,
      ylim = c(
        y_lower,
        y_upper
      ),
      clip = "off"
    ) +
    ggplot2::labs(
      x = prep$predictor,
      y = prep$outcome
    ) +
    ggplot2::theme_light(
      base_size = base_size
    ) +
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
            identical(
              style,
              "publication"
            )
          ) {
            "bold"
          } else {
            "plain"
          }
        ),

      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          colour = palette[["text"]]
        ),

      axis.text.y =
        ggplot2::element_text(
          colour = palette[["text"]]
        ),

      plot.margin =
        ggplot2::margin(
          t = 14,
          r = 12,
          b = 10,
          l = 10
        )
    )

  plot
}


#' Generate the caption for a one-predictor contrast plot
#'
#' Internal helper used by `moplot_onepred()`.
#'
#' @param model A fitted model.
#' @param predictor Character string identifying the focal predictor.
#' @param contrast_output Output returned by `contrasts_onepred()`.
#' @param sig_level Significance threshold.
#' @param show_data Whether raw observations are displayed.
#' @param show_coef Whether the coefficient panel is displayed.
#' @param coefficient_scale Scale used for the lower-panel estimates.
#' @param per_difference Whether lower-panel coefficients are divided by
#'   their scale factor.
#' @param hypothesis_check Optional table of declared hypotheses checked by
#'   `contrasts_onepred()`.
#' @param coding_notes Notes on coding columns whose signs do not describe
#'   the estimated comparisons.
#'
#' @return A character vector of length one.
#'
#' @keywords internal
#' @noRd
.moplot_onepred_caption <- function(
    model,
    predictor,
    contrast_output,
    sig_level = 0.05,
    show_data = FALSE,
    show_coef = TRUE,
    coefficient_scale = c("model", "interpreted"),
    per_difference = FALSE,
    family,
    link,
    estimate_type,
    prediction_scale,
    estimate_label,
    label_key = NULL,
    hypothesis_check = NULL,
    coding_notes = character()
) {

  coefficient_scale <- match.arg(coefficient_scale)

  family_name <- tolower(as.character(family)[1L])

  link_name <- tolower(as.character(link)[1L])

  estimate_type <- tolower(as.character(estimate_type)[1L])

  prediction_scale <- tolower(as.character(prediction_scale)[1L])

  model_information <- model_contrast_table(model,type = "all")

  variable_information <- model_information[["Variable-Level Summary"]]

  dependent_variables <- variable_information$Variable[variable_information$Role == "Dependent"]

  predictor_variables <- variable_information$Variable[variable_information$Role == "Predictor"]

  if (length(dependent_variables) != 1L) {
    stop(
      "`.moplot_onepred_caption()` requires exactly one dependent variable.",
      call. = FALSE
    )
  }

  if (!predictor %in% predictor_variables) {
    stop(
      paste0(
        "`",
        predictor,
        "` was not identified as a model predictor."
      ),
      call. = FALSE
    )
  }

  predictor_row <- variable_information[
    variable_information$Variable == predictor &
      variable_information$Role == "Predictor",
    ,
    drop = FALSE
  ]

  contrast_type <-predictor_row$Contrast_Matrix_Type[1L]

  dependent_variable <-dependent_variables[1L]

  mean_description <- if (identical(estimate_type, "observed")) {

    if (family_name %in% c("binomial", "quasibinomial")) {

      mean_description <- "observed proportions"

    } else if (
      family_name %in% c(
        "poisson",
        "quasipoisson",
        "negative binomial",
        "negative binomial 1",
        "negative binomial 2"
      )
    ) {

      mean_description <- "observed mean counts"

    } else {

      mean_description <- "observed means"
    }

  } else if (identical(estimate_type, "predicted")) {

    if (identical(prediction_scale, "link")) {

      mean_description <- paste0(
        "model-based estimates on the ",
        link_name,
        " link scale"
      )

    } else if (
      family_name %in% c("binomial", "quasibinomial")
    ) {

      mean_description <- "model-based predicted probabilities"

    } else if (
      family_name %in% c(
        "poisson",
        "quasipoisson",
        "negative binomial",
        "negative binomial 1",
        "negative binomial 2"
      )
    ) {

      mean_description <- "model-based predicted mean counts"

    } else {

      mean_description <- "model-based estimates on the response scale"
    }

  } else {

    mean_description <- "displayed estimates"
  }

  is_polynomial <- identical(contrast_type, "Polynomial")

  # The weights of the linear component are proportional to the centered
  # scores of the levels.
  equally_spaced_scores <- if (is_polynomial) {
    linear_weights <- contrast_output$Weight[
      contrast_output$Contrast == contrast_output$Contrast[1L]
    ]
    isTRUE(all.equal(diff(linear_weights), rep(diff(linear_weights)[1L], length(linear_weights) - 1L)))
  } else {
    NA
  }

  contrast_description <- switch(
    contrast_type,

    "Polynomial" = paste0(
      "Under polynomial coding, each coefficient estimates one trend ",
      "component (linear, quadratic, and so forth) across the ordered levels, ",
      if (isTRUE(equally_spaced_scores)) {
        "which are assumed to be equally spaced."
      } else {
        "spaced according to the scores of the coding, as on the horizontal axis."
      }
    ),

    "Treatment" = paste0(
      "Under treatment coding, each non-reference level is compared ",
      "with the reference category."
    ),

    "Sum" = paste0(
      "Under sum coding, each represented level is compared with the ",
      "unweighted Grand Mean of the displayed level-specific estimates."
    ),

    "Helmert" = paste0(
      "Under Helmert coding, each level is compared with the average ",
      "of the preceding levels."
    ),

    "Reverse Helmert" = paste0(
      "Under Reverse Helmert coding, each level is compared with the ",
      "average of the subsequent levels."
    ),

    "Sliding Differences" = paste0(
      "Under sliding-differences coding, each level is compared with ",
      "the immediately preceding level."
    ),

    "Custom" = paste0(
      "Under custom coding, each contrast represents a user-defined ",
      "comparison rather than a standard contrast pattern."
    ),

    paste0(
      "The focal predictor uses ",
      contrast_type,
      " contrast coding."
    )
  )

  segment_description <- switch(
    contrast_type,

    "Treatment" = paste0(
      "Dashed segments connect each non-reference level with the ",
      "reference category."
    ),

    "Sum" = paste0(
      "Dashed segments connect each represented level with the ",
      "unweighted Grand Mean."
    ),

    "Helmert" = paste0(
      "Dashed segments connect each focal level with the average of ",
      "the preceding levels."
    ),

    "Reverse Helmert" = paste0(
      "Dashed segments connect each focal level with the average of ",
      "the subsequent levels."
    ),

    "Sliding Differences" = paste0(
      "Dashed segments connect each level with the immediately ",
      "preceding level."
    ),

    "Custom" = paste0(
      "Dashed segments connect the quantities specified by each ",
      "user-defined contrast."
    ),

    "Dashed segments connect the two quantities compared within each contrast."
  )

  upper_description <- if (is_polynomial) {
    paste0(
      "Each upper panel shows one component. Points represent the ",
      mean_description,
      " for each level, and error bars represent their confidence ",
      "intervals. The solid curve shows the pattern implied by that ",
      "component alone, drawn around the mean of the estimates on the ",
      "linear-predictor scale (dotted horizontal line); crosses mark its ",
      "value at each level. Adding all components to that mean reproduces ",
      "the estimates on the linear-predictor scale."
    )
  } else {
    paste0(
      "Points in the upper panel represent the ",
      mean_description,
      " for the quantities involved in each comparison, and error bars ",
      "represent their confidence intervals. ",
      segment_description
    )
  }

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
  label_description <- sprintf(
    paste0(
      "Labels C1 and C2 and so forth identify the corresponding ",
      if (is_polynomial) "components" else "contrasts",
      "; an asterisk indicates p < %.3f."
    ),
    sig_level
  )

  raw_description <- if (isTRUE(show_data)) {
    paste0(
      "Raw observations are displayed in the background using ",
      "horizontal jittering."
    )
  } else {
    NULL
  }

  coefficient_description <- NULL

  if (isTRUE(show_coef)) {

    if (identical(coefficient_scale, "model")) {

      lower_quantity <- "model coefficient"

    } else {

      lower_quantity <- switch(
        estimate_label,

        "Coefficient" = "model coefficient",
        "Odds ratio" = "odds ratio",
        "Risk ratio" = "risk ratio",
        "Incidence rate ratio" = "incidence rate ratio",
        "Hazard ratio" = "hazard ratio",
        "Multiplicative ratio" = "multiplicative ratio",
        "Exponentiated coefficient" = "exponentiated coefficient",
        "Transformed coefficient" = "transformed coefficient",

        "transformed effect estimate"
      )
    }

    coefficient_description <- if (is_polynomial) {
      paste0(
        "The lower panel displays the ",
        lower_quantity,
        " for each component, together with its confidence interval and an ",
        "approximate Cohen's d, beneath the panel of that component."
      )
    } else {
      paste0(
        "The lower panel displays the ",
        lower_quantity,
        " for each contrast, together with its confidence interval and an ",
        "approximate Cohen's d. Each coefficient (indicated as Coef. 1, Coef. 2, and so forth) is horizontally aligned with ",
        "the center of the corresponding comparison in the upper panel."
      )
    }
  }

  # Groups with unequal weights are labelled W(level weight, ...).
  weighted_description <- if (
    any(grepl("^W\\(", contrast_output$Levels))
  ) {
    paste0(
      "Labels of the form W(A .20, C .80) denote weighted means of the ",
      "listed levels with the indicated weights, as implied by the ",
      "hypothesis matrix of the coding."
    )
  } else {
    NULL
  }

  # Each coefficient equals its scale factor s multiplied by the
  # displayed difference on the linear-predictor scale.
  scale_description <- NULL

  scale_factors <- contrast_output$Scale_factor[
    !duplicated(contrast_output$Contrast)
  ]

  if (
    length(scale_factors) > 0L &&
    any(abs(scale_factors - 1) > 1e-8)
  ) {

    scale_list <- paste0(
      "s = ",
      .format_scale_factor(scale_factors),
      " for C",
      seq_along(scale_factors),
      collapse = ", "
    )

    scale_description <- if (isTRUE(per_difference) && isTRUE(show_coef)) {
      paste0(
        "Under this coding, each model coefficient equals a scale factor s ",
        "multiplied by the difference displayed in the upper panel (",
        scale_list,
        "); the lower panel reports the coefficients divided by s, so that ",
        "each value refers to the displayed difference, whereas test ",
        "statistics, p-values, and Cohen's d are unaffected by this rescaling."
      )
    } else {
      paste0(
        "Under this coding, each model coefficient equals a scale factor s ",
        "multiplied by the difference displayed in the upper panel on the ",
        "linear-predictor scale (",
        scale_list,
        ")."
      )
    }
  }

  # Relation between the intended and the estimated comparisons. The
  # figure always shows the comparisons that the coefficients estimate.
  contrast_labels <- function(contrasts) {
    numbers <- sub("^Contrast_", "C", contrasts)
    if (length(numbers) == 1L) numbers else paste(paste(numbers[-length(numbers)], collapse = ", "), "and", numbers[length(numbers)])
  }

  hypothesis_description <- NULL

  if (!is.null(hypothesis_check) && nrow(hypothesis_check) > 0L) {

    failing <- hypothesis_check$Contrast[hypothesis_check$Status %in% c("differs", "reversed")]
    rescaled <- hypothesis_check$Contrast[hypothesis_check$Status == "rescaled"]

    hypothesis_description <- if (length(failing) == 0L) {
      paste0(
        "Each coefficient estimates the corresponding declared hypothesis",
        if (length(rescaled) > 0L) paste0(", ", contrast_labels(rescaled), " up to a positive factor"),
        "."
      )
    } else {
      paste0(
        contrast_labels(failing),
        if (length(failing) == 1L) " does" else " do",
        " not estimate the declared ",
        if (length(failing) == 1L) "hypothesis" else "hypotheses",
        "; the figure shows the comparisons that the coefficients estimate."
      )
    }

  } else if (length(coding_notes) > 0L) {

    noted <- sub(":.*$", "", coding_notes)

    hypothesis_description <- paste0(
      "The signs of the coding ",
      if (length(noted) == 1L) "column of " else "columns of ",
      contrast_labels(noted),
      " do not describe the comparisons that the coefficients estimate; ",
      "the figure shows the estimated comparisons."
    )
  }

  introduction <- paste0(
    "This plot summarizes the contrasts for the categorical predictor ",
    predictor,
    " in the fitted model for ",
    dependent_variable,
    "."
  )

  caption_parts <- c(
    introduction,
    contrast_description,
    hypothesis_description,
    upper_description,
    weighted_description,
    abbreviation_description,
    label_description,
    raw_description,
    coefficient_description,
    scale_description
  )

  caption_parts <- caption_parts[
    !is.na(caption_parts) &
      nzchar(caption_parts)
  ]

  paste(
    caption_parts,
    collapse = " "
  )
}


#' Build the lower panel of a one-predictor contrast plot
#'
#' Internal helper used by `moplot_onepred()`.
#'
#' @param prep An object returned by `.prepare_moplot_onepred()`.
#'
#' @param style Character string specifying the graphical style. Available
#'   options are `"minimal"` and `"publication"`.
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
.plot_moplot_onepred_lower <- function(
    prep,
    style = c("minimal", "publication"),
    digits = 2L,
    base_size = 13
) {

  # ==========================================================================
  # 1. ARGUMENT MATCHING AND VALIDATION
  # ==========================================================================

  style <- match.arg(style)

  if (
    !inherits(prep, "moplot_onepred_data") ||
    !is.list(prep)
  ) {
    stop(
      paste0(
        "`prep` must be an object returned by ",
        "`.prepare_moplot_onepred()`."
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

  parameter_data <- prep$parameter_data
  segment_data <- prep$segment_data

  axis_breaks <- as.numeric(prep$axis_breaks)

  required_parameter_columns <- c(
    "Contrast",
    "Contrast_number",
    "x",
    "Display_estimate",
    "Display_CI_low",
    "Display_CI_high",
    "Display_type",
    "d"
  )

  missing_parameter_columns <- setdiff(required_parameter_columns, names(parameter_data))

  if (length(missing_parameter_columns) > 0L) {
    stop(
      paste0(
        "`prep$parameter_data` is missing the following columns: ",
        paste(
          missing_parameter_columns,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  required_segment_columns <- c(
    "Contrast",
    "x_mid",
    "Significant"
  )

  missing_segment_columns <- setdiff(
    required_segment_columns,
    names(segment_data)
  )

  if (length(missing_segment_columns) > 0L) {
    stop(
      paste0(
        "`prep$segment_data` is missing the following columns: ",
        paste(
          missing_segment_columns,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  if (
    length(axis_breaks) == 0L ||
    anyNA(axis_breaks) ||
    any(!is.finite(axis_breaks))
  ) {
    stop(
      "The horizontal-axis information in `prep` is invalid.",
      call. = FALSE
    )
  }

  parameter_data$x <- as.numeric(parameter_data$x)

  parameter_data$Display_estimate <- as.numeric(parameter_data$Display_estimate)

  parameter_data$Display_CI_low <- as.numeric(parameter_data$Display_CI_low)

  parameter_data$Display_CI_high <- as.numeric(parameter_data$Display_CI_high)

  parameter_data$d <- as.numeric(parameter_data$d)

  segment_data$x_mid <- as.numeric(segment_data$x_mid)

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
  # 2. VERIFY ALIGNMENT WITH THE UPPER-PANEL MIDPOINTS
  # ==========================================================================

  upper_midpoints <- segment_data$x_mid

  if (
    length(upper_midpoints) != nrow(parameter_data) ||
    !isTRUE(
      all.equal(
        parameter_data$x,
        upper_midpoints,
        tolerance = 1e-10
      )
    )
  ) {
    stop(
      paste0(
        "The parameter positions are not aligned with the contrast ",
        "midpoints defined by the upper panel."
      ),
      call. = FALSE
    )
  }

  if (
    nrow(parameter_data) != nrow(segment_data) ||
    !identical(
      parameter_data$Contrast,
      segment_data$Contrast
    )
  ) {
    stop(
      paste0(
        "The lower-panel parameters could not be aligned with the ",
        "upper-panel contrasts."
      ),
      call. = FALSE
    )
  }

  # ==========================================================================
  # 3. ESTIMATE ABBREVIATIONS
  # ==========================================================================

  parameter_data <- .moplot_parameter_labels(
    parameter_data = parameter_data,
    per_difference = prep$per_difference,
    digits = digits
  )

  parameter_data$Contrast_label <- paste0(
    "Coef. ",
    parameter_data$Contrast_number
  )

  # ==========================================================================
  # 4. STYLE SETTINGS
  # ==========================================================================

  grayscale_palette <- c(
    midpoint = "#252525",
    segment = "#525252",
    interval = "#969696",
    guide = "#BDBDBD",
    text = "#252525",
    label_border = "#BDBDBD",
    label_fill = "white"
  )

  publication_palette <- c(
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

    parameter_data$Point_colour <- palette[["midpoint"]]

  } else {

    palette <- publication_palette

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
    is.finite(displayed_y)
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

  x_limits <- c(min(axis_breaks) - 0.5,max(axis_breaks) + 0.5)

  # ==========================================================================
  # 6. LABEL POSITIONS
  # ==========================================================================

  parameter_data$Label_x <-parameter_data$x + label_x_offset

  parameter_data$Label_y <-parameter_data$Display_estimate

  # ==========================================================================
  # 7. CREATE PLOT
  # ==========================================================================

  plot <- ggplot2::ggplot()

  plot <- plot +
    ggplot2::geom_vline(
      xintercept = upper_midpoints,
      linewidth = guide_linewidth,
      linetype = "dotted",
      colour = palette[["guide"]]
    )

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

  # Only the X markers inherit the upper-panel contrast colours.
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
    ) +
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
      labels = parameter_data$Contrast_label
    ) +
    ggplot2::coord_cartesian(
      xlim = x_limits,
      ylim = c(
        y_lower,
        y_upper
      ),
      clip = "off"
    ) +
    ggplot2::labs(
      x = "Contrast",
      y = if (
        identical(
          prep$coefficient_scale,
          "model"
        )
      ) {
        if (isTRUE(prep$per_difference)) "Model estimate / s" else "Model estimate"
      } else {
        unique(
          parameter_data$Display_type
        )[1L]
      }
    ) +
    ggplot2::theme_light(
      base_size = base_size
    ) +
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
            identical(
              style,
              "publication"
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


#' Compute recommended export dimensions for a one-predictor MoPlot
#'
#' Internal helper used by `moplot_onepred()`. Width grows with the number
#' of displayed quantities or, for polynomial contrasts, with the number of
#' component panels and of levels within each panel. Height reserves room
#' for the lower panel according to `panel_heights`.
#'
#' @param prep An object returned by `.prepare_moplot_onepred()`.
#' @param show_coef Logical. Whether the lower coefficient panel is present.
#' @param panel_heights Relative heights of the upper and lower panels.
#'
#' @return A named numeric vector containing recommended width and height
#'   in inches.
#'
#' @keywords internal
#' @noRd
.moplot_onepred_export_size <- function(
    prep,
    show_coef = TRUE,
    panel_heights = c(2, 1)
) {

  number_positions <- length(prep$axis_breaks)

  if (identical(prep$contrast_type, "Polynomial")) {

    # One panel per component, each wide enough for its levels and for
    # the coefficient label beside the marker in the lower panel.
    number_components <- nrow(prep$polynomial_data$components)

    width <- 1 + number_components * (1 + 0.35 * number_positions)

  } else {

    width <- 3.5 + 0.55 * number_positions
  }

  width <- max(6, width)

  upper_height <- 4.5

  height <- if (isTRUE(show_coef)) {
    upper_height * sum(panel_heights) / panel_heights[1L]
  } else {
    upper_height
  }

  c(
    width = round(width, 1),
    height = round(height, 1)
  )
}
