#' Plot contrasts for interacting categorical predictors
#'
#' @description
#' `moplot_multpred()` creates a graphical representation of model
#' interactions involving two or more categorical predictors.
#'
#' The function decomposes each selected interaction coefficient into the
#' elementary and higher-order contrasts implied by the contrast matrices used
#' to fit the model. The upper panel displays the observed or model-predicted
#' quantities involved in this decomposition. The optional lower panel displays
#' the corresponding model coefficient, its confidence interval, and an
#' approximate Cohen's d.
#'
#' By default, the function identifies the interaction of highest order in the
#' fitted model and displays all coefficients associated with that interaction.
#' If more than one interaction of the highest order is present, the function
#' asks the user to select one explicitly through `interaction`.
#'
#' @param model A fitted model supported by the package. Supported model
#'   classes include objects returned by [stats::lm()], [stats::glm()],
#'   [lme4::lmer()], [lme4::glmer()], [nlme::lme()], and
#'   [glmmTMB::glmmTMB()].
#'
#' @param data Optional data frame containing the variables used to fit
#'   `model`. If `NULL`, the data supplied in the model call are recovered,
#'   restricted to the rows used for fitting; the model frame is used only when
#'   they are unavailable. The original variables are needed for model-based
#'   estimates when the model contains transformed terms such as `log(x)`.
#'   Rows of a supplied data frame that were excluded from the fit, through
#'   `subset` or missing values, are dropped, matching rows by name.
#'
#' @param interaction Optional character string identifying the interaction
#'   term to display, for example `"A:B:C"`.
#'
#'   If `NULL`, the function automatically selects the unique interaction of
#'   highest order in the fitted model. If several interactions share the
#'   highest order, the function stops and asks the user to specify
#'   `interaction`.
#'
#' @param coefficients Optional character vector identifying specific model
#'   coefficients to display.
#'
#'   If `NULL`, all coefficients associated with the selected `interaction`
#'   are displayed. If `interaction = NULL` and `coefficients` is supplied,
#'   the corresponding interaction is inferred automatically, provided that
#'   all selected coefficients belong to the same interaction.
#'
#' @param interaction_info Logical. Should the function print the available
#'   interaction terms and their corresponding coefficient selections?
#'   The default is `TRUE`.
#'
#'   The printed table shows syntax that can be copied directly into
#'   `interaction` and `coefficients` in subsequent calls.
#'
#' @param split_rule Character string specifying how the top-level split
#'   variable is selected during the hierarchical decomposition. Available
#'   options are:
#'
#'   \describe{
#'     \item{`"last_variable"`}{
#'       Selects the last variable appearing in the interaction coefficient.
#'     }
#'     \item{`"first_variable"`}{
#'       Selects the first variable appearing in the interaction coefficient.
#'     }
#'     \item{`"lowest_complexity"`}{
#'       Selects the variable with the lowest contrast-complexity score.
#'     }
#'     \item{`"highest_complexity"`}{
#'       Selects the variable with the highest contrast-complexity score.
#'     }
#'   }
#'
#' @param top_split Optional character string identifying the predictor to
#'   force as the top-level split variable. When supplied, `top_split`
#'   overrides `split_rule`.
#'
#' @param sig_level Numeric value between 0 and 1 defining the significance
#'   threshold used for coefficient labels. The default is `0.05`.
#'
#' @param conf_level Numeric value between 0 and 1 specifying the confidence
#'   level used for displayed confidence intervals. The default is `0.95`.
#'
#' @param estimate_type Character string specifying how quantities in the
#'   upper panel are computed. Available options are:
#'
#'   \describe{
#'     \item{`"auto"`}{
#'       Uses observed means when the selected interaction involves all the
#'       predictors of the fixed part, the model has no random effects and no
#'       unequal prior weights, and, for generalized models,
#'       `prediction_scale = "response"`. Model-based predictions are used
#'       otherwise.
#'     }
#'     \item{`"observed"`}{
#'       Displays observed aggregated means computed directly from `data`.
#'     }
#'     \item{`"predicted"`}{
#'       Displays model-based predictions computed from the fixed effects and
#'       their variance-covariance matrix.
#'     }
#'   }
#'
#'   The default is `"auto"`.
#'
#' @param prediction_scale Character string specifying the scale used for
#'   model-based predictions when `estimate_type = "predicted"`.
#'   Available options are:
#'
#'   \describe{
#'     \item{`"response"`}{
#'       Displays predictions on the response scale.
#'     }
#'     \item{`"link"`}{
#'       Displays predictions on the linear-predictor/link scale.
#'     }
#'   }
#'
#'   The default is `"response"`. This argument is ignored when
#'   `estimate_type = "observed"`.
#'
#' @param coefficient_scale Character string specifying the quantity displayed
#'   in the lower panel. Available options are:
#'
#'   \describe{
#'     \item{`"model"`}{
#'       Displays the coefficient on its original model scale.
#'     }
#'     \item{`"interpreted"`}{
#'       Displays the interpreted coefficient when an appropriate
#'       transformation is available, such as an odds ratio or incidence rate
#'       ratio.
#'     }
#'   }
#'
#'   The default is `"model"`.
#'
#' @param per_difference Logical. Should the lower-panel coefficients be
#'   divided by their scale factor? The default is `FALSE`.
#'
#'   Each interaction coefficient equals the product of the scale factors of
#'   its components multiplied by the displayed difference of differences
#'   (for example `s = 1/3` when a Helmert contrast with `s = 1/3` interacts
#'   with a sum contrast). When `per_difference = TRUE`, the lower panel
#'   reports `beta / s` and its confidence interval; with
#'   `coefficient_scale = "interpreted"`, ratios are computed as
#'   `exp(beta / s)`. Test statistics, p-values, and Cohen's d are unaffected.
#'
#' @param hypotheses Optional named list declaring, for one or more
#'   categorical predictors, the comparisons that their contrast coefficients
#'   are intended to estimate, one matrix per predictor. Each interaction
#'   coefficient is checked through its components; the result is reported in
#'   a message and in the caption. The figure always shows the comparisons
#'   that the coefficients estimate. See [contrasts_multpred()].
#'
#' @param bounds Optional numeric vector `c(lower, upper)` giving the natural
#'   range of a response analyzed with a Gaussian (identity link) model, such
#'   as `c(1, 7)` for a Likert score or `c(0, Inf)` for a reaction time, so that
#'   confidence intervals remain within it. Ignored, with a warning, for
#'   generalized models. See [contrasts_multpred()].
#'
#' @param show_coef Logical. Should the lower coefficient panel be displayed?
#'   The default is `TRUE`.
#'
#' @param style Character string specifying the graphical representation of the
#'   interaction decomposition. Available options are:
#'
#'   \describe{
#'     \item{`"minimal"`}{
#'       Displays the observed or predicted quantities and the dashed segments
#'       representing the elementary component contrasts, without displaying
#'       intermediate delta formulas.
#'     }
#'     \item{`"tree"`}{
#'       Displays the decomposition as a hierarchical tree showing how
#'       elementary contrasts are recursively combined into higher-order
#'       contrasts and ultimately into the final interaction contrast.
#'     }
#'   }
#'
#'   The default is `"minimal"`.
#'
#' @param caption Controls whether a descriptive caption is printed in the
#'   console. It can be:
#'
#'   \itemize{
#'     \item `TRUE`, to generate and print an automatic caption;
#'     \item `FALSE` or `NULL`, to suppress the caption;
#'     \item a single character string to print a custom caption.
#'   }
#'
#'   The automatic caption describes the selected interaction, the graphical
#'   hierarchy, the contrast coding assigned to each interacting variable, the
#'   displayed estimate type, the selected graphical style, and the lower
#'   coefficient panel when present.
#'
#'   The caption is printed in the console and is not added to the graphical
#'   object.
#'
#' @param digits Non-negative integer specifying the number of decimal places
#'   used in lower-panel estimate and effect-size labels. The default is `2`.
#'
#' @param base_size Positive numeric value defining the base font size used in
#'   both panels. The default is `13`.
#'
#' @param panel_heights Numeric vector of length two specifying the relative
#'   heights of the upper and lower panels. The default is `c(2, 1)`.
#'   This argument is used only when `show_coef = TRUE`.
#'
#' @details
#' `moplot_multpred()` is a graphical interface to the hierarchical
#' decomposition performed by [contrasts_multpred()].
#'
#' Each selected interaction coefficient is decomposed according to the
#' contrast matrices assigned to the categorical predictors involved in the
#' interaction. This allows treatment, sum, Helmert, reverse Helmert,
#' sliding-difference, and interpretable custom contrasts to be represented
#' directly in terms of the factor levels or level aggregates implied by the
#' fitted model.
#'
#' The graphical hierarchy is determined by `split_rule` or `top_split`.
#' The top-level split variable defines the outermost grouping structure,
#' whereas successive predictors define progressively lower levels of the
#' decomposition. The final predictor in the resulting hierarchy defines the
#' quantities represented on the x-axis.
#'
#' When an interaction contains more than two predictors, the decomposition is
#' recursive. Elementary contrasts are first constructed from pairs of
#' displayed quantities and are then combined into progressively higher-order
#' contrasts until the complete model interaction contrast is obtained.
#'
#' For `estimate_type = "predicted"`, predictions are population-level
#' predictions based on the fixed effects. For generalized models on the
#' response scale, uncertainty is propagated using the fixed-effect
#' variance-covariance matrix and the multivariate delta method where
#' appropriate. A methodological note describing only the settings relevant to
#' the current call is printed when the function is run.
#'
#' @section Interaction and coefficient selection:
#' The arguments `interaction` and `coefficients` provide two complementary
#' levels of control.
#'
#' To display all coefficients associated with one interaction, specify only
#' `interaction`, for example:
#'
#' \preformatted{
#' interaction = "A:B:C"
#' }
#'
#' To display only selected coefficients from that interaction, specify both:
#'
#' \preformatted{
#' interaction = "A:B:C"
#' coefficients = c("A1:B2:C2", "A2:B2:C2")
#' }
#'
#' If neither argument is supplied, the unique highest-order interaction is
#' selected automatically.
#'
#' When `interaction_info = TRUE`, the function prints the interactions found
#' in the model together with valid coefficient selections, providing
#' copy-ready syntax for subsequent calls.
#'
#' @section Graphical styles:
#' The three graphical styles differ in how much of the hierarchical
#' decomposition is shown.
#'
#' `"minimal"` emphasizes the observed or predicted quantities and elementary
#' contrasts.
#'
#' `"tree"` represents the same hierarchy through branching segments,
#' emphasizing the recursive structure of the decomposition.
#'
#' The lower coefficient panel, when requested, remains horizontally aligned
#' with the corresponding contrast block in the upper panel for all three
#' styles.
#'
#' @return
#' When `show_coef = FALSE`, a `ggplot` object containing the upper interaction
#' decomposition panel.
#'
#' When `show_coef = TRUE`, a `patchwork` object containing the aligned upper and
#' lower panels.
#'
#' @examples
#' set.seed(123)
#'
#' dat <- expand.grid(
#'   A = factor(c("A1", "A2", "A3")),
#'   B = factor(c("B1", "B2")),
#'   C = factor(c("C1", "C2")),
#'   rep = seq_len(20)
#' )
#'
#' contrasts(dat$A) <- reverse_helmert(3)
#' contrasts(dat$B) <- MASS::contr.sdif(2)
#' contrasts(dat$C) <- stats::contr.treatment(2)
#'
#' dat$y <- with(
#'   dat,
#'   10 +
#'     2 * (A == "A2") +
#'     4 * (A == "A3") +
#'     2 * (B == "B2") +
#'     2 * (C == "C2") +
#'     2 * (A == "A2" & B == "B2" & C == "C2") +
#'     rnorm(nrow(dat))
#' )
#'
#' mod <- stats::lm(
#'   y ~ A * B * C,
#'   data = dat
#' )
#'
#' # Automatically display the unique highest-order interaction.
#' plot_auto <- moplot_multpred(
#'   model = mod,
#'   data = dat
#' )
#'
#' # Display all coefficients from a specific interaction.
#' plot_tree <- moplot_multpred(
#'   model = mod,
#'   interaction = "A:B:C",
#'   style = "tree"
#' )
#'
#' # Display one specific coefficient using the minimal representation.
#' plot_minimal <- moplot_multpred(
#'   model = mod,
#'   interaction = "A:B:C",
#'   coefficients = "A1:B2-1:C2",
#'   style = "minimal"
#' )
#'
#' # Use model-based predictions instead of observed means.
#' plot_predicted <- moplot_multpred(
#'   model = mod,
#'   interaction = "A:B:C",
#'   estimate_type = "predicted",
#'   prediction_scale = "response"
#' )
#'
#' # The recommended size is stored in the plot and used by moplot_export().
#' attr(plot_tree, "moplot_export_size")
#'
#' path <- moplot_export(
#'   plot_tree,
#'   filename = file.path(tempdir(), "moplot_multpred.png"),
#'   dpi = 72
#' )
#'
#' unlink(path)
#'
#' @importFrom patchwork wrap_plots
#' @export
moplot_multpred <- function(
    model,
    data=NULL,
    interaction = NULL,
    coefficients = NULL,
    interaction_info = TRUE,
    split_rule = c(
      "last_variable",
      "first_variable",
      "lowest_complexity",
      "highest_complexity"
    ),
    top_split = NULL,
    sig_level = 0.05,
    conf_level = 0.95,
    estimate_type = c("auto", "observed", "predicted"),
    prediction_scale = c("response", "link"),
    coefficient_scale = c("model", "interpreted"),
    per_difference = FALSE,
    bounds = NULL,
    hypotheses = NULL,
    show_coef = TRUE,
    style = c("minimal", "tree"),
    caption = TRUE,
    digits = 2L,
    base_size = 13,
    panel_heights = c(2, 1)
) {

  # ============================================================
  # 1. ARGUMENT MATCHING
  # ============================================================

  split_rule <- match.arg(split_rule)

  estimate_type <- match.arg(estimate_type)

  prediction_scale <- match.arg(prediction_scale)

  coefficient_scale <- match.arg(coefficient_scale)

  style <- match.arg(style)

  if (
    !is.logical(per_difference) ||
    length(per_difference) != 1L ||
    is.na(per_difference)
  ) {
    stop(
      "`per_difference` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (
    !is.logical(show_coef) ||
    length(show_coef) != 1L ||
    is.na(show_coef)
  ) {
    stop(
      "`show_coef` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  # ============================================================
  # 1.1 RESOLVE INTERACTION AND COEFFICIENT SELECTION
  # ============================================================

  selection <- .resolve_moplot_interaction_selection(
    model = model,
    interaction = interaction,
    coefficients = coefficients,
    interaction_info = interaction_info
  )

  interaction <- selection$interaction
  coefficients <- selection$coefficients

  # ============================================================
  # 2. COMPUTE INTERACTION CONTRASTS
  # ============================================================

  contrast_output <- contrasts_multpred(
    model = model,
    data = data,
    coefficients = coefficients,
    split_rule = split_rule,
    top_split = top_split,
    conf_level = conf_level,
    sig_level = sig_level,
    configuration = "all",
    estimate_type = estimate_type,
    prediction_scale = prediction_scale,
    methodological_note = FALSE,
    bounds = bounds,
    hypotheses = hypotheses
  )

  # estimate_type = "auto" is resolved by contrasts_multpred().
  estimate_type <- attr(contrast_output, "estimate_type")

  # ============================================================
  # 2.1 METHODOLOGICAL NOTE
  # ============================================================

  plot_methodological_note <- .build_multpred_methodological_note(
    model = model,
    estimate_type = estimate_type,
    prediction_scale = prediction_scale,
    split_rule = split_rule,
    top_split = top_split,
    configuration = NULL
  )

  message(plot_methodological_note)

  # ============================================================
  # 3. PREPARE UPPER-PANEL DATA
  # ============================================================

  prep <- .prepare_moplot_multpred_data(contrast_output)

  label_key <- .collect_moplot_label_key(prep = prep, max_chars = 5L)

  # ============================================================
  # 3.1 CAPTION
  # ============================================================

  if (isTRUE(caption)) {

    caption_text <- .moplot_multpred_caption(
      model = model,
      contrast_output = contrast_output,
      prep = prep,
      interaction = interaction,
      style = style,
      sig_level = sig_level,
      show_coef = show_coef,
      coefficient_scale = coefficient_scale,
      per_difference = per_difference,
      estimate_type = estimate_type,
      prediction_scale = prediction_scale,
      label_key=label_key
    )

    cat(
      "Plot caption:\n\n",
      paste0(
        strwrap(
          caption_text,
          width = 100
        ),
        collapse = "\n"
      ),
      "\n\n",
      sep = ""
    )

  } else if (
    is.character(caption) &&
    length(caption) == 1L &&
    !is.na(caption)
  ) {

    cat(
      "Plot caption:\n\n",
      paste0(
        strwrap(
          caption,
          width = 100
        ),
        collapse = "\n"
      ),
      "\n\n",
      sep = ""
    )
  }

  # ============================================================
  # 4. BUILD UPPER PANEL
  # ============================================================

  upper_plot <- .plot_moplot_multpred_upper(
    prep = prep,
    style = style,
    base_size = base_size,
    label_key=label_key,
    y_title = .outcome_name(model)
  )


  export_size <- .moplot_multpred_export_size(
    prep = prep,
    style = style,
    show_coef = show_coef
  )

  # ============================================================
  # 5. RETURN UPPER PANEL ONLY
  # ============================================================

  if (!isTRUE(show_coef)) {

    attr(
      upper_plot,
      "moplot_export_size"
    ) <- export_size

    return(
      upper_plot
    )
  }

  # ============================================================
  # 6. PREPARE MODEL PARAMETERS
  # ============================================================

  parameters <- .prepare_moplot_multpred_parameters(
    model = model,
    coefficient_data = prep$coefficient_data,
    conf_level = conf_level,
    coefficient_scale = coefficient_scale,
    scale_factors = attr(contrast_output, "scale_factors"),
    per_difference = per_difference
  )

  # ============================================================
  # 7. BUILD LOWER PANEL
  # ============================================================

  lower_plot <- .plot_moplot_multpred_lower(
    prep = prep,
    parameters = parameters,
    style = style,
    digits = digits,
    base_size = base_size
  )

  # ============================================================
  # 8. ASSEMBLE PANELS
  # ============================================================

  final_plot <- .assemble_moplot_multpred(
    upper_plot = upper_plot,
    lower_plot = lower_plot,
    panel_heights = panel_heights
  )

  attr(
    final_plot,
    "moplot_export_size"
  ) <- export_size

  final_plot
}
