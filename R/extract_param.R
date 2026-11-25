#' Extract and interpret model parameters
#'
#' @description
#' `extract_param()` extracts model coefficients and inferential statistics
#' using `parameters::model_parameters()`. It returns a compact,
#' non-redundant parameter table and, when appropriate, transforms coefficients
#' and confidence intervals into a more directly interpretable metric.
#'
#' Depending on the model family and link function, the interpreted estimate
#' may represent an odds ratio, risk ratio, incidence rate ratio, hazard
#' ratio, or another multiplicative ratio. For models for which no meaningful automatic
#' transformation is available, the original coefficient is retained.
#'
#' The function also computes an approximate Cohen's d for each coefficient:
#' from the log odds ratio for logistic models, and from the coefficient-level
#' t or z statistic otherwise.
#'
#' @param model A fitted statistical model supported by
#'   `parameters::model_parameters()`.
#'
#' @param conf_level A single numeric value between 0 and 1 specifying the confidence
#'   level. The default is `0.95`.
#'
#' @param include_intercept Logical. Should the intercept be included in the
#'   returned table? The default is `FALSE`.
#'
#' @param effects Character string specifying the model effects to extract for
#'   models that distinguish between fixed and random effects. The default is
#'   `"fixed"`.
#'
#' @param component Character string specifying the model component to extract
#'   for models with multiple components, such as zero-inflated or hurdle
#'   models. The default is `"conditional"`.
#'
#' @param transform Character string specifying how model coefficients should
#'   be transformed for interpretation. Possible values are:
#'
#'   \describe{
#'     \item{`"auto"`}{Automatically exponentiate coefficients from models
#'       with log or logit links and from Cox proportional-hazards models.}
#'     \item{`"none"`}{Retain the original coefficient and confidence interval
#'       as the interpreted estimate.}
#'     \item{`"exponentiate"`}{Exponentiate the coefficient and confidence
#'       interval regardless of the detected model family or link.}
#'   }
#'
#'   The default is `"auto"`.
#'
#' @details
#' The `Estimate` column always contains the coefficient on the scale on which
#' the model was fitted. For generalized models, this is generally the scale of
#' the linear predictor.
#'
#' When `transform = "auto"`, the following transformations are applied:
#'
#' \describe{
#'   \item{Identity link}{No transformation. `Interpreted_Estimate` is equal
#'     to `Estimate`.}
#'   \item{Binomial model with logit link}{The coefficient is exponentiated
#'     and interpreted as an odds ratio.}
#'   \item{Binomial model with log link}{The coefficient is exponentiated and
#'     interpreted as a risk ratio.}
#'   \item{Poisson or negative-binomial model with log link}{The coefficient
#'     is exponentiated and interpreted as an incidence rate ratio.}
#'   \item{Other models with a log link}{The coefficient is exponentiated and
#'     interpreted as a multiplicative ratio.}
#'   \item{Cox proportional-hazards model}{The coefficient is exponentiated
#'     and interpreted as a hazard ratio.}
#' }
#'
#' No automatic exponentiation is performed for links such as `probit`,
#' `cauchit`, or `cloglog`, because exponentiating an individual coefficient
#' does not generally produce a standard effect measure with a direct
#' interpretation.
#'
#' The interpreted estimates are not necessarily predictions on the response
#' scale. For example, an odds ratio is a ratio of odds rather than a predicted
#' probability. Predicted probabilities or expected responses require the
#' intercept and values of the other model predictors.
#'
#' Confidence intervals and p-values use the same inferential method as the
#' significance labels of the plotting functions, so that an interval excluding
#' zero always corresponds to a significant coefficient:
#'
#' \itemize{
#'   \item linear mixed models fitted with lme4 or lmerTest use Satterthwaite
#'     degrees of freedom when lmerTest is installed, as in their p-values;
#'   \item generalized linear (mixed) models use Wald intervals, which
#'     correspond to the Wald z tests reported by `summary()`; profile
#'     intervals can disagree with these tests near the significance threshold;
#'   \item other models use the default method of
#'     `parameters::model_parameters()`.
#' }
#'
#' Approximate Cohen's d values are computed as follows, as reported in
#' `d_method`:
#'
#' \itemize{
#'   \item for models with a logit link, the log odds ratio is converted as
#'     \eqn{d = \beta \sqrt{3} / \pi} (Chinn, 2000), which treats the
#'     logistic distribution of the latent response as the reference
#'     standard deviation; the same conversion is applied to the confidence
#'     limits;
#'   \item t statistics are converted using `effectsize::t_to_d()` with the
#'     coefficient-specific degrees of freedom;
#'   \item z statistics are converted using `effectsize::z_to_d()` with the
#'     number of observations returned by `stats::nobs()`.
#' }
#'
#' The conversions based on t or z do not depend on how the contrast is
#' scaled. The conversion from the log odds ratio refers to the coefficient
#' as coded: for a contrast whose coefficient equals a scale factor \eqn{s}
#' multiplied by the compared difference (for example Helmert contrasts; see
#' [contrasts_onepred()]), the d of the compared difference is \eqn{d / s}.
#' The plotting functions apply this rescaling.
#'
#' These conversions are approximations and assume conditions implicit in the
#' conversion formulas. In particular, the conversion to Cohen's d assumes two
#' groups of approximately equal size, and the conversion from the log odds
#' ratio assumes a logistic latent response. Directly estimated effect sizes should
#' be preferred when the necessary data and an appropriate estimator are
#' available.
#'
#' @return
#' A data frame with one row for each extracted model parameter and the
#' following columns:
#'
#' \describe{
#'   \item{Parameter}{Name of the model parameter.}
#'   \item{Estimate}{Untransformed model coefficient on the model scale.}
#'   \item{SE}{Standard error of the untransformed coefficient.}
#'   \item{CI}{Confidence level used for the interval.}
#'   \item{CI_low}{Lower confidence limit of `Estimate`.}
#'   \item{CI_high}{Upper confidence limit of `Estimate`.}
#'   \item{Statistic}{Coefficient-level test statistic.}
#'   \item{Statistic_type}{Type of test statistic, usually `"t"` or `"z"`.}
#'   \item{df}{Coefficient-level degrees of freedom, when available.}
#'   \item{p}{Coefficient-level p-value, when available.}
#'   \item{Interpreted_Estimate}{Coefficient transformed into the
#'     interpretation indicated by `Estimate_type`.}
#'   \item{Interpreted_CI_low}{Lower confidence limit of the interpreted
#'     estimate.}
#'   \item{Interpreted_CI_high}{Upper confidence limit of the interpreted
#'     estimate.}
#'   \item{Estimate_type}{Interpretation of `Interpreted_Estimate`, such as
#'     `"Coefficient"`, `"Odds ratio"`, `"Risk ratio"`,
#'     `"Incidence rate ratio"`, `"Hazard ratio"`, or `"Multiplicative ratio"`.}
#'   \item{d}{Approximate Cohen's d.}
#'   \item{d_CI_low}{Lower confidence limit of approximate Cohen's d.}
#'   \item{d_CI_high}{Upper confidence limit of approximate Cohen's d.}
#'   \item{d_method}{Conversion used for `d`: `"log odds ratio"`, `"t"`, or
#'     `"z"`; `NA` when no conversion was possible.}
#' }
#'
#' @note
#' The coefficient transformation is based on the model family and link
#' detected through `stats::family()` whenever available. Some model classes
#' may use non-standard structures or parameterizations and may therefore
#' require `transform = "none"` or `transform = "exponentiate"` to override
#' automatic behaviour.
#'
#' In models with multiple components, the interpretation may depend on the
#' component selected through `component`.
#'
#' @examples
#' set.seed(123)
#'
#' # Linear regression
#' dat_lm <- data.frame(
#'   y = rnorm(120),
#'   x = rnorm(120),
#'   group = factor(rep(c("A", "B", "C"), each = 40))
#' )
#'
#' contrasts(dat_lm$group) <- contr.sum(3)
#'
#' mod_lm <- lm(y ~ x + group, data = dat_lm)
#'
#' extract_param(mod_lm)
#'
#' # Binomial logistic regression
#' dat_glm <- data.frame(
#'   y = rbinom(200, size = 1, prob = 0.40),
#'   x = rnorm(200),
#'   group = factor(sample(c("A", "B", "C"), 200, replace = TRUE))
#' )
#'
#' contrasts(dat_glm$group) <- contr.sum(3)
#'
#' mod_glm <- glm(
#'   y ~ x + group,
#'   data = dat_glm,
#'   family = binomial(link = "logit")
#' )
#'
#' extract_param(mod_glm)
#'
#' # Disable automatic coefficient transformation
#' extract_param(mod_glm, transform = "none")
#'
#' @importFrom parameters model_parameters
#' @references
#' Chinn, S. (2000). A simple method for converting an odds ratio to effect
#' size for use in meta-analysis. *Statistics in Medicine, 19*(22),
#' 3127-3131.
#'
#' @importFrom effectsize t_to_d z_to_d
#' @importFrom stats family nobs
#' @export

extract_param <- function(
    model,
    conf_level = 0.95,
    include_intercept = FALSE,
    effects = "fixed",
    component = "conditional",
    transform = c("auto", "none", "exponentiate")
) {

  # ==========================================================================
  # Argument checks
  # ==========================================================================

  transform <- match.arg(transform)

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

  if (
    !is.logical(include_intercept) ||
    length(include_intercept) != 1L ||
    is.na(include_intercept)
  ) {
    stop(
      "`include_intercept` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (
    !is.character(effects) ||
    length(effects) != 1L ||
    is.na(effects) ||
    !nzchar(effects)
  ) {
    stop(
      "`effects` must be a single non-empty character string.",
      call. = FALSE
    )
  }

  if (
    !is.character(component) ||
    length(component) != 1L ||
    is.na(component) ||
    !nzchar(component)
  ) {
    stop(
      "`component` must be a single non-empty character string.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # Internal helpers
  # ==========================================================================

  first_existing_column <- function(data, candidates) {

    available <- candidates[candidates %in% names(data)]

    if (length(available) == 0L) {
      return(NULL)
    }

    available[1L]
  }

  extract_numeric_column <- function(data, column, n) {

    if (is.null(column)) {
      return(rep(NA_real_, n))
    }

    suppressWarnings(
      as.numeric(data[[column]])
    )
  }

  remove_intercept <- function(data) {

    if (
      !include_intercept &&
      "Parameter" %in% names(data)
    ) {
      data <- data[
        data$Parameter != "(Intercept)",
        ,
        drop = FALSE
      ]
    }

    data
  }

  # ==========================================================================
  # Detect model family and link
  # ==========================================================================

  family_object <- tryCatch(
    stats::family(model),
    error = function(e) NULL
  )

  family_name <- if (
    !is.null(family_object) &&
    !is.null(family_object$family)
  ) {
    tolower(as.character(family_object$family)[1L])
  } else {
    NA_character_
  }

  link_name <- if (
    !is.null(family_object) &&
    !is.null(family_object$link)
  ) {
    tolower(as.character(family_object$link)[1L])
  } else {
    NA_character_
  }

  is_cox_model <- inherits(model, "coxph")

  # ==========================================================================
  # Determine transformation and interpretation
  # ==========================================================================

  exponentiate_estimate <- switch(
    transform,
    none = FALSE,
    exponentiate = TRUE,
    auto = {
      is_cox_model ||
        isTRUE(link_name %in% c("log", "logit"))
    }
  )

  if (!exponentiate_estimate) {

    estimate_type <- "Coefficient"

  } else if (is_cox_model) {

    estimate_type <- "Hazard ratio"

  } else if (identical(link_name, "logit")) {

    estimate_type <- "Odds ratio"

  } else if (
    identical(link_name, "log") &&
    family_name %in% c("binomial", "quasibinomial")
  ) {

    estimate_type <- "Risk ratio"

  } else if (
    identical(link_name, "log") &&
    grepl(
      "poisson|negative binomial|negative-binomial|nbinom",
      family_name
    )
  ) {

    estimate_type <- "Incidence rate ratio"

  } else if (identical(link_name, "log")) {

    estimate_type <- "Multiplicative ratio"

  } else {

    estimate_type <- "Exponentiated coefficient"
  }

  # ==========================================================================
  # Extract untransformed parameters
  # ==========================================================================

  # Intervals and p-values follow the inferential method of the
  # significance labels used by the plotting functions (see
  # .extract_coef_table()): Satterthwaite degrees of freedom for linear
  # mixed models when lmerTest is available, and Wald intervals for
  # generalized models, matching their Wald z tests. Profile intervals,
  # the default of parameters for glm, can disagree with these tests near
  # the significance threshold.
  ci_method <- if (inherits(model, "lmerMod") &&
                   requireNamespace("lmerTest", quietly = TRUE)) {
    "satterthwaite"
  } else if (!is.na(family_name) &&
             !(family_name == "gaussian" && identical(link_name, "identity"))) {
    "wald"
  } else {
    NULL
  }

  parameter_arguments <- list(
    model,
    ci = conf_level,
    effects = effects,
    component = component,
    exponentiate = FALSE,
    pretty_names = FALSE,
    verbose = FALSE
  )

  if (!is.null(ci_method)) {
    parameter_arguments$ci_method <- ci_method
  }

  raw_parameters <- tryCatch(
    do.call(parameters::model_parameters, parameter_arguments),
    error = function(e) {
      stop(
        paste0(
          "Model parameters could not be extracted. ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  raw_parameters <- as.data.frame(raw_parameters)
  raw_parameters <- remove_intercept(raw_parameters)
  rownames(raw_parameters) <- NULL

  if (!"Parameter" %in% names(raw_parameters)) {
    stop(
      paste0(
        "The table returned by `parameters::model_parameters()` ",
        "does not contain a `Parameter` column."
      ),
      call. = FALSE
    )
  }

  n_parameters <- nrow(raw_parameters)

  # ==========================================================================
  # Identify columns returned by model_parameters()
  # ==========================================================================

  estimate_column <- first_existing_column(
    raw_parameters,
    c("Coefficient", "Estimate")
  )

  se_column <- first_existing_column(
    raw_parameters,
    c("SE", "Std_Error", "Std. Error")
  )

  ci_column <- first_existing_column(
    raw_parameters,
    "CI"
  )

  ci_low_column <- first_existing_column(
    raw_parameters,
    c("CI_low", "CI_Low", "Conf.Low")
  )

  ci_high_column <- first_existing_column(
    raw_parameters,
    c("CI_high", "CI_High", "Conf.High")
  )

  p_column <- first_existing_column(
    raw_parameters,
    c("p", "p_value", "P")
  )

  df_column <- first_existing_column(
    raw_parameters,
    c("df_error", "df", "df_residual")
  )

  # ==========================================================================
  # Identify test statistic and statistic type
  # ==========================================================================

  if ("t" %in% names(raw_parameters)) {

    statistic_column <- "t"
    statistic_type <- rep("t", n_parameters)

  } else if ("z" %in% names(raw_parameters)) {

    statistic_column <- "z"
    statistic_type <- rep("z", n_parameters)

  } else if ("Statistic" %in% names(raw_parameters)) {

    statistic_column <- "Statistic"
    statistic_type <- rep(NA_character_, n_parameters)

    statistic_attribute <- attr(
      raw_parameters$Statistic,
      "statistic"
    )

    if (
      !is.null(statistic_attribute) &&
      grepl("^t$", statistic_attribute, ignore.case = TRUE)
    ) {
      statistic_type[] <- "t"
    }

    if (
      !is.null(statistic_attribute) &&
      grepl("^z$", statistic_attribute, ignore.case = TRUE)
    ) {
      statistic_type[] <- "z"
    }

    if (all(is.na(statistic_type))) {

      if (!is.null(df_column)) {

        df_values <- suppressWarnings(
          as.numeric(raw_parameters[[df_column]])
        )

        finite_df <- is.finite(df_values)

        statistic_type[finite_df] <- "t"
        statistic_type[!finite_df] <- "z"

      } else {

        statistic_type[] <- "z"
      }
    }

  } else {

    statistic_column <- NULL
    statistic_type <- rep(NA_character_, n_parameters)
  }

  # ==========================================================================
  # Build the non-redundant output table
  # ==========================================================================

  output <- data.frame(
    Parameter = as.character(raw_parameters$Parameter),

    Estimate = extract_numeric_column(
      raw_parameters,
      estimate_column,
      n_parameters
    ),

    SE = extract_numeric_column(
      raw_parameters,
      se_column,
      n_parameters
    ),

    CI = if (!is.null(ci_column)) {
      extract_numeric_column(
        raw_parameters,
        ci_column,
        n_parameters
      )
    } else {
      rep(conf_level, n_parameters)
    },

    CI_low = extract_numeric_column(
      raw_parameters,
      ci_low_column,
      n_parameters
    ),

    CI_high = extract_numeric_column(
      raw_parameters,
      ci_high_column,
      n_parameters
    ),

    Statistic = extract_numeric_column(
      raw_parameters,
      statistic_column,
      n_parameters
    ),

    Statistic_type = statistic_type,

    df = extract_numeric_column(
      raw_parameters,
      df_column,
      n_parameters
    ),

    p = extract_numeric_column(
      raw_parameters,
      p_column,
      n_parameters
    ),

    stringsAsFactors = FALSE
  )

  # ==========================================================================
  # Extract interpreted estimates
  # ==========================================================================

  # Exponentiated estimates and limits are exp() of the coefficient and
  # of its limits; a second call to model_parameters() is not needed.
  if (exponentiate_estimate) {

    output$Interpreted_Estimate <- exp(output$Estimate)
    output$Interpreted_CI_low <- exp(output$CI_low)
    output$Interpreted_CI_high <- exp(output$CI_high)

  } else {

    output$Interpreted_Estimate <- output$Estimate
    output$Interpreted_CI_low <- output$CI_low
    output$Interpreted_CI_high <- output$CI_high
  }

  output$Estimate_type <- rep(
    estimate_type,
    nrow(output)
  )

  # ==========================================================================
  # Extract sample size for z-to-d conversions
  # ==========================================================================

  n_observations <- tryCatch(
    as.numeric(stats::nobs(model)),
    error = function(e) NA_real_
  )

  if (
    length(n_observations) != 1L ||
    is.na(n_observations) ||
    !is.finite(n_observations) ||
    n_observations <= 0
  ) {
    n_observations <- NA_real_
  }

  # ==========================================================================
  # Compute approximate Cohen's d
  # ==========================================================================

  output$d <- NA_real_
  output$d_CI_low <- NA_real_
  output$d_CI_high <- NA_real_
  output$d_method <- NA_character_

  # Logistic models: conversion of the log odds ratio (Chinn, 2000),
  # d = beta * sqrt(3) / pi, applied to the estimate and to its limits.
  if (identical(link_name, "logit")) {

    logit_to_d <- sqrt(3) / pi

    output$d <- output$Estimate * logit_to_d
    output$d_CI_low <- output$CI_low * logit_to_d
    output$d_CI_high <- output$CI_high * logit_to_d
    output$d_method[is.finite(output$d)] <- "log odds ratio"
  }

  for (i in seq_len(nrow(output))) {

    if (!is.na(output$d_method[i])) {
      next
    }

    statistic_value <- output$Statistic[i]
    current_statistic_type <- output$Statistic_type[i]
    df_value <- output$df[i]

    if (
      is.na(statistic_value) ||
      !is.finite(statistic_value) ||
      is.na(current_statistic_type)
    ) {
      next
    }

    d_result <- tryCatch(
      {
        if (
          identical(current_statistic_type, "t") &&
          !is.na(df_value) &&
          is.finite(df_value) &&
          df_value > 0
        ) {

          effectsize::t_to_d(
            t = statistic_value,
            df_error = df_value,
            ci = conf_level
          )

        } else if (
          identical(current_statistic_type, "z") &&
          !is.na(n_observations) &&
          is.finite(n_observations)
        ) {

          effectsize::z_to_d(
            z = statistic_value,
            n = n_observations,
            ci = conf_level
          )

        } else {

          NULL
        }
      },
      error = function(e) NULL
    )

    if (is.null(d_result)) {
      next
    }

    d_result <- as.data.frame(d_result)

    d_column <- first_existing_column(
      d_result,
      c("d", "Cohens_d")
    )

    d_ci_low_column <- first_existing_column(
      d_result,
      c("CI_low", "CI_Low", "Conf.Low")
    )

    d_ci_high_column <- first_existing_column(
      d_result,
      c("CI_high", "CI_High", "Conf.High")
    )

    if (!is.null(d_column)) {
      output$d[i] <- suppressWarnings(
        as.numeric(d_result[[d_column]][1L])
      )
      output$d_method[i] <- current_statistic_type
    }

    if (!is.null(d_ci_low_column)) {
      output$d_CI_low[i] <- suppressWarnings(
        as.numeric(d_result[[d_ci_low_column]][1L])
      )
    }

    if (!is.null(d_ci_high_column)) {
      output$d_CI_high[i] <- suppressWarnings(
        as.numeric(d_result[[d_ci_high_column]][1L])
      )
    }
  }

  # ==========================================================================
  # Final column order
  # ==========================================================================

  output <- output[
    ,
    c(
      "Parameter",
      "Estimate",
      "SE",
      "CI",
      "CI_low",
      "CI_high",
      "Statistic",
      "Statistic_type",
      "df",
      "p",
      "Interpreted_Estimate",
      "Interpreted_CI_low",
      "Interpreted_CI_high",
      "Estimate_type",
      "d",
      "d_CI_low",
      "d_CI_high",
      "d_method"
    ),
    drop = FALSE
  ]

  rownames(output) <- NULL

  output
}
