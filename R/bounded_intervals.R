# ============================================================
# CONFIDENCE INTERVALS FOR OBSERVED MEANS WITHIN THE NATURAL RANGE
# ============================================================
#
# Symmetric intervals for observed means can leave the natural range
# of the response: a proportion close to 0 may receive a negative
# lower bound. Observed (raw) estimates therefore use intervals that
# respect the range of the response:
#
#   - proportions (binary outcomes or binomial models): Wilson score
#     interval; for averages of several cells, the Wilson interval is
#     computed with the effective sample size p (1 - p) / SE^2
#     (Korn and Graubard), which reduces to the classical Wilson
#     interval for a single cell;
#   - positive responses (count and positive-continuous families):
#     interval on the log scale, exp(log(m) +/- z SE / m); a mean of
#     zero receives the exact Poisson upper bound -log(alpha / 2) / n;
#   - other responses: symmetric interval, mean +/- critical * SE.
# ============================================================


#' Validate user-supplied bounds of the response
#'
#' Internal helper. `bounds` declares the natural range of a response
#' analyzed with a model that does not know it, such as a Likert score
#' modelled with a Gaussian model. Generalized models already respect the
#' range of their family, so `bounds` is ignored for them with a warning.
#'
#' @param bounds `NULL` or a numeric vector `c(lower, upper)`; either value
#'   may be infinite, but not both.
#' @param model The fitted model.
#'
#' @return The validated bounds, or `NULL`.
#'
#' @keywords internal
#' @noRd
.validate_bounds <- function(bounds,
                             model) {

  if (is.null(bounds)) {
    return(NULL)
  }

  if (!is.numeric(bounds) ||
      length(bounds) != 2L ||
      anyNA(bounds) ||
      bounds[1L] >= bounds[2L] ||
      all(is.infinite(bounds))) {

    stop(
      paste0(
        "`bounds` must be NULL or a numeric vector c(lower, upper) with ",
        "lower < upper; one of the two values may be infinite, ",
        "for example c(1, 7) or c(0, Inf)."
      ),
      call. = FALSE
    )
  }

  if (.is_generalized_model(model)) {

    warning(
      paste0(
        "`bounds` is ignored for generalized models: their intervals already ",
        "respect the range of the response implied by the model family."
      ),
      call. = FALSE
    )

    return(NULL)
  }

  as.numeric(bounds)
}


#' Confidence intervals within user-supplied bounds
#'
#' Internal helper. The interval is computed on a scale on which the bounds
#' are infinitely far and transformed back:
#'
#'   - two finite bounds [a, b]: logit of (m - a) / (b - a);
#'   - lower bound only, [a, Inf): log of (m - a);
#'   - upper bound only, (-Inf, b]: log of (b - m).
#'
#' The standard error on the transformed scale is obtained with the delta
#' method. Estimates on or outside a bound, or with a null standard error,
#' keep the symmetric interval restricted to the bounds.
#'
#' @param estimate Estimates.
#' @param se Standard errors on the scale of the estimates.
#' @param bounds Validated bounds `c(lower, upper)`.
#' @param critical Critical value (normal or t quantile), recycled.
#'
#' @return A list with `ci_low` and `ci_high`.
#'
#' @keywords internal
#' @noRd
.bounded_interval <- function(estimate,
                              se,
                              bounds,
                              critical) {

  lower <- bounds[1L]
  upper <- bounds[2L]

  critical <- rep_len(critical, length(estimate))

  ci_low <- estimate - critical * se
  ci_high <- estimate + critical * se

  inside <- is.finite(estimate) &
    is.finite(se) &
    se > 0 &
    estimate > lower &
    estimate < upper

  if (is.finite(lower) && is.finite(upper)) {

    width <- upper - lower
    position <- (estimate[inside] - lower) / width
    logit_se <- se[inside] / (width * position * (1 - position))

    ci_low[inside] <- lower + width * stats::plogis(stats::qlogis(position) - critical[inside] * logit_se)
    ci_high[inside] <- lower + width * stats::plogis(stats::qlogis(position) + critical[inside] * logit_se)

  } else if (is.finite(lower)) {

    distance <- estimate[inside] - lower
    log_se <- se[inside] / distance

    ci_low[inside] <- lower + exp(log(distance) - critical[inside] * log_se)
    ci_high[inside] <- lower + exp(log(distance) + critical[inside] * log_se)

  } else {

    distance <- upper - estimate[inside]
    log_se <- se[inside] / distance

    ci_low[inside] <- upper - exp(log(distance) + critical[inside] * log_se)
    ci_high[inside] <- upper - exp(log(distance) - critical[inside] * log_se)
  }

  list(
    ci_low = pmin(pmax(ci_low, lower), upper),
    ci_high = pmin(pmax(ci_high, lower), upper)
  )
}


#' Check that observed values lie within the declared bounds
#'
#' @keywords internal
#' @noRd
.check_values_within_bounds <- function(values,
                                        bounds) {

  if (is.null(bounds)) {
    return(invisible(NULL))
  }

  values <- values[is.finite(values)]

  if (any(values < bounds[1L] | values > bounds[2L])) {

    stop(
      paste0(
        "Some values of the response lie outside `bounds` = c(",
        bounds[1L], ", ", bounds[2L], "). Check the data or the bounds."
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' Warn when estimates lie outside the declared bounds
#'
#' @keywords internal
#' @noRd
.warn_estimates_outside_bounds <- function(estimates,
                                           bounds) {

  if (is.null(bounds)) {
    return(invisible(NULL))
  }

  estimates <- estimates[is.finite(estimates)]

  if (any(estimates <= bounds[1L] | estimates >= bounds[2L])) {

    warning(
      paste0(
        "Some estimates lie on or outside `bounds`: the model implies values ",
        "at or beyond the limits of the response scale. Their intervals were ",
        "restricted to the bounds; consider a model suited to the scale of ",
        "the response."
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' Natural range of the response for observed intervals
#'
#' @param family_object Model family, or `NULL`.
#' @param values Observed response values.
#'
#' @return One of `"proportion"`, `"positive"`, `"unbounded"`.
#'
#' @keywords internal
#' @noRd
.response_range_type <- function(family_object,
                                  values) {

  values <- values[!is.na(values)]

  if (length(values) == 0L || !is.numeric(values)) {
    return("unbounded")
  }

  family_name <- if (is.null(family_object)) "" else tolower(family_object$family)

  is_binomial_family <- family_name %in% c("binomial", "quasibinomial")

  if ((is_binomial_family || all(values %in% c(0, 1))) &&
      all(values >= 0 & values <= 1)) {

    return("proportion")
  }

  is_positive_family <- grepl(
    "poisson|negative binomial|nbinom|gamma|inverse.gaussian|tweedie|genpois|compois",
    family_name
  )

  if (is_positive_family && all(values >= 0)) {
    return("positive")
  }

  "unbounded"
}


#' Confidence intervals for observed means within the natural range
#'
#' @param estimate Observed means (or weighted averages of means).
#' @param se Their standard errors.
#' @param n Number of observations; for a single cell, the cell size. For an
#'   average of cells, the total number of observations, used only when the
#'   effective sample size cannot be computed (estimate at a bound).
#' @param range_type Output of `.response_range_type()`.
#' @param conf_level Confidence level.
#' @param critical Critical value for the symmetric interval of unbounded
#'   responses (for example a t quantile). Defaults to the normal quantile.
#' @param single_cell Logical. For proportions, should `n` be used directly
#'   (single cells) rather than the effective sample size (averages)?
#'
#' @return A list with `ci_low` and `ci_high`.
#'
#' @keywords internal
#' @noRd
.observed_interval <- function(estimate,
                               se,
                               n,
                               range_type,
                               conf_level = 0.95,
                               critical = NULL,
                               single_cell = TRUE,
                               bounds = NULL) {

  z_value <- stats::qnorm(1 - (1 - conf_level) / 2)

  if (is.null(critical)) {
    critical <- z_value
  }

  # User-supplied bounds take precedence over the automatic range.
  if (!is.null(bounds)) {
    return(.bounded_interval(estimate, se, bounds, critical))
  }

  n <- rep_len(n, length(estimate))
  critical <- rep_len(critical, length(estimate))

  ci_low <- estimate - critical * se
  ci_high <- estimate + critical * se

  if (identical(range_type, "proportion")) {

    # Effective sample size of an average of proportions.
    sample_size <- n

    if (!single_cell) {

      effective <- estimate * (1 - estimate) / se^2

      usable <- is.finite(effective) & estimate > 0 & estimate < 1 & se > 0

      sample_size[usable] <- effective[usable]
    }

    valid <- is.finite(estimate) & is.finite(sample_size) & sample_size > 0

    denominator <- 1 + z_value^2 / sample_size
    centre <- (estimate + z_value^2 / (2 * sample_size)) / denominator
    half_width <- z_value / denominator *
      sqrt(estimate * (1 - estimate) / sample_size + z_value^2 / (4 * sample_size^2))

    ci_low[valid] <- pmax(0, centre[valid] - half_width[valid])
    ci_high[valid] <- pmin(1, centre[valid] + half_width[valid])

    # Exact bounds at the edges, avoiding floating-point residues.
    ci_low[valid & estimate == 0] <- 0
    ci_high[valid & estimate == 1] <- 1
  }

  if (identical(range_type, "positive")) {

    on_log_scale <- is.finite(estimate) & estimate > 0 & is.finite(se)

    ci_low[on_log_scale] <- exp(log(estimate[on_log_scale]) - z_value * se[on_log_scale] / estimate[on_log_scale])
    ci_high[on_log_scale] <- exp(log(estimate[on_log_scale]) + z_value * se[on_log_scale] / estimate[on_log_scale])

    # All observations equal to zero: exact Poisson upper bound.
    at_zero <- is.finite(estimate) & estimate == 0 & is.finite(n) & n > 0

    ci_low[at_zero] <- 0
    ci_high[at_zero] <- -log((1 - conf_level) / 2) / n[at_zero]
  }

  list(
    ci_low = ci_low,
    ci_high = ci_high
  )
}
