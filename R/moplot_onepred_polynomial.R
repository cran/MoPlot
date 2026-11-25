# ============================================================
# POLYNOMIAL CONTRASTS IN moplot_onepred()
# ============================================================
#
# Each polynomial component involves every level through its weight, so
# it is displayed in its own panel as the curve it implies across the
# ordered levels, instead of as a segment between two quantities.


#' Plotmath labels of model parameters
#'
#' Internal helper shared by the lower panels of `moplot_onepred()` and
#' `moplot_multpred()`. Labels are plotmath expressions, drawn with
#' `parse = TRUE`, so that beta is rendered on every graphics device;
#' devices such as [grDevices::pdf()] cannot encode the Unicode character.
#'
#' @param display_type The `Display_type` column of the parameter data.
#' @param estimate Displayed estimates.
#' @param d Cohen's d values.
#' @param per_difference Logical. Were the coefficients divided by their
#'   scale factor?
#' @param digits Number of decimal places.
#'
#' @return A list with `estimate_label` (the plotmath name of each
#'   estimate) and `text_label` (the two-line plotmath label).
#'
#' @keywords internal
#' @noRd
.parameter_label_expressions <- function(display_type,
                                         estimate,
                                         d,
                                         per_difference = FALSE,
                                         digits = 2L) {

  estimate_label <- vapply(
    display_type,
    function(estimate_type) {
      switch(
        estimate_type,
        "Model estimate" = "beta",
        "Coefficient" = "beta",
        "Odds ratio" = "\"OR\"",
        "Risk ratio" = "\"RR\"",
        "Incidence rate ratio" = "\"IRR\"",
        "Hazard ratio" = "\"HR\"",
        "Multiplicative ratio" = "\"MR\"",
        "Exponentiated coefficient" = "exp(beta)",
        "\"Estimate\""
      )
    },
    character(1),
    USE.NAMES = FALSE
  )

  # Rescaled model coefficients are labelled explicitly. Ratios keep
  # their usual label: exp(beta / s) is the ratio for the displayed
  # comparison itself.
  if (isTRUE(per_difference)) {
    estimate_label[estimate_label == "beta"] <- "beta/s"
  }

  format_value <- function(values) {
    ifelse(
      is.na(values),
      "NA",
      formatC(values, format = "f", digits = digits)
    )
  }

  # Numbers are quoted so that their decimals are kept as formatted.
  text_label <- paste0(
    "displaystyle(atop(",
    estimate_label,
    " == \"",
    format_value(estimate),
    "\", d == \"",
    format_value(d),
    "\"))"
  )

  list(
    estimate_label = estimate_label,
    text_label = text_label
  )
}


#' Format the lower-panel labels of model parameters
#'
#' Internal helper shared by the lower panels of `moplot_onepred()`.
#'
#' @param parameter_data The `parameter_data` element of a
#'   `"moplot_onepred_data"` object.
#' @param per_difference Logical. Were the coefficients divided by their
#'   scale factor?
#' @param digits Number of decimal places.
#'
#' @return `parameter_data` with the columns `Estimate_label` and
#'   `Text_label` added.
#'
#' @keywords internal
#' @noRd
.moplot_parameter_labels <- function(parameter_data,
                                     per_difference = FALSE,
                                     digits = 2L) {

  labels <- .parameter_label_expressions(
    display_type = parameter_data$Display_type,
    estimate = parameter_data$Display_estimate,
    d = parameter_data$d,
    per_difference = per_difference,
    digits = digits
  )

  parameter_data$Estimate_label <- labels$estimate_label
  parameter_data$Text_label <- labels$text_label

  parameter_data
}


#' Name a polynomial component
#'
#' @param contrast_ids Contrast identifiers returned by `contrasts_onepred()`,
#'   such as `"Contrast_2_.Q"`.
#'
#' @return A character vector: `"Linear"`, `"Quadratic"`, `"Cubic"`,
#'   `"Quartic"`, `"Quintic"`, then `"Degree 6"` and so forth.
#'
#' @keywords internal
#' @noRd
.polynomial_component_names <- function(contrast_ids) {

  degrees <- as.integer(sub("^Contrast_([0-9]+)_.*$", "\\1", contrast_ids))

  named_degrees <- c("Linear", "Quadratic", "Cubic", "Quartic", "Quintic")

  ifelse(
    degrees <= length(named_degrees),
    named_degrees[pmin(degrees, length(named_degrees))],
    paste("Degree", degrees)
  )
}


#' Horizontal positions of the levels under polynomial coding
#'
#' The linear column of a polynomial coding is proportional to the centered
#' scores of the levels. Positions keep the relative spacing of the scores
#' and span 1 to k, as the equally spaced positions of the other layouts.
#'
#' @param linear_weights Weights of the linear component, in level order.
#'
#' @return A numeric vector of positions, in level order.
#'
#' @keywords internal
#' @noRd
.polynomial_level_positions <- function(linear_weights) {

  number_levels <- length(linear_weights)

  weight_range <- range(linear_weights)

  positions <- 1 + (number_levels - 1) *
    (linear_weights - weight_range[1L]) / diff(weight_range)

  # Equal spacing is exact rather than up to rounding.
  if (isTRUE(all.equal(positions, seq_len(number_levels)))) {
    positions <- as.numeric(seq_len(number_levels))
  }

  positions
}


#' Curves and labels of the polynomial components
#'
#' Each component curve is the baseline plus the component coefficient
#' multiplied by the orthogonal polynomial, evaluated on a fine grid between
#' the levels and transformed to the scale of the displayed estimates. At
#' the levels it equals `Polynomial_fit` returned by `contrasts_onepred()`.
#'
#' @param contrast_data Output of `contrasts_onepred()` with a column `x`.
#' @param coefficients Model coefficients of the components, in order.
#' @param significant Logical vector, one value per component.
#' @param linkfun,linkinv Functions mapping the displayed scale to the
#'   linear-predictor scale and back.
#' @param grid_size Number of grid points of each curve.
#'
#' @return A list with elements `components` (one row per component with
#'   `Contrast`, `Contrast_number`, `Component`, `Significant`, `Label`),
#'   `level_fits` (one row per component and level with `x` and `Fit`),
#'   `curves` (`Contrast`, `x`, `y`), and `baseline` (on the displayed
#'   scale).
#'
#' @keywords internal
#' @noRd
.polynomial_component_data <- function(contrast_data,
                                       coefficients,
                                       significant,
                                       linkfun = identity,
                                       linkinv = identity,
                                       grid_size = 200L) {

  contrast_ids <- unique(contrast_data$Contrast)

  number_components <- length(contrast_ids)

  component_names <- .polynomial_component_names(contrast_ids)

  components <- data.frame(
    Contrast = contrast_ids,
    Contrast_number = seq_len(number_components),
    Component = component_names,
    Significant = significant,
    Label = paste0(
      "C",
      seq_len(number_components),
      ifelse(significant, "*", ""),
      "  ",
      component_names
    ),
    stringsAsFactors = FALSE
  )

  level_fits <- data.frame(
    Contrast = contrast_data$Contrast,
    x = contrast_data$x,
    Fit = contrast_data$Polynomial_fit,
    stringsAsFactors = FALSE
  )

  # The baseline is the same for every component: the fitted value minus
  # the contribution of the component, on the linear-predictor scale.
  baseline_values <- linkfun(contrast_data$Polynomial_fit) -
    coefficients[match(contrast_data$Contrast, contrast_ids)] *
    contrast_data$Weight

  baseline <- if (all(is.finite(baseline_values))) {
    mean(baseline_values)
  } else {
    NA_real_
  }

  if (!is.finite(baseline)) {

    return(
      list(
        components = components,
        level_fits = level_fits,
        curves = data.frame(Contrast = character(), x = numeric(), y = numeric()),
        baseline = NA_real_
      )
    )
  }

  first_component <- contrast_data[contrast_data$Contrast == contrast_ids[1L], , drop = FALSE]

  level_positions <- first_component$x

  grid <- seq(min(level_positions), max(level_positions), length.out = grid_size)

  # Orthonormal polynomials on the level positions, which reproduce the
  # coding columns up to their sign, evaluated between the levels.
  basis <- stats::poly(level_positions, degree = number_components)

  grid_basis <- stats::predict(basis, grid)

  curves <- do.call(
    rbind,
    lapply(
      seq_len(number_components),
      function(component_index) {

        weights <- contrast_data$Weight[contrast_data$Contrast == contrast_ids[component_index]]

        alignment <- sum(basis[, component_index] * weights)

        component_grid <- if (abs(abs(alignment) - 1) < 1e-6) {
          sign(alignment) * grid_basis[, component_index]
        } else {
          # Not an orthonormal polynomial in the positions: join the
          # weights of the levels instead.
          stats::approx(level_positions, weights, xout = grid)$y
        }

        data.frame(
          Contrast = contrast_ids[component_index],
          x = grid,
          y = linkinv(baseline + coefficients[component_index] * component_grid),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  list(
    components = components,
    level_fits = level_fits,
    curves = curves,
    baseline = linkinv(baseline)
  )
}


#' Graphical settings of the polynomial panels
#'
#' @param style `"minimal"` or `"publication"`.
#' @param number_components Number of polynomial components.
#'
#' @return A list of colours and sizes.
#'
#' @keywords internal
#' @noRd
.polynomial_style <- function(style,
                              number_components) {

  if (identical(style, "minimal")) {

    list(
      component_colours = rep("#525252", number_components),
      marker_colours = rep("#252525", number_components),
      estimate_point_size = 2.8,
      estimate_linewidth = 0.60,
      curve_linewidth = 0.90,
      marker_size = 3.6,
      marker_stroke = 1.0,
      raw_point_size = 1.5,
      raw_alpha = 0.25,
      coefficient_size = 4.5,
      coefficient_stroke = 1.0,
      interval_linewidth = 0.70,
      guide_linewidth = 0.40,
      label_size = 3.1,
      label_border_width = 0.20,
      face = "plain"
    )

  } else {

    colours <- .moplot_contrast_palette(number_components)

    list(
      component_colours = colours,
      marker_colours = colours,
      estimate_point_size = 3.5,
      estimate_linewidth = 0.80,
      curve_linewidth = 1.10,
      marker_size = 4.2,
      marker_stroke = 1.2,
      raw_point_size = 1.8,
      raw_alpha = 0.30,
      coefficient_size = 5.5,
      coefficient_stroke = 1.2,
      interval_linewidth = 0.85,
      guide_linewidth = 0.45,
      label_size = 3.3,
      label_border_width = 0.22,
      face = "bold"
    )
  }
}


#' Shared theme of the polynomial panels
#'
#' @keywords internal
#' @noRd
.polynomial_theme <- function(base_size,
                              face) {

  ggplot2::theme_light(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92"),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = face),
      axis.text = ggplot2::element_text(colour = "#252525"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(colour = "#252525", face = face, size = ggplot2::rel(1))
    )
}


#' Build the upper panel for polynomial contrasts
#'
#' Internal helper used by `moplot_onepred()`. One facet per component
#' shows the displayed estimates with their confidence intervals, the curve
#' implied by the component alone around the baseline, and the value of
#' the curve at each level.
#'
#' @param prep An object returned by `.prepare_moplot_onepred()` for a
#'   predictor with polynomial contrasts.
#' @param show_data Logical. Should raw observations be displayed?
#' @param style `"minimal"` or `"publication"`.
#' @param base_size Base font size.
#' @param label_key Abbreviation key returned by
#'   `.abbreviate_moplot_labels()`.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd
.plot_moplot_onepred_polynomial_upper <- function(prep,
                                                  show_data = FALSE,
                                                  style = c("minimal", "publication"),
                                                  base_size = 13,
                                                  label_key = NULL) {

  style <- match.arg(style)

  polynomial_data <- prep$polynomial_data

  if (is.null(polynomial_data)) {
    stop(
      "`prep` does not contain polynomial components.",
      call. = FALSE
    )
  }

  components <- polynomial_data$components

  settings <- .polynomial_style(style, nrow(components))

  facet_levels <- components$Contrast

  facet_labels <- stats::setNames(components$Label, components$Contrast)

  as_facet <- function(contrasts) factor(contrasts, levels = facet_levels)

  # Every facet repeats the level estimates and, optionally, the raw data.
  repeat_by_facet <- function(rows) {
    if (nrow(rows) == 0L) {
      return(cbind(rows, Contrast = factor(character(), levels = facet_levels)))
    }
    do.call(
      rbind,
      lapply(facet_levels, function(contrast) cbind(rows, Contrast = as_facet(contrast)))
    )
  }

  estimate_data <- prep$contrast_data

  estimate_data$Contrast <- as_facet(estimate_data$Contrast)

  level_fits <- polynomial_data$level_fits

  marker_colours <- stats::setNames(settings$marker_colours, facet_levels)

  level_fits$Colour <- unname(marker_colours[level_fits$Contrast])

  level_fits$Contrast <- as_facet(level_fits$Contrast)

  component_colours <- stats::setNames(settings$component_colours, facet_levels)

  curves <- polynomial_data$curves

  curves$Colour <- unname(component_colours[curves$Contrast])

  curves$Contrast <- as_facet(curves$Contrast)

  axis_labels <- .apply_moplot_label_key(
    labels = as.character(prep$axis_labels),
    label_key = label_key
  )

  raw_data <- prep$raw_data[
    is.finite(prep$raw_data$x) & is.finite(prep$raw_data$y),
    c("x", "y"),
    drop = FALSE
  ]

  plot <- ggplot2::ggplot()

  if (isTRUE(show_data) && nrow(raw_data) > 0L) {

    plot <- plot +
      ggplot2::geom_jitter(
        data = repeat_by_facet(raw_data),
        mapping = ggplot2::aes(x = x, y = y),
        width = 0.06,
        height = 0,
        alpha = settings$raw_alpha,
        size = settings$raw_point_size,
        colour = "#D9D9D9",
        inherit.aes = FALSE
      )
  }

  if (is.finite(polynomial_data$baseline)) {

    plot <- plot +
      ggplot2::geom_hline(
        yintercept = polynomial_data$baseline,
        linewidth = settings$guide_linewidth,
        linetype = "dotted",
        colour = "#BDBDBD"
      )
  }

  plot <- plot +
    ggplot2::geom_line(
      data = curves,
      mapping = ggplot2::aes(x = x, y = y, colour = Colour),
      linewidth = settings$curve_linewidth,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_errorbar(
      data = estimate_data,
      mapping = ggplot2::aes(x = x, ymin = CI_low, ymax = CI_high),
      width = 0.08,
      linewidth = settings$estimate_linewidth,
      colour = "#737373",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = estimate_data,
      mapping = ggplot2::aes(x = x, y = Estimate),
      size = settings$estimate_point_size,
      colour = "#737373",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = level_fits[is.finite(level_fits$Fit), , drop = FALSE],
      mapping = ggplot2::aes(x = x, y = Fit, colour = Colour),
      shape = 4,
      size = settings$marker_size,
      stroke = settings$marker_stroke,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::facet_wrap(
      ggplot2::vars(Contrast),
      nrow = 1L,
      labeller = ggplot2::as_labeller(facet_labels)
    ) +
    ggplot2::scale_x_continuous(
      breaks = as.numeric(prep$axis_breaks),
      labels = axis_labels,
      expand = ggplot2::expansion(add = 0.4)
    ) +
    ggplot2::labs(
      x = prep$predictor,
      y = prep$outcome
    ) +
    .polynomial_theme(base_size, settings$face) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = ggplot2::margin(t = 8, r = 12, b = 10, l = 10)
    )

  plot
}


#' Build the lower panel for polynomial contrasts
#'
#' Internal helper used by `moplot_onepred()`. Each coefficient is shown
#' beneath the facet of its component, with its confidence interval and an
#' approximate Cohen's d.
#'
#' @inheritParams .plot_moplot_onepred_polynomial_upper
#' @param digits Number of decimal places used in the labels.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd
.plot_moplot_onepred_polynomial_lower <- function(prep,
                                                  style = c("minimal", "publication"),
                                                  digits = 2L,
                                                  base_size = 13) {

  style <- match.arg(style)

  components <- prep$polynomial_data$components

  settings <- .polynomial_style(style, nrow(components))

  parameter_data <- .moplot_parameter_labels(
    parameter_data = prep$parameter_data,
    per_difference = FALSE,
    digits = digits
  )

  parameter_data$Contrast <- factor(parameter_data$Contrast, levels = components$Contrast)

  marker_colours <- stats::setNames(settings$marker_colours, components$Contrast)

  displayed_y <- c(
    parameter_data$Display_estimate,
    parameter_data$Display_CI_low,
    parameter_data$Display_CI_high,
    prep$null_value
  )

  displayed_y <- displayed_y[is.finite(displayed_y)]

  y_span <- diff(range(displayed_y))

  if (!is.finite(y_span) || y_span == 0) {
    y_span <- max(abs(displayed_y), 1)
  }

  y_limits <- range(displayed_y) + c(-0.12, 0.12) * y_span

  y_title <- if (identical(prep$coefficient_scale, "model")) {
    "Model estimate"
  } else {
    unique(parameter_data$Display_type)[1L]
  }

  ggplot2::ggplot(parameter_data) +
    ggplot2::geom_hline(
      yintercept = prep$null_value,
      linewidth = settings$guide_linewidth,
      linetype = "dotted",
      colour = "#BDBDBD"
    ) +
    ggplot2::geom_linerange(
      mapping = ggplot2::aes(x = 0, ymin = Display_CI_low, ymax = Display_CI_high),
      linewidth = settings$interval_linewidth,
      colour = "#969696"
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = 0, y = Display_estimate, colour = Contrast),
      shape = 4,
      size = settings$coefficient_size,
      stroke = settings$coefficient_stroke,
      show.legend = FALSE
    ) +
    ggplot2::geom_label(
      mapping = ggplot2::aes(x = 0.12, y = Display_estimate, label = Text_label),
      parse = TRUE,
      hjust = 0,
      vjust = 0.5,
      size = settings$label_size,
      fill = "white",
      colour = "#252525",
      linewidth = settings$label_border_width,
      label.r = grid::unit(0.12, "lines"),
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = marker_colours, guide = "none") +
    ggplot2::facet_wrap(ggplot2::vars(Contrast), nrow = 1L) +
    ggplot2::scale_x_continuous(limits = c(-0.5, 0.9), breaks = NULL) +
    ggplot2::coord_cartesian(ylim = y_limits, clip = "off") +
    ggplot2::labs(x = NULL, y = y_title) +
    .polynomial_theme(base_size, settings$face) +
    ggplot2::theme(
      strip.text = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 4, r = 12, b = 10, l = 10)
    )
}
