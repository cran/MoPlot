test_that("moplot_onepred() handles Reverse Helmert contrasts", {

  set.seed(123)

  dat <- data.frame(
    group = factor(rep(c("A", "B", "C", "D"), each = 20))
  )

  dat$y <- c(
    rnorm(20, 10, 2),
    rnorm(20, 12, 2),
    rnorm(20, 15, 2),
    rnorm(20, 18, 2)
  )

  contrasts(dat$group) <- reverse_helmert(4)

  mod <- stats::lm(
    y ~ group,
    data = dat
  )

  plot <- moplot_onepred(
    model = mod,
    predictor = "group",
    caption = FALSE
  )

  expect_s3_class(
    plot,
    "patchwork"
  )
})


test_that("moplot_onepred() displays raw data with show_data = TRUE", {

  set.seed(123)

  dat <- data.frame(
    group = factor(rep(c("A", "B", "C", "D"), each = 20))
  )

  dat$y <- c(
    rnorm(20, 10, 2),
    rnorm(20, 12, 2),
    rnorm(20, 15, 2),
    rnorm(20, 18, 2)
  )

  contrasts(dat$group) <- reverse_helmert(4)

  mod <- stats::lm(
    y ~ group,
    data = dat
  )

  plot <- moplot_onepred(
    model = mod,
    predictor = "group",
    show_data = TRUE,
    caption = FALSE
  )

  expect_s3_class(
    plot,
    "patchwork"
  )
})


test_that("moplot_onepred() publication style assigns distinct colours to contrasts", {

  set.seed(123)

  dat <- data.frame(
    group = factor(
      rep(c("A", "B", "C", "D"), each = 20),
      levels = c("A", "B", "C", "D")
    )
  )

  dat$y <- c(
    rnorm(20, 10, 1),
    rnorm(20, 12, 1),
    rnorm(20, 14, 1),
    rnorm(20, 16, 1)
  )

  contrasts(dat$group) <- reverse_helmert(4)

  mod <- stats::lm(
    y ~ group,
    data = dat
  )

  prep <- .prepare_moplot_onepred(
    model = mod,
    predictor = "group",
    estimate_type = "observed",
    prediction_scale = "response",
    coefficient_scale = "model"
  )

  upper <- .plot_moplot_onepred_upper(
    prep = prep,
    style = "publication"
  )

  build <- ggplot2::ggplot_build(
    upper
  )

  segment_layer_index <- which(
    vapply(
      build$data,
      function(layer) {
        all(c("x", "y", "xend", "yend") %in% names(layer))
      },
      logical(1)
    )
  )

  expect_length(
    segment_layer_index,
    1L
  )

  segment_layer <- build$data[[segment_layer_index]]

  expect_gt(
    length(unique(segment_layer$colour)),
    1L
  )
})

test_that("moplot_onepred() minimal style keeps contrast colours in grayscale", {

  set.seed(123)

  dat <- data.frame(
    group = factor(
      rep(c("A", "B", "C"), each = 20)
    )
  )

  dat$y <- rnorm(
    nrow(dat)
  )

  contrasts(dat$group) <- reverse_helmert(3)

  mod <- stats::lm(
    y ~ group,
    data = dat
  )

  prep <- .prepare_moplot_onepred(
    model = mod,
    predictor = "group",
    estimate_type = "observed",
    prediction_scale = "response",
    coefficient_scale = "model"
  )

  upper <- .plot_moplot_onepred_upper(
    prep = prep,
    style = "minimal"
  )

  build <- ggplot2::ggplot_build(upper)

  segment_layer_index <- which(
    vapply(
      build$data,
      function(layer) {
        all(c("x", "y", "xend", "yend") %in% names(layer))
      },
      logical(1)
    )
  )

  expect_length(
    segment_layer_index,
    1L
  )

  segment_layer <- build$data[[segment_layer_index]]

  expect_equal(
    length(unique(segment_layer$colour)),
    1L
  )
})


test_that("moplot_onepred() handles sliding-difference contrasts", {

  set.seed(123)

  dat <- data.frame(group = factor(rep(c("A", "B", "C", "D"), each = 20)))
  dat$y <- stats::rnorm(80, c(10, 12, 15, 18)[dat$group], 2)

  contrasts(dat$group) <- MASS::contr.sdif(4)

  plot <- suppressMessages(
    moplot_onepred(
      model = stats::lm(y ~ group, data = dat),
      predictor = "group",
      caption = FALSE
    )
  )

  expect_s3_class(plot, "patchwork")
})


test_that("per_difference expresses coefficients as the displayed differences", {

  set.seed(123)

  dat <- data.frame(group = factor(rep(c("A", "B", "C", "D"), each = 20)))
  dat$y <- stats::rnorm(80, c(10, 12, 15, 18)[dat$group], 2)

  contrasts(dat$group) <- stats::contr.helmert(4)

  model <- stats::lm(y ~ group, data = dat)

  prepare <- function(per_difference) {
    suppressMessages(
      .prepare_moplot_onepred(
        model = model,
        predictor = "group",
        estimate_type = "observed",
        per_difference = per_difference
      )
    )
  }

  coefficients <- prepare(FALSE)
  differences <- prepare(TRUE)

  segments <- differences$segment_data
  displayed <- segments$yend - segments$y

  expect_equal(differences$parameter_data$Display_estimate, displayed, tolerance = 1e-8)
  expect_equal(coefficients$parameter_data$Display_estimate, displayed / (2:4), tolerance = 1e-8)

  # Rescaling does not change the inferential quantities.
  expect_equal(differences$parameter_data$d, coefficients$parameter_data$d)

  expect_equal(
    differences$parameter_data$Display_CI_high - differences$parameter_data$Display_CI_low,
    (coefficients$parameter_data$Display_CI_high - coefficients$parameter_data$Display_CI_low) * (2:4)
  )
})


test_that("per_difference rescales ratios on the log scale", {

  set.seed(123)

  dat <- data.frame(group = factor(rep(c("A", "B", "C"), each = 60)))
  dat$success <- stats::rbinom(180, 1, c(0.2, 0.4, 0.7)[dat$group])

  contrasts(dat$group) <- stats::contr.helmert(3)

  model <- stats::glm(success ~ group, family = stats::binomial, data = dat)

  prep <- suppressMessages(
    .prepare_moplot_onepred(
      model = model,
      predictor = "group",
      estimate_type = "predicted",
      coefficient_scale = "interpreted",
      per_difference = TRUE
    )
  )

  expect_equal(
    prep$parameter_data$Display_estimate,
    exp(unname(stats::coef(model)[-1]) / c(1 / 2, 1 / 3)),
    tolerance = 1e-8
  )
})


test_that("moplot_onepred() handles transformed covariates and factor outcomes", {

  set.seed(91)

  dat <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 30)),
    x = stats::rexp(90) + 0.1
  )
  dat$y <- stats::rnorm(90, as.integer(dat$group) + log(dat$x))
  dat$answer <- factor(ifelse(stats::runif(90) < c(0.2, 0.5, 0.7)[dat$group], "yes", "no"))

  plot <- suppressMessages(
    moplot_onepred(stats::lm(y ~ group + log(x), data = dat), "group", caption = FALSE)
  )
  expect_s3_class(plot, "patchwork")

  plot <- suppressMessages(
    moplot_onepred(stats::lm(y ~ group + log(x), data = dat), "group", data = dat, caption = FALSE)
  )
  expect_s3_class(plot, "patchwork")

  binomial_model <- stats::glm(answer ~ group, family = stats::binomial, data = dat)

  prep <- suppressMessages(
    .prepare_moplot_onepred(binomial_model, "group", estimate_type = "observed")
  )

  expect_setequal(unique(prep$raw_data$y), c(0, 1))
})


polynomial_dose_data <- function() {
  set.seed(42)
  dat <- data.frame(
    dose = factor(rep(c("0mg", "10mg", "20mg", "40mg"), each = 30),
                  levels = c("0mg", "10mg", "20mg", "40mg"))
  )
  dat$y <- stats::rnorm(120, c(2.0, 3.4, 4.1, 4.3)[dat$dose], 1)
  dat
}


test_that("moplot_onepred() displays one panel per polynomial component", {

  dat <- polynomial_dose_data()
  contrasts(dat$dose) <- stats::contr.poly(4)
  mod <- stats::lm(y ~ dose, data = dat)

  prep <- suppressMessages(.prepare_moplot_onepred(mod, "dose"))

  expect_identical(prep$position_layout, "polynomial")
  expect_null(prep$segment_data)
  expect_identical(as.numeric(prep$axis_breaks), c(1, 2, 3, 4))

  components <- prep$polynomial_data$components

  expect_identical(components$Component, c("Linear", "Quadratic", "Cubic"))
  expect_identical(components$Label, c("C1*  Linear", "C2*  Quadratic", "C3  Cubic"))

  # Each curve passes through Polynomial_fit at the levels.
  curves <- prep$polynomial_data$curves
  level_fits <- prep$polynomial_data$level_fits

  for (contrast in components$Contrast) {
    curve <- curves[curves$Contrast == contrast, ]
    fits <- level_fits[level_fits$Contrast == contrast, ]
    expect_equal(
      stats::approx(curve$x, curve$y, xout = fits$x)$y,
      fits$Fit,
      tolerance = 1e-3
    )
  }

  # Baseline plus all components reproduces the estimates.
  first <- prep$contrast_data[prep$contrast_data$Contrast == components$Contrast[1L], ]
  expect_equal(prep$polynomial_data$baseline, mean(first$Estimate))

  plot <- suppressMessages(moplot_onepred(mod, "dose", caption = FALSE))

  expect_s3_class(plot, "patchwork")

  upper <- plot[[1L]]
  lower <- plot[[2L]]

  expect_s3_class(upper$facet, "FacetWrap")
  expect_s3_class(lower$facet, "FacetWrap")

  expect_identical(
    ggplot2::ggplot_build(upper)$layout$layout$PANEL,
    factor(1:3)
  )
})


test_that("polynomial plots use the publication palette and grayscale", {

  dat <- polynomial_dose_data()
  contrasts(dat$dose) <- stats::contr.poly(4)
  mod <- stats::lm(y ~ dose, data = dat)

  curve_colours <- function(style) {
    plot <- suppressMessages(moplot_onepred(mod, "dose", style = style, caption = FALSE))
    built <- ggplot2::ggplot_build(plot[[1L]])
    line_layer <- which(vapply(plot[[1L]]$layers, function(layer) inherits(layer$geom, "GeomLine"), logical(1)))
    unique(built$data[[line_layer]]$colour)
  }

  expect_identical(curve_colours("publication"), .moplot_contrast_palette(3))
  expect_identical(curve_colours("minimal"), "#525252")
})


test_that("polynomial plots space the levels by the scores of the coding", {

  dat <- polynomial_dose_data()
  contrasts(dat$dose) <- stats::contr.poly(4, scores = c(0, 10, 20, 40))
  mod <- stats::lm(y ~ dose, data = dat)

  prep <- suppressMessages(.prepare_moplot_onepred(mod, "dose"))

  expect_equal(as.numeric(prep$axis_breaks), c(0, 10, 20, 40) / 40 * 3 + 1)

  # The assignment keeps capture.output() from printing the plot.
  caption <- utils::capture.output(
    plot <- suppressMessages(moplot_onepred(mod, "dose", show_coef = FALSE))
  )

  expect_true(any(grepl("spaced according to the scores", caption)))
})


test_that("polynomial plots of logistic models stay on the probability scale", {

  set.seed(7)
  dat <- data.frame(freq = factor(rep(c("a", "b", "c", "d", "e"), each = 60)))
  dat$success <- stats::rbinom(300, 1, stats::plogis(c(-1.6, -0.5, 0.4, 0.6, 0.3))[dat$freq])
  contrasts(dat$freq) <- stats::contr.poly(5)
  mod <- stats::glm(success ~ freq, family = stats::binomial, data = dat)

  prep <- suppressMessages(
    .prepare_moplot_onepred(mod, "freq", coefficient_scale = "interpreted")
  )

  expect_true(all(prep$polynomial_data$curves$y > 0 & prep$polynomial_data$curves$y < 1))
  expect_identical(prep$polynomial_data$components$Component[4L], "Quartic")
  expect_identical(prep$null_value, 1)

  plot <- suppressMessages(
    moplot_onepred(mod, "freq", coefficient_scale = "interpreted", caption = FALSE)
  )

  expect_s3_class(plot, "patchwork")
})


test_that("per_difference is ignored with a warning for polynomial contrasts", {

  dat <- polynomial_dose_data()
  contrasts(dat$dose) <- stats::contr.poly(4)
  mod <- stats::lm(y ~ dose, data = dat)

  expect_warning(
    plot <- suppressMessages(moplot_onepred(mod, "dose", per_difference = TRUE, caption = FALSE)),
    "`per_difference` is ignored for polynomial contrasts"
  )

  expect_s3_class(plot, "patchwork")
})


test_that("the caption describes polynomial components", {

  dat <- polynomial_dose_data()
  contrasts(dat$dose) <- stats::contr.poly(4)
  mod <- stats::lm(y ~ dose, data = dat)

  caption <- paste(
    utils::capture.output(plot <- suppressMessages(moplot_onepred(mod, "dose"))),
    collapse = " "
  )

  expect_match(caption, "Under polynomial coding")
  expect_match(caption, "assumed to be equally spaced")
  expect_match(caption, "beneath the panel of that component")
  expect_false(grepl("Dashed segments", caption))
})
