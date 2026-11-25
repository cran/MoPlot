test_that("contrasts_multpred() handles Reverse Helmert contrasts correctly", {

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2"))
  )

  dat <- dat[rep(seq_len(nrow(dat)), each = 20), ]

  dat$y <- rnorm(
    nrow(dat),
    mean = with(
      dat,
      ifelse(
        A == "A1", 10,
        ifelse(A == "A2", 12, 15)
      ) +
        ifelse(B == "B2", 2, 0)
    ),
    sd = 2
  )

  contrasts(dat$A) <- reverse_helmert(3)
  contrasts(dat$B) <- stats::contr.treatment(2)

  mod <- stats::lm(
    y ~ A * B,
    data = dat
  )

  out <- contrasts_multpred(
    model = mod,
    data = dat
  )

  expect_true(is.data.frame(out))

  expect_equal(
    sort(unique(out$contrast_id)),
    c(1, 2)
  )

  expect_true(
    all(c(
      "A1",
      "Mean(A2,A3)",
      "A2",
      "A3"
    ) %in% out$A)
  )

  expect_equal(
    unique(out$split_var),
    "B"
  )

  expect_equal(
    sort(unique(out$B)),
    c("B1", "B2")
  )
})

test_that("contrasts_multpred() recovers analysis data from model", {

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2")),
    rep = seq_len(20)
  )

  dat$y <- with(
    dat,
    10 +
      2 * (A == "A2") +
      4 * (A == "A3") +
      2 * (B == "B2") +
      2 * (A == "A2" & B == "B2")
  ) +
    stats::rnorm(
      nrow(dat),
      0,
      1
    )

  contrasts(dat$A) <- reverse_helmert(3)
  contrasts(dat$B) <- MASS::contr.sdif(2)

  mod <- stats::lm(
    y ~ A * B,
    data = dat
  )

  out_explicit <- contrasts_multpred(
    model = mod,
    data = dat,
    configuration = "all",
    methodological_note = FALSE
  )

  out_recovered <- contrasts_multpred(
    model = mod,
    configuration = "all",
    methodological_note = FALSE
  )

  expect_equal(
    out_recovered,
    out_explicit
  )
})


test_that("contrasts_multpred() reports interactions with numeric variables clearly", {

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2")),
    rep = seq_len(20)
  )

  dat$x <- stats::rnorm(nrow(dat))
  dat$y <- with(dat, 10 + 2 * (A == "A2") + x) + stats::rnorm(nrow(dat))

  numeric_only <- stats::lm(y ~ A * x, data = dat)

  expect_error(
    contrasts_multpred(numeric_only, methodological_note = FALSE),
    "No interactions among factors.*numeric variables \\(x\\)"
  )

  mixed <- stats::lm(y ~ A * B * x, data = dat)

  expect_error(
    contrasts_multpred(mixed, interaction = "A:x", methodological_note = FALSE),
    "`A:x` involves numeric variables \\(x\\)"
  )

  expect_error(
    contrasts_multpred(mixed, coefficients = "AA2:x", methodological_note = FALSE),
    "do not belong to interactions among factors: AA2:x"
  )

  # The factor-only interaction is selected by default; its coefficients
  # refer to x = 0, whereas the predictions fix x at its mean.
  expect_warning(
    out <- suppressMessages(
      contrasts_multpred(mixed, methodological_note = FALSE)
    ),
    "do not reconstruct"
  )

  expect_setequal(unique(out$coefficient), c("AA2:BB2", "AA3:BB2"))
})


test_that("contrasts_multpred() restricts supplied data to the rows used for fitting", {

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2")),
    rep = seq_len(20)
  )

  dat$x <- stats::rnorm(nrow(dat))
  dat$y <- with(dat, 10 + 2 * (A == "A2") + 2 * (A == "A2" & B == "B2")) +
    stats::rnorm(nrow(dat))

  mod <- stats::lm(y ~ A * B, data = dat, subset = x > 0)

  out_explicit <- contrasts_multpred(
    model = mod,
    data = dat,
    configuration = "all",
    methodological_note = FALSE
  )

  out_recovered <- contrasts_multpred(
    model = mod,
    configuration = "all",
    methodological_note = FALSE
  )

  expect_equal(out_explicit, out_recovered)
})


# ============================================================
# HYPOTHESIS MATRIX, REFERENCE GRID, AND SCALE FACTORS
# ============================================================

simulate_two_way_data <- function(n_per_cell = 30) {
  set.seed(21)
  data <- expand.grid(
    A = factor(c("a1", "a2", "a3")),
    B = factor(c("b1", "b2")),
    replicate = seq_len(n_per_cell)
  )
  data$replicate <- NULL
  # The covariate distribution differs between cells.
  data$x <- stats::rnorm(nrow(data), mean = 2 * (data$A == "a3" & data$B == "b2"))
  data$y <- stats::rnorm(
    nrow(data),
    mean = (data$A == "a2") + 1.5 * (data$A == "a3") * (data$B == "b2") + 0.8 * data$x
  )
  data$success <- stats::rbinom(nrow(data), 1, stats::plogis(-0.5 + (data$A == "a3") * (data$B == "b2")))
  data
}

# For a two-way coefficient the rows are ordered as positive/positive,
# positive/negative, negative/positive, negative/negative.
displayed_interactions <- function(result) {
  vapply(
    unique(result$coefficient),
    function(coefficient_name) {
      rows <- result[result$coefficient == coefficient_name, ]
      rows$scale_factor[1] *
        ((rows$estimate[1] - rows$estimate[2]) - (rows$estimate[3] - rows$estimate[4]))
    },
    numeric(1)
  )
}

quiet_mult <- function(...) {
  suppressMessages(contrasts_multpred(..., methodological_note = FALSE))
}


test_that("raw interaction contrasts reconstruct the coefficients for every coding", {

  data <- simulate_two_way_data()

  codings <- list(
    treatment_base_3 = stats::contr.treatment(3, base = 3),
    sum = stats::contr.sum(3),
    helmert = stats::contr.helmert(3),
    sliding = MASS::contr.sdif(3),
    custom = cbind(c(-3, 1, 2) / 6, c(0, -1, 1) / 2)
  )

  for (coding_name in names(codings)) {

    stats::contrasts(data$A) <- codings[[coding_name]]
    stats::contrasts(data$B) <- stats::contr.sum(2)

    model <- stats::lm(y ~ A * B, data = data)

    result <- quiet_mult(model = model, data = data, split_rule = "first_variable")

    expect_equal(
      unname(displayed_interactions(result)),
      unname(stats::coef(model)[grep(":", names(stats::coef(model)))]),
      tolerance = 1e-8,
      info = coding_name
    )

    expect_true(all(attr(result, "identity_holds")), info = coding_name)
  }
})


test_that("predicted interaction contrasts are not biased by unbalanced covariates", {

  data <- simulate_two_way_data()
  stats::contrasts(data$A) <- stats::contr.helmert(3)

  model <- stats::lm(y ~ A * B + x, data = data)

  result <- quiet_mult(
    model = model,
    data = data,
    estimate_type = "predicted",
    split_rule = "first_variable"
  )

  expect_equal(
    unname(displayed_interactions(result)),
    unname(stats::coef(model)[grep(":", names(stats::coef(model)))]),
    tolerance = 1e-8
  )

  expect_equal(unname(attr(result, "scale_factors")), c(1 / 2, 1 / 3))
})


test_that("predicted interaction contrasts reconstruct glm coefficients on the link scale", {

  data <- simulate_two_way_data()
  stats::contrasts(data$A) <- stats::contr.helmert(3)

  model <- stats::glm(success ~ A * B + x, family = stats::binomial, data = data)

  result <- quiet_mult(
    model = model,
    data = data,
    estimate_type = "predicted",
    prediction_scale = "link",
    split_rule = "first_variable"
  )

  expect_equal(
    unname(displayed_interactions(result)),
    unname(stats::coef(model)[grep(":", names(stats::coef(model)))]),
    tolerance = 1e-8
  )
})


test_that("lower-order coefficients are evaluated where they are defined", {

  set.seed(22)

  data <- expand.grid(
    A = factor(c("a1", "a2")),
    B = factor(c("b1", "b2")),
    C = factor(c("c1", "c2", "c3")),
    replicate = 1:20
  )
  data$y <- stats::rnorm(nrow(data), (data$A == "a2") * (data$B == "b2") * as.integer(data$C))

  stats::contrasts(data$A) <- stats::contr.sum(2)
  stats::contrasts(data$B) <- stats::contr.sum(2)

  # With treatment-coded C, the A:B coefficient refers to C = c1.
  result <- quiet_mult(
    model = stats::lm(y ~ A * B * C, data = data),
    data = data,
    coefficients = "A1:B1",
    estimate_type = "predicted"
  )

  expect_true(attr(result, "identity_holds")[["A1:B1"]])
})


test_that("non-orthogonal custom codings produce weighted labels and a note", {

  data <- simulate_two_way_data()
  stats::contrasts(data$A) <- cbind(c(-3, 1, 2) / 6, c(0, -1, 1) / 2)

  expect_message(
    result <- contrasts_multpred(
      model = stats::lm(y ~ A * B, data = data),
      data = data,
      configuration = "all",
      methodological_note = FALSE
    ),
    "suggest a3 vs a2"
  )

  expect_true("W(a1 .20, a3 .80)" %in% result$A)
  expect_true(any(startsWith(result$configuration_semantic, "W(a1 .20, a3 .80)_")))
})


test_that("an interaction with a numeric predictor triggers a warning", {

  data <- simulate_two_way_data()

  expect_warning(
    quiet_mult(
      model = stats::lm(y ~ A * B * x, data = data),
      data = data,
      interaction = "A:B",
      estimate_type = "predicted"
    ),
    "do not reconstruct"
  )
})


test_that("moplot_multpred() per_difference divides coefficients by their scale", {

  data <- simulate_two_way_data()
  stats::contrasts(data$A) <- stats::contr.helmert(3)
  stats::contrasts(data$B) <- stats::contr.sum(2)

  model <- stats::lm(y ~ A * B, data = data)

  output <- quiet_mult(model = model, data = data, configuration = "all")
  prep <- .prepare_moplot_multpred_data(output)

  parameters <- .prepare_moplot_multpred_parameters(
    model = model,
    coefficient_data = prep$coefficient_data,
    scale_factors = attr(output, "scale_factors"),
    per_difference = TRUE
  )

  coefficients <- stats::coef(model)[parameters$parameter_data$Parameter]

  expect_equal(
    parameters$parameter_data$Display_estimate,
    unname(coefficients / attr(output, "scale_factors")[parameters$parameter_data$Parameter])
  )

  plot <- suppressMessages(
    moplot_multpred(model = model, data = data, per_difference = TRUE, caption = FALSE, interaction_info = FALSE)
  )

  expect_s3_class(plot, "patchwork")
})


test_that("raw proportions in interaction decompositions stay within [0, 1]", {

  set.seed(23)

  data <- expand.grid(A = factor(c("a1", "a2", "a3")), B = factor(c("b1", "b2")), replicate = 1:25)
  data$success <- stats::rbinom(nrow(data), 1, ifelse(data$A == "a1", 0.01, 0.4))
  stats::contrasts(data$A) <- stats::contr.helmert(3)

  result <- quiet_mult(
    model = stats::glm(success ~ A * B, family = stats::binomial, data = data),
    data = data
  )

  expect_true(all(result$ci_low >= 0 & result$ci_high <= 1))
})


test_that("estimate_type = 'auto' follows the rule of contrasts_onepred()", {

  set.seed(3)
  data <- expand.grid(
    A = factor(c("a1", "a2", "a3")),
    B = factor(c("b1", "b2")),
    subject = factor(1:20)
  )
  data$x <- stats::rnorm(nrow(data))
  data$y <- stats::rnorm(nrow(data), as.integer(data$A) * (data$B == "b2"))
  stats::contrasts(data$A) <- stats::contr.helmert(3)
  stats::contrasts(data$B) <- stats::contr.sum(2)

  resolved <- function(model) {
    attr(
      suppressMessages(contrasts_multpred(model, data = data, methodological_note = FALSE)),
      "estimate_type"
    )
  }

  # The interaction involves all the predictors: observed means.
  expect_identical(resolved(stats::lm(y ~ A * B, data = data)), "observed")

  # Another predictor in the fixed part: model-based predictions.
  expect_identical(resolved(stats::lm(y ~ A * B + x, data = data)), "predicted")

  # Random effects: model-based predictions.
  expect_identical(
    resolved(suppressMessages(lme4::lmer(y ~ A * B + (1 | subject), data = data))),
    "predicted"
  )

  # The former value "raw" is no longer accepted.
  expect_error(
    contrasts_multpred(stats::lm(y ~ A * B, data = data), data = data, estimate_type = "raw"),
    "should be one of"
  )
})
