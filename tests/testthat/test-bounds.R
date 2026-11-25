simulate_likert_data <- function() {
  set.seed(101)
  data <- data.frame(group = factor(rep(c("A", "B", "C"), each = 8)))
  data$score <- c(sample(3:5, 8, TRUE), sample(4:6, 8, TRUE), c(7, 7, 7, 7, 7, 6, 7, 5))
  data
}


test_that(".bounded_interval() keeps intervals within the bounds", {

  interval <- .bounded_interval(estimate = 6.6, se = 0.26, bounds = c(1, 7), critical = 1.96)

  expect_lt(interval$ci_high, 7)
  expect_gt(interval$ci_low, 1)

  # Asymmetric near the upper bound.
  expect_gt(6.6 - interval$ci_low, interval$ci_high - 6.6)
})


test_that(".bounded_interval() is close to the symmetric interval at mid-scale", {

  interval <- .bounded_interval(estimate = 4, se = 0.1, bounds = c(1, 7), critical = 1.96)

  expect_equal(interval$ci_low, 4 - 1.96 * 0.1, tolerance = 1e-3)
  expect_equal(interval$ci_high, 4 + 1.96 * 0.1, tolerance = 1e-3)
})


test_that(".bounded_interval() handles one-sided bounds", {

  lower_only <- .bounded_interval(estimate = 0.15, se = 0.1, bounds = c(0, Inf), critical = 1.96)

  expect_gt(lower_only$ci_low, 0)
  expect_equal(lower_only$ci_low * lower_only$ci_high, 0.15^2, tolerance = 1e-8)

  upper_only <- .bounded_interval(estimate = 99, se = 1, bounds = c(-Inf, 100), critical = 1.96)

  expect_lt(upper_only$ci_high, 100)
})


test_that("bounds apply to observed and predicted estimates", {

  data <- simulate_likert_data()
  model <- stats::lm(score ~ group, data = data)

  unbounded <- suppressMessages(
    contrasts_onepred(model, "group", data, "score", estimate_type = "observed")
  )

  expect_gt(max(unbounded$CI_high), 7)

  for (estimate_type in c("observed", "predicted")) {

    result <- suppressMessages(
      contrasts_onepred(model, "group", data, "score", estimate_type = estimate_type, bounds = c(1, 7))
    )

    expect_true(all(result$CI_low >= 1 & result$CI_high <= 7), info = estimate_type)
    expect_equal(result$Estimate, unbounded$Estimate, tolerance = 1e-8, info = estimate_type)
    expect_identical(attr(result, "bounds"), c(1, 7))
  }
})


test_that("bounds are validated and checked against the data", {

  data <- simulate_likert_data()
  model <- stats::lm(score ~ group, data = data)

  expect_error(
    contrasts_onepred(model, "group", data, "score", bounds = c(7, 1)),
    "lower < upper"
  )

  expect_error(
    contrasts_onepred(model, "group", data, "score", bounds = c(-Inf, Inf)),
    "lower < upper"
  )

  expect_error(
    suppressMessages(contrasts_onepred(model, "group", data, "score", estimate_type = "observed", bounds = c(1, 6))),
    "outside `bounds`"
  )

  data$high <- data$score > 5

  expect_warning(
    suppressMessages(
      contrasts_onepred(stats::glm(high ~ group, family = stats::binomial, data = data), "group", data, "high", bounds = c(0, 1))
    ),
    "ignored for generalized models"
  )
})


test_that("estimates outside the bounds trigger a warning", {

  expect_warning(.warn_estimates_outside_bounds(c(3, 7.2), c(1, 7)), "outside `bounds`")
  expect_silent(.warn_estimates_outside_bounds(c(3, 6.9), c(1, 7)))
})


test_that("contrasts_multpred() and the plotting functions accept bounds", {

  set.seed(102)

  data <- expand.grid(A = factor(c("a1", "a2")), B = factor(c("b1", "b2")), replicate = 1:8)
  data$score <- pmin(7, pmax(1, round(stats::rnorm(32, 5 + 1.8 * (data$A == "a2") * (data$B == "b2")))))

  model <- stats::lm(score ~ A * B, data = data)

  for (estimate_type in c("observed", "predicted")) {

    result <- suppressMessages(
      contrasts_multpred(model, data = data, estimate_type = estimate_type, bounds = c(1, 7), methodological_note = FALSE)
    )

    expect_true(all(result$ci_low >= 1 & result$ci_high <= 7), info = estimate_type)
  }

  expect_s3_class(
    suppressMessages(moplot_multpred(model, data = data, bounds = c(1, 7), caption = FALSE, interaction_info = FALSE)),
    "patchwork"
  )

  likert <- simulate_likert_data()

  expect_s3_class(
    suppressMessages(moplot_onepred(stats::lm(score ~ group, data = likert), "group", bounds = c(1, 7), caption = FALSE)),
    "patchwork"
  )
})
