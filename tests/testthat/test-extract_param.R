test_that("linear mixed models use Satterthwaite degrees of freedom", {

  set.seed(111)

  # A between-subject predictor: 24 subjects with 10 observations each.
  data <- data.frame(id = factor(rep(1:24, each = 10)))
  data$group <- factor(rep(rep(c("A", "B", "C"), each = 8), each = 10))
  data$y <- stats::rnorm(240, c(0, 0.5, 1)[data$group]) + stats::rnorm(24)[data$id]

  model <- lme4::lmer(y ~ group + (1 | id), data = data)

  result <- extract_param(model)
  satterthwaite <- summary(lmerTest::as_lmerModLmerTest(model))$coefficients[-1, ]

  expect_equal(result$df, unname(satterthwaite[, "df"]), tolerance = 1e-6)
  expect_equal(result$p, .extract_p_values(model, result$Parameter)$table$p_value, tolerance = 1e-6)
})


test_that("confidence intervals of generalized models agree with their p-values", {

  set.seed(112)

  agreement <- vapply(
    seq_len(60),
    function(replicate) {
      data <- data.frame(x = factor(rep(c("A", "B"), each = 30)))
      data$success <- stats::rbinom(60, 1, c(0.15, 0.35)[data$x])
      result <- extract_param(stats::glm(success ~ x, family = stats::binomial, data = data))
      excludes_zero <- result$CI_low > 0 | result$CI_high < 0
      identical(excludes_zero, result$p < 0.05)
    },
    logical(1)
  )

  expect_true(all(agreement))
})


test_that("interpreted estimates are exponentiated coefficients and limits", {

  set.seed(113)

  data <- data.frame(group = factor(rep(c("A", "B", "C"), each = 50)))
  data$count <- stats::rpois(150, c(1, 2, 3)[data$group])

  result <- extract_param(stats::glm(count ~ group, family = stats::poisson, data = data))

  expect_identical(unique(result$Estimate_type), "Incidence rate ratio")
  expect_equal(result$Interpreted_Estimate, exp(result$Estimate))
  expect_equal(result$Interpreted_CI_low, exp(result$CI_low))
  expect_equal(result$Interpreted_CI_high, exp(result$CI_high))
})


test_that("binomial models with a log link report risk ratios", {

  set.seed(114)

  data <- data.frame(group = factor(rep(c("A", "B"), each = 100)))
  data$success <- stats::rbinom(200, 1, c(0.2, 0.3)[data$group])

  model <- stats::glm(success ~ group, family = stats::binomial(link = "log"), data = data)

  expect_identical(unique(extract_param(model)$Estimate_type), "Risk ratio")
})


test_that("logistic models convert the log odds ratio to d (Chinn, 2000)", {

  set.seed(115)

  data <- data.frame(group = factor(rep(c("A", "B", "C"), each = 60)))
  data$success <- stats::rbinom(180, 1, stats::plogis(-0.5 + c(0, 0.4, 0.8)[data$group]))
  stats::contrasts(data$group) <- stats::contr.helmert(3)

  model <- stats::glm(success ~ group, family = stats::binomial, data = data)

  result <- extract_param(model)

  expect_equal(result$d, result$Estimate * sqrt(3) / pi)
  expect_equal(result$d_CI_low, result$CI_low * sqrt(3) / pi)
  expect_true(all(result$d_method == "log odds ratio"))

  # The plot reports d for the displayed difference: d / s.
  prep <- suppressMessages(.prepare_moplot_onepred(model, "group", estimate_type = "predicted"))

  expect_equal(prep$parameter_data$d, result$d / c(1 / 2, 1 / 3))
})


test_that("other models convert t or z statistics to d", {

  set.seed(116)

  data <- data.frame(group = factor(rep(c("A", "B"), each = 40)))
  data$y <- stats::rnorm(80, c(0, 0.5)[data$group])
  data$count <- stats::rpois(80, c(2, 3)[data$group])

  linear <- extract_param(stats::lm(y ~ group, data = data))
  poisson <- extract_param(stats::glm(count ~ group, family = stats::poisson, data = data))

  expect_identical(linear$d_method, "t")
  expect_equal(linear$d, 2 * linear$Statistic / sqrt(linear$df))
  expect_identical(poisson$d_method, "z")
})
