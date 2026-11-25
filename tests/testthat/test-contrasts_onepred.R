simulate_one_predictor_data <- function(k = 4, n_per_level = 30) {
  set.seed(11)
  data <- data.frame(
    group = factor(rep(LETTERS[seq_len(k)], each = n_per_level)),
    group_size = stats::rnorm(k * n_per_level)
  )
  data$y <- stats::rnorm(nrow(data), mean = c(0, 0.3, 1.5, 2)[seq_len(k)][data$group])
  data$success <- stats::rbinom(nrow(data), 1, c(0.05, 0.03, 0.5, 0.6)[seq_len(k)][data$group])
  data
}

quiet_contrasts <- function(...) {
  suppressMessages(contrasts_onepred(...))
}

# Each coefficient must equal Scale_factor times the displayed difference.
displayed_coefficients <- function(result) {
  positive <- result[seq(1, nrow(result), by = 2), ]
  negative <- result[seq(2, nrow(result), by = 2), ]
  positive$Scale_factor * (positive$Estimate - negative$Estimate)
}


test_that("displayed differences reconstruct lm coefficients for every coding", {

  data <- simulate_one_predictor_data()

  codings <- list(
    treatment = stats::contr.treatment(4),
    treatment_base_3 = stats::contr.treatment(4, base = 3),
    sum = stats::contr.sum(4),
    sum_first_omitted = rbind(c(-1, -1, -1), diag(3)),
    helmert = stats::contr.helmert(4),
    reverse_helmert = reverse_helmert(4),
    sliding = MASS::contr.sdif(4),
    custom = cbind(c(-3, 1, 2, 0) / 6, c(0, -1, 1, 0) / 2, c(1, 1, 1, -3))
  )

  for (coding_name in names(codings)) {

    stats::contrasts(data$group) <- codings[[coding_name]]

    model <- stats::lm(y ~ group + group_size, data = data)

    result <- quiet_contrasts(
      model = model,
      predictor = "group",
      outcome = "y",
      data = data,
      estimate_type = "predicted"
    )

    expect_equal(
      displayed_coefficients(result),
      # The three group coefficients follow the intercept.
      unname(stats::coef(model)[2:4]),
      tolerance = 1e-8,
      info = coding_name
    )
  }
})


test_that("coefficients are matched by term, not by name prefix", {

  data <- simulate_one_predictor_data(k = 3)

  model <- stats::lm(y ~ group_size + group, data = data)

  result <- quiet_contrasts(
    model = model,
    predictor = "group",
    outcome = "y",
    data = data,
    estimate_type = "predicted"
  )

  p_values <- summary(model)$coefficients[c("groupB", "groupC"), "Pr(>|t|)"]

  expect_identical(
    result$Significance[c(1, 3)],
    ifelse(p_values < 0.05, "sig", "notsig"),
    ignore_attr = TRUE
  )
})


test_that("treatment contrasts with a non-first reference are returned", {

  data <- simulate_one_predictor_data(k = 3)
  stats::contrasts(data$group) <- stats::contr.treatment(3, base = 3)

  result <- quiet_contrasts(
    model = stats::lm(y ~ group, data = data),
    predictor = "group",
    outcome = "y",
    data = data,
    estimate_type = "observed"
  )

  expect_identical(result$Levels, c("A", "C", "B", "C"))
  expect_identical(attr(result, "contrast_type"), "Treatment")
})


test_that("sum contrasts with any omitted level are compared with the Grand Mean", {

  data <- simulate_one_predictor_data(k = 3)
  stats::contrasts(data$group) <- rbind(c(-1, -1), c(1, 0), c(0, 1))

  model <- stats::lm(y ~ group, data = data)

  result <- quiet_contrasts(
    model = model,
    predictor = "group",
    outcome = "y",
    data = data,
    estimate_type = "observed"
  )

  expect_identical(result$Levels, c("B", "Grand Mean", "C", "Grand Mean"))
  expect_equal(displayed_coefficients(result), unname(stats::coef(model)[-1]))
  expect_identical(attr(result, "contrast_type"), "Sum")
})


test_that("non-orthogonal custom contrasts show weights and a divergence note", {

  data <- simulate_one_predictor_data(k = 3)
  stats::contrasts(data$group) <- cbind(c(-3, 1, 2) / 6, c(0, -1, 1) / 2)

  expect_message(
    result <- contrasts_onepred(
      model = stats::lm(y ~ group, data = data),
      predictor = "group",
      outcome = "y",
      data = data,
      estimate_type = "observed"
    ),
    "suggest C vs B"
  )

  expect_identical(result$Levels[3:4], c("W(A .20, C .80)", "B"))
  expect_equal(unique(result$Scale_factor), c(4 / 3, 10 / 9))
  expect_match(attr(result, "disclaimer"), "derived from a hypothesis matrix")
})


test_that("single-level estimates keep their link-scale confidence intervals", {

  data <- simulate_one_predictor_data()
  stats::contrasts(data$group) <- stats::contr.helmert(4)

  result <- quiet_contrasts(
    model = stats::glm(success ~ group + group_size, family = stats::binomial, data = data),
    predictor = "group",
    outcome = "success",
    data = data,
    estimate_type = "predicted",
    prediction_scale = "response"
  )

  single_levels <- result$Levels %in% LETTERS[1:4]

  expect_true(all(result$CI_low[single_levels] >= 0))
  expect_true(all(result$CI_high[single_levels] <= 1))
})


test_that("glm coefficients are reconstructed on the link scale", {

  data <- simulate_one_predictor_data()
  stats::contrasts(data$group) <- stats::contr.helmert(4)

  model <- stats::glm(success ~ group + group_size, family = stats::binomial, data = data)

  result <- quiet_contrasts(
    model = model,
    predictor = "group",
    outcome = "success",
    data = data,
    estimate_type = "predicted",
    prediction_scale = "link"
  )

  expect_equal(displayed_coefficients(result), unname(stats::coef(model)[2:4]), tolerance = 1e-8)
})


test_that("an interaction with a numeric predictor triggers a warning", {

  data <- simulate_one_predictor_data(k = 3)

  expect_warning(
    quiet_contrasts(
      model = stats::lm(y ~ group * group_size, data = data),
      predictor = "group",
      outcome = "y",
      data = data,
      estimate_type = "predicted"
    ),
    "do not reconstruct"
  )
})


test_that("contrasts in data must match the fitted model", {

  data <- simulate_one_predictor_data(k = 3)
  model <- stats::lm(y ~ group, data = data)

  # Different data with an observed level that the model never saw.
  data$group <- factor(data$group, levels = c("A", "B", "C", "D"))
  data$group[1:5] <- "D"
  rownames(data) <- paste0("other_", seq_len(nrow(data)))

  expect_error(
    quiet_contrasts(model = model, predictor = "group", outcome = "y", data = data),
    "must be those used to fit the model"
  )
})


test_that("polynomial output carries the same attributes as other codings", {

  data <- simulate_one_predictor_data()
  data$ordered_group <- ordered(data$group)

  result <- quiet_contrasts(
    model = stats::lm(y ~ ordered_group, data = data),
    predictor = "ordered_group",
    outcome = "y",
    data = data,
    estimate_type = "predicted"
  )

  expect_identical(attr(result, "contrast_type"), "Polynomial")
  expect_identical(attr(result, "estimate_type"), "predicted")
  expect_true(all(c("Weight", "Polynomial_fit", "Display_Levels") %in% names(result)))
})


test_that("model_contrast_table() uses the shared classification", {

  data <- simulate_one_predictor_data(k = 3)
  stats::contrasts(data$group) <- rbind(c(-1, -1), c(1, 0), c(0, 1))

  variables <- model_contrast_table(stats::lm(y ~ group, data = data), type = "x")

  expect_identical(
    variables$Contrast_Matrix_Type[variables$Variable == "group"],
    "Sum"
  )
})


test_that("other factors are averaged or held where the coefficients are defined", {

  set.seed(31)

  data <- expand.grid(
    group = factor(c("A", "B", "C")),
    other = factor(c("u", "v")),
    replicate = 1:25
  )
  data$x <- stats::rnorm(nrow(data), mean = data$other == "v")
  data$y <- stats::rnorm(
    nrow(data),
    as.integer(data$group) + 2 * (data$other == "v") + (data$group == "C") * (data$other == "v") + data$x
  )

  for (coding in list(stats::contr.treatment(2), stats::contr.sum(2))) {

    stats::contrasts(data$other) <- coding

    model <- stats::lm(y ~ group * other + x, data = data)

    expect_silent(
      result <- quiet_contrasts(
        model = model,
        predictor = "group",
        outcome = "y",
        data = data,
        estimate_type = "predicted"
      )
    )

    expect_equal(displayed_coefficients(result), unname(stats::coef(model)[2:3]), tolerance = 1e-8)
  }

  # Without an interaction, the other factor is averaged with equal weights.
  model <- stats::lm(y ~ group + other, data = data)

  result <- quiet_contrasts(
    model = model,
    predictor = "group",
    outcome = "y",
    data = data,
    estimate_type = "predicted"
  )

  balanced <- stats::predict(
    model,
    newdata = expand.grid(group = factor("A", levels = c("A", "B", "C")), other = factor(c("u", "v")))
  )

  expect_equal(result$Estimate[result$Levels == "A"][1], mean(balanced))
})


test_that("response-scale intervals stay within the range of the response", {

  set.seed(41)

  data <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 40)),
    other = factor(rep(c("u", "v"), 60)),
    z = stats::rnorm(120)
  )
  data$success <- stats::rbinom(120, 1, c(0.03, 0.02, 0.5)[data$group])
  data$count <- stats::rpois(120, c(0.2, 0.1, 3)[data$group])
  stats::contrasts(data$group) <- stats::contr.helmert(3)

  binomial_result <- quiet_contrasts(
    model = stats::glm(success ~ group + other, family = stats::binomial, data = data),
    predictor = "group",
    outcome = "success",
    data = data,
    estimate_type = "predicted"
  )

  expect_true(all(binomial_result$CI_low > 0 & binomial_result$CI_high < 1))
  expect_true(all(binomial_result$CI_low < binomial_result$Estimate))
  expect_true(all(binomial_result$CI_high > binomial_result$Estimate))

  poisson_result <- quiet_contrasts(
    model = stats::glm(count ~ group + z, family = stats::poisson, data = data),
    predictor = "group",
    outcome = "count",
    data = data,
    estimate_type = "predicted"
  )

  expect_true(all(poisson_result$CI_low > 0))
})


test_that("the link-scale interval of an average reduces to the single-prediction interval", {

  set.seed(42)

  data <- data.frame(group = factor(rep(c("A", "B"), each = 50)), z = stats::rnorm(100))
  data$success <- stats::rbinom(100, 1, c(0.2, 0.6)[data$group])

  model <- stats::glm(success ~ group + z, family = stats::binomial, data = data)

  predictions <- .compute_fixed_predictions(
    model,
    data.frame(group = factor("B", levels = c("A", "B")), z = 0),
    prediction_scale = "response"
  )

  combined <- .combine_fixed_predictions(predictions, weights = 1, prediction_scale = "response")

  expect_equal(combined$ci_low, predictions$ci_low)
  expect_equal(combined$ci_high, predictions$ci_high)
})


test_that("intervals are ordered for decreasing links", {

  set.seed(43)

  data <- data.frame(group = factor(rep(c("A", "B", "C"), each = 40)), z = stats::rnorm(120))
  data$y <- stats::rgamma(120, shape = 2, rate = 2 / c(1, 2, 4)[data$group])

  result <- quiet_contrasts(
    model = stats::glm(y ~ group + z, family = stats::Gamma, data = data),
    predictor = "group",
    outcome = "y",
    data = data,
    estimate_type = "predicted"
  )

  expect_true(all(result$CI_low < result$Estimate & result$Estimate < result$CI_high))
})


test_that("observed proportions use Wilson intervals within [0, 1]", {

  data <- data.frame(group = factor(rep(c("A", "B", "C"), each = 40)))
  data$success <- c(rep(0, 40), rep(c(1, 0, 0, 0), 10), rep(c(1, 0), 20))
  stats::contrasts(data$group) <- stats::contr.helmert(3)

  result <- quiet_contrasts(
    model = stats::glm(success ~ group, family = stats::binomial, data = data),
    predictor = "group",
    outcome = "success",
    data = data,
    estimate_type = "observed"
  )

  expect_true(all(result$CI_low >= 0 & result$CI_high <= 1))

  level_a <- result[result$Levels == "A", ][1, ]
  expect_identical(level_a$CI_low, 0)
  expect_gt(level_a$CI_high, 0)

  level_b <- result[result$Levels == "B", ][1, ]
  expect_equal(
    c(level_b$CI_low, level_b$CI_high),
    as.numeric(stats::prop.test(10, 40, correct = FALSE)$conf.int)
  )
})


test_that("observed counts use positive intervals", {

  data <- data.frame(group = factor(rep(c("A", "B"), each = 30)))
  data$count <- c(rep(0, 30), stats::rpois(30, 2))

  result <- quiet_contrasts(
    model = stats::glm(count ~ group, family = stats::poisson, data = data),
    predictor = "group",
    outcome = "count",
    data = data,
    estimate_type = "observed"
  )

  level_a <- result[result$Levels == "A", ][1, ]

  expect_identical(level_a$CI_low, 0)
  expect_equal(level_a$CI_high, -log(0.025) / 30)
  expect_true(all(result$CI_low >= 0))
})


test_that("observed intervals of unbounded responses are unchanged", {

  data <- simulate_one_predictor_data(k = 3)

  result <- quiet_contrasts(
    model = stats::lm(y ~ group, data = data),
    predictor = "group",
    outcome = "y",
    data = data,
    estimate_type = "observed"
  )

  expect_equal(result$CI_high - result$Estimate, result$Estimate - result$CI_low)
})


test_that("estimate_type = 'auto' follows the model and the requested scale", {

  set.seed(61)

  data <- data.frame(
    group = factor(rep(c("A", "B", "C"), 60)),
    id = factor(rep(1:30, each = 6)),
    z = stats::rnorm(180)
  )
  data$y <- stats::rnorm(180, as.integer(data$group))
  data$success <- stats::rbinom(180, 1, 0.4)

  resolved <- function(model, outcome, prediction_scale = "response") {
    attr(
      suppressWarnings(
        quiet_contrasts(
          model = model,
          predictor = "group",
          outcome = outcome,
          data = data,
          prediction_scale = prediction_scale
        )
      ),
      "estimate_type"
    )
  }

  glm_model <- stats::glm(success ~ group, family = stats::binomial, data = data)

  expect_identical(resolved(stats::lm(y ~ group, data = data), "y"), "observed")
  expect_identical(resolved(glm_model, "success"), "observed")
  expect_identical(resolved(glm_model, "success", prediction_scale = "link"), "predicted")
  expect_identical(resolved(stats::lm(y ~ group + z, data = data), "y"), "predicted")
  expect_identical(resolved(lme4::lmer(y ~ group + (1 | id), data = data), "y"), "predicted")
})


test_that("nlme::lme models are supported", {

  set.seed(62)

  data <- data.frame(
    group = factor(rep(c("A", "B", "C"), 60)),
    id = factor(rep(1:30, each = 6)),
    unused = stats::rnorm(180)
  )
  data$y <- stats::rnorm(180, as.integer(data$group)) + stats::rnorm(30)[data$id]
  stats::contrasts(data$group) <- stats::contr.helmert(3)

  model <- nlme::lme(y ~ group, random = ~ 1 | id, data = data)

  variables <- model_contrast_table(model, type = "x&y")

  expect_setequal(variables$Variable[variables$Role == "Predictor"], c("group", "id"))
  expect_identical(variables$Variable[variables$Role == "Dependent"], "y")

  result <- quiet_contrasts(
    model = model,
    predictor = "group",
    outcome = "y",
    data = data,
    estimate_type = "predicted"
  )

  expect_equal(displayed_coefficients(result), unname(nlme::fixef(model)[2:3]), tolerance = 1e-8)
})


simulate_ordered_data <- function() {
  set.seed(71)
  data <- data.frame(
    dose = ordered(rep(c("low", "mid", "high"), each = 50), levels = c("low", "mid", "high")),
    z = stats::rnorm(150, mean = 5)
  )
  data$y <- stats::rnorm(150, c(1, 3, 4)[as.integer(data$dose)] + 0.5 * data$z)
  data$success <- stats::rbinom(150, 1, c(0.1, 0.4, 0.5)[as.integer(data$dose)])
  data
}

polynomial_components <- function(result) {
  split(result, result$Contrast)
}


test_that("polynomial curves are centred on the displayed estimates", {

  data <- simulate_ordered_data()

  result <- quiet_contrasts(
    model = stats::lm(y ~ dose + z, data = data),
    predictor = "dose",
    outcome = "y",
    data = data,
    estimate_type = "predicted"
  )

  for (component in polynomial_components(result)) {
    expect_equal(mean(component$Polynomial_fit), mean(component$Estimate))
  }

  expect_identical(result$Weight[result$Levels == "mid" & endsWith(result$Contrast, ".L")], 0)
})


test_that("polynomial components reconstruct the level estimates on the link scale", {

  data <- simulate_ordered_data()

  result <- quiet_contrasts(
    model = stats::glm(success ~ dose, family = stats::binomial, data = data),
    predictor = "dose",
    outcome = "success",
    data = data,
    estimate_type = "predicted",
    prediction_scale = "link"
  )

  components <- polynomial_components(result)
  estimates <- components[[1]]$Estimate
  baseline <- mean(estimates)

  reconstructed <- baseline + Reduce(
    `+`,
    lapply(components, function(component) component$Polynomial_fit - baseline)
  )

  expect_equal(reconstructed, estimates)
})


test_that("polynomial curves are on the response scale with response-scale estimates", {

  data <- simulate_ordered_data()

  result <- quiet_contrasts(
    model = stats::glm(success ~ dose, family = stats::binomial, data = data),
    predictor = "dose",
    outcome = "success",
    data = data,
    estimate_type = "predicted",
    prediction_scale = "response"
  )

  expect_true(all(result$Polynomial_fit > 0 & result$Polynomial_fit < 1))
})


test_that("models without an intercept are rejected with an explanation", {

  data <- simulate_ordered_data()

  expect_error(
    quiet_contrasts(
      model = stats::lm(y ~ 0 + dose, data = data),
      predictor = "dose",
      outcome = "y",
      data = data
    ),
    "requires a model with an intercept"
  )
})


simulate_outcome_data <- function() {
  set.seed(81)
  data <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 40)),
    x = stats::rexp(120) + 0.1
  )
  data$y <- exp(stats::rnorm(120, as.integer(data$group) / 2))
  data$answer <- factor(
    ifelse(stats::runif(120) < c(0.2, 0.4, 0.6)[data$group], "yes", "no"),
    levels = c("no", "yes")
  )
  data$successes <- stats::rbinom(120, 5, 0.4)
  data
}


test_that("factor and logical outcomes of binomial models are converted to 0/1", {

  data <- simulate_outcome_data()

  factor_result <- quiet_contrasts(
    model = stats::glm(answer ~ group, family = stats::binomial, data = data),
    predictor = "group",
    outcome = "answer",
    data = data,
    estimate_type = "observed"
  )

  observed_proportions <- tapply(data$answer == "yes", data$group, mean)

  expect_equal(
    factor_result$Estimate[match(c("A", "B", "C"), factor_result$Levels)],
    as.numeric(observed_proportions)
  )

  data$answer_logical <- data$answer == "yes"

  logical_result <- quiet_contrasts(
    model = stats::glm(answer_logical ~ group, family = stats::binomial, data = data),
    predictor = "group",
    outcome = "answer_logical",
    data = data,
    estimate_type = "observed"
  )

  expect_equal(logical_result$Estimate, factor_result$Estimate)
})


test_that("observed estimates use the transformed response of the model", {

  data <- simulate_outcome_data()

  result <- quiet_contrasts(
    model = stats::lm(log(y) ~ group, data = data),
    predictor = "group",
    outcome = "y",
    data = data,
    estimate_type = "observed"
  )

  expect_equal(
    result$Estimate[result$Levels == "A"][1],
    mean(log(data$y[data$group == "A"]))
  )
})


test_that("two-column binomial responses support predicted estimates only", {

  data <- simulate_outcome_data()

  model <- stats::glm(cbind(successes, 5 - successes) ~ group, family = stats::binomial, data = data)

  expect_error(
    quiet_contrasts(model = model, predictor = "group", outcome = "successes", data = data, estimate_type = "observed"),
    "several columns"
  )

  result <- quiet_contrasts(
    model = model,
    predictor = "group",
    outcome = "successes",
    data = data,
    estimate_type = "predicted"
  )

  expect_true(all(result$Estimate > 0 & result$Estimate < 1))
})


test_that("the analysis data are recovered with fitted rows and contrasts", {

  data <- simulate_outcome_data()
  data$y[1:3] <- NA
  stats::contrasts(data$group) <- stats::contr.sum(3)

  model <- stats::lm(y ~ group + log(x), data = data)

  # The data are modified after fitting.
  stats::contrasts(data$group) <- stats::contr.treatment(3)

  analysis_data <- .extract_analysis_data(model)

  expect_true("x" %in% names(analysis_data))
  expect_identical(nrow(analysis_data), 117L)
  expect_equal(unname(stats::contrasts(analysis_data$group)), unname(stats::contr.sum(3)))
})


test_that("contrasts_onepred() recovers data and outcome from the model", {

  data <- simulate_one_predictor_data()
  stats::contrasts(data$group) <- stats::contr.helmert(4)
  model <- stats::lm(y ~ group, data = data)

  implicit <- quiet_contrasts(model, "group")
  explicit <- quiet_contrasts(model, "group", data = data, outcome = "y")

  expect_equal(implicit, explicit)
})


test_that("supplied data are restricted to the rows used for fitting", {

  data <- simulate_one_predictor_data()
  data$y[c(3, 50)] <- NA
  model <- stats::lm(y ~ group, data = data, subset = group_size > 0)

  implicit <- quiet_contrasts(model, "group")
  explicit <- quiet_contrasts(model, "group", data = data)

  expect_equal(explicit, implicit)
  expect_equal(displayed_coefficients(explicit), unname(stats::coef(model)[-1]))
  expect_false(grepl("do not reconstruct", attr(explicit, "disclaimer")))
})


test_that("contrasts given in the model call are used", {

  data <- simulate_one_predictor_data()
  data$cluster <- factor(rep(seq_len(10), length.out = nrow(data)))

  models <- list(
    lm = stats::lm(y ~ group, data = data, contrasts = list(group = "contr.sum")),
    lmer = lme4::lmer(y ~ group + (1 | cluster), data = data, contrasts = list(group = MASS::contr.sdif(4)))
  )

  for (model_name in names(models)) {

    model <- models[[model_name]]

    implicit <- quiet_contrasts(model, "group", estimate_type = "predicted")
    explicit <- quiet_contrasts(model, "group", data = data, estimate_type = "predicted")

    expect_equal(explicit, implicit, info = model_name)
    expect_equal(
      displayed_coefficients(implicit),
      unname(.extract_fixef(model)[-1]),
      tolerance = 1e-6,
      info = model_name
    )
  }

  expect_identical(
    model_contrast_table(models$lm, type = "x")$Contrast_Matrix_Type,
    "Sum"
  )
})


test_that("factor levels without observations are dropped, as in fitting", {

  data <- simulate_one_predictor_data(k = 3)
  data$group <- factor(data$group, levels = c("A", "B", "C", "Z"))

  model <- stats::lm(y ~ group, data = data)

  implicit <- quiet_contrasts(model, "group")
  explicit <- quiet_contrasts(model, "group", data = data)

  expect_equal(explicit, implicit)
  expect_setequal(implicit$Levels, c("A", "B", "C"))
  expect_equal(displayed_coefficients(implicit), unname(stats::coef(model)[-1]))

  # A level emptied by 'subset' is dropped as well.
  emptied <- stats::lm(y ~ group, data = data, subset = group != "C")

  expect_setequal(quiet_contrasts(emptied, "group", data = data)$Levels, c("A", "B"))
})


test_that("auto uses model-based estimates with unequal prior weights", {

  data <- simulate_one_predictor_data()
  data$w <- stats::runif(nrow(data), 0.2, 3)

  weighted <- stats::lm(y ~ group, data = data, weights = w)
  result <- quiet_contrasts(weighted, "group")

  expect_match(attr(result, "disclaimer"), "^Estimates are model-based")
  expect_equal(displayed_coefficients(result), unname(stats::coef(weighted)[-1]))

  # Constant weights leave the estimates unchanged.
  data$w <- 2
  constant <- stats::lm(y ~ group, data = data, weights = w)

  expect_match(attr(quiet_contrasts(constant, "group"), "disclaimer"), "^Estimates are observed")
  expect_match(attr(quiet_contrasts(stats::glm(success ~ group, binomial, data), "group"), "disclaimer"), "^Estimates are observed")
})


test_that("divergence notes are printed once and kept in the disclaimer", {

  data <- simulate_one_predictor_data(k = 3)
  stats::contrasts(data$group) <- cbind(c(-1, 0.5, 0.5), c(-1, 0, 1))
  model <- stats::lm(y ~ group, data = data)

  messages <- character()
  result <- withCallingHandlers(
    contrasts_onepred(model, "group"),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  printed <- paste(messages, collapse = " ")
  occurrences <- lengths(regmatches(printed, gregexpr("coding column suggest", printed, fixed = TRUE)))

  expect_identical(occurrences, 2L)
  expect_match(attr(result, "disclaimer"), "coding column suggest")
})
