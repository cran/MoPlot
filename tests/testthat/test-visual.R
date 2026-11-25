# Visual regression tests. vdiffr compares each figure with a reviewed
# SVG snapshot in tests/testthat/_snaps/visual/ and skips these tests on
# CRAN, where font rendering differs between systems.

expect_plot <- function(title, plot) {
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(title, function() print(plot))
}

quiet_plot <- function(expression) {
  suppressMessages(suppressWarnings(expression))
}

one_predictor_data <- function() {
  set.seed(123)
  data <- data.frame(group = factor(rep(c("A", "B", "C", "D"), each = 30)))
  data$y <- stats::rnorm(120, c(0, 0.4, 0.9, 1)[data$group], 0.9)
  data
}

interaction_data <- function() {
  set.seed(21)
  data <- expand.grid(
    A = factor(c("a1", "a2", "a3")),
    B = factor(c("b1", "b2")),
    C = factor(c("c1", "c2")),
    replicate = 1:15
  )
  data$y <- stats::rnorm(
    nrow(data),
    (data$A == "a2") + 1.5 * (data$A == "a3") * (data$B == "b2") + (data$C == "c2") * (data$B == "b2")
  )
  stats::contrasts(data$A) <- stats::contr.helmert(3)
  stats::contrasts(data$B) <- stats::contr.sum(2)
  stats::contrasts(data$C) <- stats::contr.sum(2)
  data
}


# ============================================================
# moplot_onepred()
# ============================================================

test_that("moplot_onepred() figures are stable", {

  data <- one_predictor_data()

  stats::contrasts(data$group) <- stats::contr.helmert(4)
  helmert <- stats::lm(y ~ group, data = data)

  expect_plot(
    "onepred helmert minimal",
    quiet_plot(moplot_onepred(helmert, "group", caption = FALSE))
  )

  expect_plot(
    "onepred helmert per difference",
    quiet_plot(moplot_onepred(helmert, "group", per_difference = TRUE, caption = FALSE))
  )

  stats::contrasts(data$group) <- stats::contr.sum(4)

  expect_plot(
    "onepred sum publication with data",
    quiet_plot(
      moplot_onepred(
        stats::lm(y ~ group, data = data),
        "group",
        style = "publication",
        show_data = TRUE,
        caption = FALSE
      )
    )
  )

  custom_data <- data.frame(group = factor(rep(c("Control", "Treatment1", "Treatment2"), each = 40)))
  set.seed(5)
  custom_data$y <- stats::rnorm(120, c(1, 2, 4)[custom_data$group])
  stats::contrasts(custom_data$group) <- cbind(c(-3, 1, 2) / 6, c(0, -1, 1) / 2)

  expect_plot(
    "onepred custom weighted labels",
    quiet_plot(moplot_onepred(stats::lm(y ~ group, data = custom_data), "group", caption = FALSE))
  )
})


test_that("moplot_onepred() polynomial figures are stable", {

  data <- one_predictor_data()

  stats::contrasts(data$group) <- stats::contr.poly(4)

  expect_plot(
    "onepred polynomial publication with data",
    quiet_plot(
      moplot_onepred(
        stats::lm(y ~ group, data = data),
        "group",
        style = "publication",
        show_data = TRUE,
        caption = FALSE
      )
    )
  )

  stats::contrasts(data$group) <- stats::contr.poly(4, scores = c(0, 10, 20, 40))

  expect_plot(
    "onepred polynomial scores minimal",
    quiet_plot(moplot_onepred(stats::lm(y ~ group, data = data), "group", caption = FALSE))
  )

  set.seed(7)
  binary_data <- data.frame(freq = factor(rep(c("a", "b", "c", "d", "e"), each = 60)))
  binary_data$success <- stats::rbinom(300, 1, stats::plogis(c(-1.6, -0.5, 0.4, 0.6, 0.3))[binary_data$freq])
  stats::contrasts(binary_data$freq) <- stats::contr.poly(5)

  expect_plot(
    "onepred polynomial logistic odds ratios",
    quiet_plot(
      moplot_onepred(
        stats::glm(success ~ freq, family = stats::binomial, data = binary_data),
        "freq",
        coefficient_scale = "interpreted",
        style = "publication",
        caption = FALSE
      )
    )
  )
})


# ============================================================
# moplot_multpred()
# ============================================================

test_that("moplot_multpred() figures are stable", {

  data <- interaction_data()

  two_way <- stats::lm(y ~ A * B, data = data)
  three_way <- stats::lm(y ~ A * B * C, data = data)

  expect_plot(
    "multpred two way minimal",
    quiet_plot(moplot_multpred(two_way, data = data, caption = FALSE, interaction_info = FALSE))
  )

  expect_plot(
    "multpred two way tree",
    quiet_plot(moplot_multpred(two_way, data = data, style = "tree", caption = FALSE, interaction_info = FALSE))
  )

  expect_plot(
    "multpred three way tree",
    quiet_plot(moplot_multpred(three_way, data = data, style = "tree", caption = FALSE, interaction_info = FALSE))
  )

  expect_plot(
    "multpred two way per difference",
    quiet_plot(
      moplot_multpred(two_way, data = data, per_difference = TRUE, caption = FALSE, interaction_info = FALSE)
    )
  )
})
