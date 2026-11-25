# Declared hypotheses: the comparisons that the coefficients are intended
# to estimate, checked against the comparisons implied by the coding.

trial_data <- function() {
  set.seed(2)
  data <- data.frame(group = factor(rep(c("Control", "T1", "T2"), each = 30), levels = c("Control", "T1", "T2")))
  data$score <- stats::rnorm(90, mean = c(10, 12, 15)[data$group], sd = 3)
  data
}

intended <- function() {
  rbind("Treat-Control" = c(Control = -1, T1 = 0.5, T2 = 0.5),
        "T2-Control" = c(Control = -1, T1 = 0, T2 = 1))
}

collect_messages <- function(expression) {
  messages <- character()
  value <- withCallingHandlers(
    expression,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, messages = messages)
}


test_that("hypotheses written into the coding columns are reported as not estimated", {

  data <- trial_data()
  stats::contrasts(data$group) <- t(intended())

  result <- collect_messages(
    contrasts_onepred(stats::lm(score ~ group, data = data), "group", hypotheses = intended())
  )

  check <- attr(result$value, "hypothesis_check")

  expect_identical(check$Status, c("differs", "differs"))
  expect_identical(check$Hypothesis, c("Treat-Control", "T2-Control"))
  expect_true(any(grepl("does not estimate the declared hypothesis 'Treat-Control'", result$messages)))
  expect_true(any(grepl("it estimates 4/3 x [T1 - Mean(Control, T2)]", result$messages, fixed = TRUE)))

  # The exact check replaces the note based on the signs of the coding.
  expect_false(any(grepl("signs of the coding column", result$messages)))
})


test_that("a coding derived from the hypothesis matrix is confirmed", {

  data <- trial_data()
  coding <- MASS::ginv(intended())
  dimnames(coding) <- list(levels(data$group), rownames(intended()))
  stats::contrasts(data$group) <- coding

  model <- stats::lm(score ~ group, data = data)
  result <- collect_messages(contrasts_onepred(model, "group", hypotheses = intended()))

  expect_identical(attr(result$value, "hypothesis_check")$Status, c("matches", "matches"))
  expect_true(all(grepl("estimates the declared hypothesis", result$messages[grepl("^Contrast_", result$messages)])))

  # The model is read, never modified.
  means <- as.numeric(tapply(data$score, data$group, mean))
  expect_equal(unname(stats::coef(model)[-1]), c(mean(means[2:3]) - means[1], means[3] - means[1]))
})


test_that("rescaled and reversed hypotheses are recognized", {

  data <- trial_data()
  stats::contrasts(data$group) <- stats::contr.helmert(3)

  model <- stats::lm(score ~ group, data = data)

  rescaled <- attr(
    suppressMessages(contrasts_onepred(model, "group", hypotheses = rbind(c(-1, 1, 0), c(-0.5, -0.5, 1)))),
    "hypothesis_check"
  )

  expect_identical(rescaled$Status, c("rescaled", "rescaled"))
  expect_equal(rescaled$Ratio, c(1 / 2, 1 / 3))
  expect_identical(rescaled$Hypothesis, c("H1", "H2"))

  reversed <- attr(
    suppressMessages(contrasts_onepred(model, "group", hypotheses = rbind(c(1, -1, 0), c(-0.5, -0.5, 1)))),
    "hypothesis_check"
  )

  expect_identical(reversed$Status, c("reversed", "rescaled"))
})


test_that("declared hypotheses are validated and reordered by level", {

  data <- trial_data()
  stats::contrasts(data$group) <- stats::contr.treatment(3)
  model <- stats::lm(score ~ group, data = data)

  expect_error(
    contrasts_onepred(model, "group", hypotheses = rbind(c(-1, 1, 0))),
    "one row per coefficient"
  )

  expect_error(
    contrasts_onepred(model, "group", hypotheses = rbind(c(1, 1, 0), c(-1, 0, 1))),
    "summing to zero"
  )

  expect_error(
    contrasts_onepred(model, "group", hypotheses = rbind(c(A = -1, B = 1, C = 0), c(A = -1, B = 0, C = 1))),
    "column names"
  )

  # Columns given in another order are matched by level name.
  shuffled <- rbind(c(T2 = 0, Control = -1, T1 = 1), c(T2 = 1, Control = -1, T1 = 0))

  check <- attr(suppressMessages(contrasts_onepred(model, "group", hypotheses = shuffled)), "hypothesis_check")

  expect_identical(check$Status, c("matches", "matches"))
})


test_that("the caption reports the result of the check", {

  data <- trial_data()

  caption_of <- function(model, ...) {
    output <- capture.output(suppressMessages(moplot_onepred(model, "group", ...)))
    paste(trimws(output), collapse = " ")
  }

  stats::contrasts(data$group) <- t(intended())
  written_into_coding <- stats::lm(score ~ group, data = data)

  expect_match(caption_of(written_into_coding, hypotheses = intended()), "C1 and C2 do not estimate the declared hypotheses")
  expect_match(caption_of(written_into_coding), "signs of the coding columns of C1 and C2")

  coding <- MASS::ginv(intended())
  dimnames(coding) <- list(levels(data$group), rownames(intended()))
  stats::contrasts(data$group) <- coding

  expect_match(
    caption_of(stats::lm(score ~ group, data = data), hypotheses = intended()),
    "Each coefficient estimates the corresponding declared hypothesis"
  )
})


# ============================================================
# Interactions
# ============================================================

interaction_cells <- function() {
  cells <- data.frame(
    treatment = factor(rep(c("Low", "Medium", "High"), 2), levels = c("Low", "Medium", "High")),
    context = factor(rep(c("Baseline", "Enhanced"), each = 3)),
    mu = c(5, 7, 3, 6, 10, 5)
  )
  data <- cells[rep(seq_len(6), each = 5), ]
  data$outcome <- data$mu + rep(c(-1, -0.5, 0, 0.5, 1), 6)
  data
}

treatment_hypotheses <- function() {
  rbind("Medium-Low" = c(Low = -1, Medium = 1, High = 0),
        "High-Medium" = c(Low = 0, Medium = -1, High = 1))
}

context_hypotheses <- function() {
  rbind("Enhanced-Baseline" = c(Baseline = -1, Enhanced = 1))
}


test_that("interaction coefficients inherit the check of their components", {

  data <- interaction_cells()
  stats::contrasts(data$treatment) <- MASS::ginv(treatment_hypotheses())
  stats::contrasts(data$context) <- MASS::contr.sdif(2)

  result <- collect_messages(
    contrasts_multpred(
      stats::lm(outcome ~ treatment * context, data = data),
      data = data,
      hypotheses = list(treatment = treatment_hypotheses(), context = context_hypotheses()),
      methodological_note = FALSE
    )
  )

  check <- attr(result$value, "hypothesis_check")

  expect_identical(check$Status, c("matches", "matches"))
  expect_identical(check$Hypothesis, c("Medium-Low x Enhanced-Baseline", "High-Medium x Enhanced-Baseline"))
  expect_true(any(grepl("Variable treatment, contrast 1 estimates the declared hypothesis", result$messages, fixed = TRUE)))
})


test_that("interaction coefficients are reported when a component differs", {

  data <- interaction_cells()
  stats::contrasts(data$treatment) <- t(treatment_hypotheses())
  stats::contrasts(data$context) <- MASS::contr.sdif(2)

  result <- collect_messages(
    contrasts_multpred(
      stats::lm(outcome ~ treatment * context, data = data),
      data = data,
      hypotheses = list(treatment = treatment_hypotheses(), context = context_hypotheses()),
      methodological_note = FALSE
    )
  )

  expect_identical(attr(result$value, "hypothesis_check")$Status, c("differs", "differs"))
  expect_true(any(grepl("does not estimate the declared hypothesis 'Medium-Low'", result$messages, fixed = TRUE)))
})


test_that("undeclared components are reported as not checked", {

  data <- interaction_cells()
  stats::contrasts(data$treatment) <- MASS::ginv(treatment_hypotheses())

  check <- attr(
    suppressMessages(
      contrasts_multpred(
        stats::lm(outcome ~ treatment * context, data = data),
        data = data,
        hypotheses = list(treatment = treatment_hypotheses()),
        methodological_note = FALSE
      )
    ),
    "hypothesis_check"
  )

  expect_identical(check$Status, c("partially declared", "partially declared"))
  expect_match(check$Hypothesis[1], "context (not declared)", fixed = TRUE)
})


test_that("hypotheses for interactions are validated", {

  data <- interaction_cells()
  model <- stats::lm(outcome ~ treatment * context, data = data)

  expect_error(
    contrasts_multpred(model, data = data, hypotheses = treatment_hypotheses(), methodological_note = FALSE),
    "named list"
  )

  expect_error(
    contrasts_multpred(model, data = data, hypotheses = list(dose = treatment_hypotheses()), methodological_note = FALSE),
    "not categorical predictors"
  )
})


test_that("contrasts_onepred() also accepts a named list", {

  data <- trial_data()
  stats::contrasts(data$group) <- t(intended())

  check <- attr(
    suppressMessages(
      contrasts_onepred(stats::lm(score ~ group, data = data), "group", hypotheses = list(group = intended()))
    ),
    "hypothesis_check"
  )

  expect_identical(check$Status, c("differs", "differs"))
})


# ============================================================
# Layout of shared quantities
# ============================================================

test_that("a level shared on the same side of two contrasts gets its own pairs", {

  data <- trial_data()
  coding <- MASS::ginv(intended())
  dimnames(coding) <- list(levels(data$group), rownames(intended()))
  stats::contrasts(data$group) <- coding

  # Both comparisons have Control on their negative side.
  same_side <- suppressMessages(.prepare_moplot_onepred(stats::lm(score ~ group, data = data), "group"))

  expect_identical(same_side$position_layout, "paired")

  # Under sliding differences the shared level lies on opposite sides.
  stats::contrasts(data$group) <- MASS::contr.sdif(3)

  chain <- suppressMessages(.prepare_moplot_onepred(stats::lm(score ~ group, data = data), "group"))

  expect_identical(chain$position_layout, "reuse")
})
