named_contrasts <- function(contrast_matrix, levels = LETTERS[seq_len(nrow(contrast_matrix))]) {
  contrast_matrix <- as.matrix(contrast_matrix)
  rownames(contrast_matrix) <- levels
  contrast_matrix
}

standard_contrasts <- function(k) {
  list(
    "Treatment" = stats::contr.treatment(k),
    "Sum" = stats::contr.sum(k),
    "Helmert" = stats::contr.helmert(k),
    "Reverse Helmert" = reverse_helmert(k),
    "Sliding Differences" = MASS::contr.sdif(k),
    "Polynomial" = stats::contr.poly(k)
  )
}

non_orthogonal_custom <- named_contrasts(
  cbind(c1 = c(-3, 1, 2) / 6, c2 = c(0, -1, 1) / 2)
)


# ============================================================
# HYPOTHESIS MATRIX
# ============================================================

test_that(".hypothesis_matrix() recovers treatment comparisons", {

  hypothesis_matrix <- .hypothesis_matrix(named_contrasts(stats::contr.treatment(3)))

  expect_equal(
    unname(hypothesis_matrix),
    rbind(c(-1, 1, 0), c(-1, 0, 1))
  )

  expect_identical(colnames(hypothesis_matrix), c("A", "B", "C"))
})


test_that(".hypothesis_matrix() rows sum to zero for any full-rank coding", {

  set.seed(1)

  for (k in 2:6) {
    for (repetition in 1:5) {
      hypothesis_matrix <- .hypothesis_matrix(matrix(stats::rnorm(k * (k - 1)), k))
      expect_equal(unname(rowSums(hypothesis_matrix)), rep(0, k - 1))
    }
  }
})


test_that(".hypothesis_matrix() removes numerical residues", {

  hypothesis_matrix <- .hypothesis_matrix(stats::contr.treatment(5, base = 3))

  expect_true(all(hypothesis_matrix %in% c(-1, 0, 1)))
})


test_that(".hypothesis_matrix() stops for rank-deficient or incomplete codings", {

  expect_error(
    .hypothesis_matrix(cbind(c(1, -1, 0), c(2, -2, 0))),
    "not full rank"
  )

  expect_error(
    .hypothesis_matrix(cbind(c(1, 1, 1), c(-1, 0, 1))),
    "not full rank"
  )

  expect_error(
    .hypothesis_matrix(matrix(c(-1, 0, 1, 0), 4)),
    "k - 1"
  )
})


# ============================================================
# CLASSIFICATION
# ============================================================

test_that(".classify_contrast() recognizes canonical matrices", {

  for (k in 3:6) {
    for (contrast_type in names(standard_contrasts(k))) {
      expect_identical(
        .classify_contrast(standard_contrasts(k)[[contrast_type]])$type,
        contrast_type,
        info = paste(contrast_type, k)
      )
    }
  }
})


test_that(".classify_contrast() keeps the canonical labels for two levels", {

  expect_identical(.classify_contrast(stats::contr.treatment(2))$type, "Treatment")
  expect_identical(.classify_contrast(stats::contr.sum(2))$type, "Sum")
  expect_identical(.classify_contrast(stats::contr.helmert(2))$type, "Helmert")
  expect_identical(.classify_contrast(MASS::contr.sdif(2))$type, "Sliding Differences")
})


test_that(".classify_contrast() recognizes treatment contrasts with any reference", {

  for (reference in 1:4) {
    classification <- .classify_contrast(stats::contr.treatment(4, base = reference))
    expect_identical(classification$type, "Treatment")
    expect_identical(classification$reference, reference)
  }

  expect_identical(.classify_contrast(stats::contr.SAS(4))$reference, 4L)

  # Rescaled and reordered treatment columns.
  expect_identical(
    .classify_contrast(2 * stats::contr.treatment(4, base = 2)[, c(3, 1, 2)])$type,
    "Treatment"
  )
})


test_that(".classify_contrast() recognizes sum contrasts with any omitted level", {

  sum_first_omitted <- rbind(c(-1, -1), c(1, 0), c(0, 1))

  expect_identical(.classify_contrast(sum_first_omitted)$type, "Sum")
  expect_identical(.classify_contrast(stats::contr.sum(4)[, c(2, 3, 1)])$type, "Sum")
  expect_identical(.classify_contrast(0.5 * stats::contr.sum(4))$type, "Sum")
})


test_that(".classify_contrast() recognizes rescaled ordered systems", {

  expect_identical(.classify_contrast(6 * stats::contr.helmert(4))$type, "Helmert")
  expect_identical(.classify_contrast(0.5 * reverse_helmert(4))$type, "Reverse Helmert")
  expect_identical(.classify_contrast(2 * MASS::contr.sdif(4))$type, "Sliding Differences")
})


test_that(".classify_contrast() returns Custom for other codings", {

  expect_identical(.classify_contrast(non_orthogonal_custom)$type, "Custom")
  expect_identical(.classify_contrast(cbind(c(1, -1, 0), c(2, -2, 0)))$type, "Custom")
})


# ============================================================
# DECOMPOSITION
# ============================================================

test_that("every decomposition reconstructs its hypothesis row", {

  set.seed(2)

  codings <- c(
    standard_contrasts(4)[names(standard_contrasts(4)) != "Polynomial"],
    list(
      "Custom 1" = non_orthogonal_custom,
      "Random 1" = matrix(stats::rnorm(12), 4),
      "Random 2" = matrix(stats::rnorm(20), 5),
      "Sum first omitted" = rbind(c(-1, -1), c(1, 0), c(0, 1))
    )
  )

  for (coding_name in names(codings)) {

    coding <- named_contrasts(codings[[coding_name]])
    hypotheses <- .contrast_hypotheses(coding)

    for (contrast_index in seq_along(hypotheses$contrasts)) {

      contrast <- hypotheses$contrasts[[contrast_index]]

      expect_gt(contrast$scale, 0)
      expect_equal(sum(contrast$positive_weights), 1)
      expect_equal(sum(contrast$negative_weights), 1)
      expect_true(all(contrast$positive_weights >= 0))
      expect_true(all(contrast$negative_weights >= 0))

      expect_equal(
        contrast$scale * (contrast$positive_weights - contrast$negative_weights),
        hypotheses$hypothesis_matrix[contrast_index, ],
        info = paste(coding_name, contrast_index)
      )
    }
  }
})


test_that("treatment contrasts with a non-first reference are decomposed", {

  hypotheses <- .contrast_hypotheses(named_contrasts(stats::contr.treatment(3, base = 3)))

  labels <- vapply(
    hypotheses$contrasts,
    function(contrast) paste(contrast$positive_label, "vs", contrast$negative_label),
    character(1)
  )

  expect_identical(unname(labels), c("A vs C", "B vs C"))
  expect_equal(unname(vapply(hypotheses$contrasts, `[[`, numeric(1), "scale")), c(1, 1))
})


test_that("sum contrasts are presented as one level versus the Grand Mean", {

  hypotheses <- .contrast_hypotheses(named_contrasts(rbind(c(-1, -1), c(1, 0), c(0, 1))))

  expect_identical(hypotheses$type, "Sum")

  expect_identical(
    unname(vapply(hypotheses$contrasts, `[[`, character(1), "positive_label")),
    c("B", "C")
  )

  expect_true(all(vapply(hypotheses$contrasts, function(x) x$negative_label == "Grand Mean", logical(1))))
  expect_true(all(vapply(hypotheses$contrasts, `[[`, logical(1), "overlapping")))
  expect_equal(unname(vapply(hypotheses$contrasts, `[[`, numeric(1), "scale")), c(1, 1))
})


test_that("Helmert scale factors are made explicit", {

  hypotheses <- .contrast_hypotheses(named_contrasts(stats::contr.helmert(4)))

  expect_equal(
    unname(vapply(hypotheses$contrasts, `[[`, numeric(1), "scale")),
    c(1 / 2, 1 / 3, 1 / 4)
  )

  expect_identical(hypotheses$contrasts[[3]]$negative_label, "Mean(A, B, C)")
})


test_that("polynomial contrasts are classified but not decomposed", {

  hypotheses <- .contrast_hypotheses(named_contrasts(stats::contr.poly(4)))

  expect_identical(hypotheses$type, "Polynomial")
  expect_null(hypotheses$contrasts)
})


test_that("non-orthogonal custom contrasts expose unequal weights and divergence", {

  hypotheses <- .contrast_hypotheses(non_orthogonal_custom)

  first <- hypotheses$contrasts[[1]]
  second <- hypotheses$contrasts[[2]]

  # c1 compares A with the mean of B and C, as its signs suggest.
  expect_identical(first$positive_label, "Mean(B, C)")
  expect_identical(first$negative_label, "A")
  expect_equal(first$scale, 4 / 3)
  expect_false(first$divergent)

  # c2 looks like C vs B, but the coefficient also involves A.
  expect_identical(second$positive_label, "W(A .20, C .80)")
  expect_identical(second$negative_label, "B")
  expect_equal(second$scale, 10 / 9)
  expect_true(second$divergent)
  expect_match(second$note, "suggest C vs B")
})


test_that("recognized types are never flagged as divergent", {

  for (coding in standard_contrasts(4)[names(standard_contrasts(4)) != "Polynomial"]) {
    hypotheses <- .contrast_hypotheses(named_contrasts(coding))
    expect_false(any(vapply(hypotheses$contrasts, `[[`, logical(1), "divergent")))
  }
})


# ============================================================
# IDENTITY WITH FITTED MODELS
# ============================================================

focal_coefficients <- function(model, data, predictor) {
  model_matrix <- stats::model.matrix(stats::formula(model), data = data)
  term_index <- which(attr(stats::terms(model), "term.labels") == predictor)
  stats::coef(model)[attr(model_matrix, "assign") == term_index]
}

simulate_identity_data <- function(k = 4, n_per_level = 30) {
  set.seed(3)
  data <- data.frame(
    group = factor(rep(LETTERS[seq_len(k)], each = n_per_level)),
    covariate = stats::rnorm(k * n_per_level)
  )
  data$y <- stats::rnorm(nrow(data), mean = seq_len(k)[data$group] + 0.5 * data$covariate)
  data$success <- stats::rbinom(nrow(data), 1, stats::plogis(-1 + 0.4 * seq_len(k)[data$group]))
  data
}

link_predictions <- function(model, data, predictor) {
  levels <- levels(data[[predictor]])
  newdata <- data[rep(1L, length(levels)), , drop = FALSE]
  newdata$covariate <- mean(data$covariate)
  newdata[[predictor]] <- factor(levels, levels = levels)
  .compute_fixed_predictions(model, newdata, prediction_scale = "link")$estimate
}

test_that("decompositions reconstruct lm coefficients with a covariate", {

  data <- simulate_identity_data()

  codings <- c(
    standard_contrasts(4),
    list(
      "Treatment base 3" = stats::contr.treatment(4, base = 3),
      "Custom" = cbind(c(-3, 1, 2, 0) / 6, c(0, -1, 1, 0) / 2, c(1, 1, 1, -3))
    )
  )

  for (coding_name in names(codings)) {

    stats::contrasts(data$group) <- codings[[coding_name]]
    model <- stats::lm(y ~ group + covariate, data = data)

    hypotheses <- .contrast_hypotheses(stats::contrasts(data$group), levels(data$group))

    expect_true(
      all(.check_hypothesis_identity(
        hypotheses,
        focal_coefficients(model, data, "group"),
        link_predictions(model, data, "group")
      )),
      info = coding_name
    )
  }
})


test_that("decompositions reconstruct glm coefficients on the link scale", {

  data <- simulate_identity_data()

  for (coding_name in c("Sum", "Helmert", "Sliding Differences")) {

    stats::contrasts(data$group) <- standard_contrasts(4)[[coding_name]]
    model <- stats::glm(success ~ group + covariate, family = stats::binomial, data = data)

    hypotheses <- .contrast_hypotheses(stats::contrasts(data$group), levels(data$group))

    expect_true(
      all(.check_hypothesis_identity(
        hypotheses,
        focal_coefficients(model, data, "group"),
        link_predictions(model, data, "group")
      )),
      info = coding_name
    )
  }
})


test_that("observed means reconstruct coefficients of a one-factor lm", {

  data <- simulate_identity_data()
  stats::contrasts(data$group) <- cbind(c(-3, 1, 2, 0) / 6, c(0, -1, 1, 0) / 2, c(1, 1, 1, -3))

  model <- stats::lm(y ~ group, data = data)
  hypotheses <- .contrast_hypotheses(stats::contrasts(data$group), levels(data$group))

  expect_true(all(.check_hypothesis_identity(
    hypotheses,
    focal_coefficients(model, data, "group"),
    tapply(data$y, data$group, mean)
  )))
})


test_that("polynomial contrasts with unequally spaced scores are recognized", {

  scored <- stats::contr.poly(4, scores = c(0, 10, 20, 40))

  expect_identical(.classify_contrast(scored)$type, "Polynomial")

  # Location and scale of the scores do not matter.
  expect_identical(
    .classify_contrast(stats::contr.poly(4, scores = c(1, 2, 3, 5) * 7 + 2))$type,
    "Polynomial"
  )

  # A matrix whose first column is not followed by its orthogonal
  # polynomials is not polynomial.
  altered <- scored
  altered[, 3L] <- rev(altered[, 3L])

  expect_identical(.classify_contrast(altered)$type, "Custom")
})
