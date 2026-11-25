test_that("model_contrast_table() detects Reverse Helmert contrasts", {

  set.seed(123)

  dat <- data.frame(
    group = factor(rep(c("A", "B", "C", "D"), each = 20)),
    y = rnorm(80)
  )

  contrasts(dat$group) <- reverse_helmert(4)

  mod <- stats::lm(y ~ group, data = dat)

  out <- model_contrast_table(mod)

  model_summary <- out[["Model-Level Summary"]]
  variable_summary <- out[["Variable-Level Summary"]]

  expect_equal(
    model_summary$Value[
      model_summary$Metric == "Contrast Matrix Types"
    ],
    "Reverse Helmert"
  )

  expect_equal(
    variable_summary$Contrast_Matrix_Type[
      variable_summary$Variable == "group"
    ],
    "Reverse Helmert"
  )
})

test_that("model_contrast_table() recognizes Treatment coding with a non-first reference level", {

  dat <- data.frame(
    group = factor(
      rep(c("A", "B", "C"), each = 10),
      levels = c("A", "B", "C")
    ),
    y = rnorm(30)
  )

  contrasts(dat$group) <- stats::contr.treatment(
    3,
    base = 2
  )

  mod <- stats::lm(
    y ~ group,
    data = dat
  )

  out <- model_contrast_table(
    mod,
    type = "all"
  )

  variable_summary <- out[["Variable-Level Summary"]]

  expect_equal(
    variable_summary$Contrast_Matrix_Type[
      variable_summary$Variable == "group"
    ],
    "Treatment"
  )
})
