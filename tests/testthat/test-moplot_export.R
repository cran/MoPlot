export_data <- function(levels = c("A", "B", "C", "D")) {
  set.seed(123)
  data <- data.frame(group = factor(rep(levels, each = 30), levels = levels))
  data$y <- stats::rnorm(nrow(data), as.integer(data$group) / 3)
  data
}

quiet <- function(expression) {
  suppressMessages(suppressWarnings(expression))
}


test_that("moplot_onepred() stores recommended export dimensions", {

  data <- export_data()
  stats::contrasts(data$group) <- stats::contr.helmert(4)
  model <- stats::lm(y ~ group, data = data)

  combined <- quiet(moplot_onepred(model, "group", caption = FALSE))
  upper_only <- quiet(moplot_onepred(model, "group", show_coef = FALSE, caption = FALSE))
  taller_lower <- quiet(moplot_onepred(model, "group", panel_heights = c(1, 1), caption = FALSE))

  size <- attr(combined, "moplot_export_size")

  expect_named(size, c("width", "height"))
  expect_true(all(size > 0))

  # The coefficient panel adds height in proportion to panel_heights.
  expect_equal(attr(upper_only, "moplot_export_size")[["height"]], 4.5)
  expect_equal(size[["height"]], 6.8)
  expect_equal(attr(taller_lower, "moplot_export_size")[["height"]], 9)
  expect_identical(attr(upper_only, "moplot_export_size")[["width"]], size[["width"]])
})


test_that("recommended widths grow with the displayed quantities", {

  width_for <- function(levels, coding) {
    data <- export_data(levels)
    stats::contrasts(data$group) <- coding(length(levels))
    model <- stats::lm(y ~ group, data = data)
    attr(quiet(moplot_onepred(model, "group", caption = FALSE)), "moplot_export_size")[["width"]]
  }

  # Treatment contrasts display two quantities per contrast.
  expect_lt(width_for(LETTERS[1:3], stats::contr.treatment), width_for(LETTERS[1:6], stats::contr.treatment))

  # Never narrower than 6 inches.
  expect_gte(width_for(LETTERS[1:3], stats::contr.helmert), 6)

  # Polynomial plots need one panel per component.
  expect_lt(width_for(LETTERS[1:4], stats::contr.poly), width_for(LETTERS[1:5], stats::contr.poly))
  expect_equal(width_for(LETTERS[1:4], stats::contr.poly), 8.2)
})


test_that("moplot_export() uses the recommended dimensions of moplot_onepred()", {

  data <- export_data()
  stats::contrasts(data$group) <- stats::contr.poly(4)
  plot <- quiet(moplot_onepred(stats::lm(y ~ group, data = data), "group", caption = FALSE))

  file <- tempfile(fileext = ".png")
  on.exit(unlink(file))

  expect_message(
    moplot_export(plot, file, dpi = 50),
    "at 8.2 x 6.8 inches"
  )

  expect_true(file.exists(file))

  # Explicit dimensions override the recommendation.
  expect_message(
    moplot_export(plot, file, width = 5, height = 4, dpi = 50),
    "at 5 x 4 inches"
  )
})


test_that("moplot_export() requires dimensions for plots without a recommendation", {

  plot <- ggplot2::ggplot()

  expect_error(
    moplot_export(plot, tempfile(fileext = ".png")),
    "Recommended export dimensions were not found"
  )
})
