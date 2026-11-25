# Text is measured with the Cairo device used by moplot_export(); pdf()
# cannot encode symbols such as beta and warns about each of them.
open_test_device <- function(width, height) {
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(tempfile(fileext = ".pdf"), width = width, height = height)
  } else {
    grDevices::pdf(NULL, width = width, height = height)
  }
}


test_that(".fit_label_scale() shrinks labels only when neighbours would touch", {

  # Two labels 1.5 inches wide whose centres are 2 inches apart fit.
  expect_identical(.fit_label_scale(x = c(0, 2), y = c(1, 1), widths = c(1.5, 1.5)), 1)

  # One inch apart they need (1 - 0.06) / 1.5 of their size.
  expect_equal(.fit_label_scale(x = c(0, 1), y = c(1, 1), widths = c(1.5, 1.5)), 0.94 / 1.5)

  # Labels in different rows do not constrain each other horizontally.
  expect_identical(.fit_label_scale(x = c(0, 0.1), y = c(1, 2), widths = c(1.5, 1.5)), 1)

  # A row must fit below the next one: 0.25 inches tall, 0.2 inches apart.
  expect_equal(
    .fit_label_scale(x = c(0, 0), y = c(1, 1.2), widths = c(1, 1), heights = c(0.25, 0.25)),
    (0.2 - 0.03) / 0.25
  )

  # The scale never drops below min_scale.
  expect_identical(.fit_label_scale(x = c(0, 0.1), y = c(1, 1), widths = c(2, 2)), 0.6)
})


test_that("group labels are fitted to the size of the device", {

  labels <- grid::gTree(
    label_data = data.frame(
      x = c(0.30, 0.55),
      y = c(0.9, 0.9),
      label = c("Delayed", "Mean(Delayed,Immediate)"),
      size = 3.87,
      fontface = "bold",
      colour = "black",
      alpha = NA,
      vjust = 0
    ),
    padding = 0.06,
    min_scale = 0.6,
    cl = "moplot_fit_text"
  )

  drawn_font_size <- function(width) {
    open_test_device(width = width, height = 4)
    on.exit(grDevices::dev.off())
    grid::pushViewport(grid::viewport())
    fitted <- grid::makeContent(labels)
    unique(fitted$children[[1L]]$gp$fontsize)
  }

  full_size <- 3.87 * ggplot2::.pt

  expect_equal(drawn_font_size(12), full_size)
  expect_lt(drawn_font_size(5), full_size)
  expect_gte(drawn_font_size(5), 0.6 * full_size)
})


test_that("moplot_multpred() draws its group labels with the fitted geom", {

  set.seed(6)
  data <- expand.grid(
    encoding = factor(c("Shallow", "Deep", "Self-ref")),
    delay = factor(c("Immediate", "Delayed")),
    id = 1:10
  )
  data$recall <- stats::rnorm(nrow(data))
  stats::contrasts(data$encoding) <- stats::contr.helmert(3)
  stats::contrasts(data$delay) <- stats::contr.sum(2)

  model <- stats::lm(recall ~ encoding * delay, data = data)

  for (style in c("minimal", "tree")) {

    plot <- suppressMessages(suppressWarnings(
      moplot_multpred(model, data = data, style = style, caption = FALSE, interaction_info = FALSE)
    ))

    geoms <- vapply(plot[[1L]]$layers, function(layer) class(layer$geom)[1L], character(1))

    expect_true("GeomMoplotFitText" %in% geoms, info = style)

    # Drawing on a small device does not fail.
    open_test_device(width = 5, height = 4)
    expect_no_error(print(plot))
    grDevices::dev.off()
  }
})


test_that(".fit_label_scale() accounts for justification and vertical bands", {

  # Left-justified labels extend only to the right of their anchor.
  expect_equal(
    .fit_label_scale(x = c(0, 1), y = c(1, 1), widths = c(1.5, 1.5), hjust = c(0, 0)),
    0.94 / 1.5
  )

  # Labels extending away from each other never collide.
  expect_identical(
    .fit_label_scale(x = c(0, 1), y = c(1, 1), widths = c(1.5, 1.5), hjust = c(1, 0)),
    1
  )

  # A label hanging below its anchor (vjust = 1) must stay above ymin.
  expect_equal(
    .fit_label_scale(x = 0, y = 1, widths = 1, heights = 0.2, vjust = 1, ymin = 0.8, ymax = 1),
    (0.2 - 0.03) / 0.2
  )

  # A centred label needs half its height on each side.
  expect_equal(
    .fit_label_scale(x = 0, y = 1, widths = 1, heights = 0.4, vjust = 0.5, ymin = 0.8, ymax = 1.5),
    (0.2 - 0.03) / 0.2
  )
})


fitted_scales <- function(plot, width, height) {
  open_test_device(width = width, height = height)
  on.exit(grDevices::dev.off())
  print(plot)
  forced <- grid::grid.force()
  scales <- list()
  for (name in unique(grid::grid.ls(print = FALSE)$name)) {
    grob <- tryCatch(grid::grid.get(name), error = function(e) NULL)
    if (inherits(grob, "moplot_fit_text") && length(grob$children) > 0L) {
      parsed <- isTRUE(grob$parse)
      scales[[if (parsed) "formulas" else "labels"]] <- unique(
        grob$children[[1L]]$gp$fontsize / (grob$label_data$size * ggplot2::.pt)
      )[1L]
    }
  }
  scales
}

tree_data <- function(ways) {
  set.seed(4)
  factors <- stats::setNames(
    lapply(LETTERS[seq_len(ways)], function(name) factor(paste0(tolower(name), 1:2))),
    LETTERS[seq_len(ways)]
  )
  data <- do.call(expand.grid, c(factors, list(replicate = 1:12)))
  data$y <- stats::rnorm(nrow(data), 0.5 * (data$A == "a2"))
  for (variable in LETTERS[seq_len(ways)]) {
    stats::contrasts(data[[variable]]) <- stats::contr.sum(2)
  }
  data
}


test_that("tree formulas sit inside their bands and fit at the recommended size", {

  for (ways in 3:4) {

    data <- tree_data(ways)
    formula <- stats::as.formula(paste("y ~", paste(LETTERS[seq_len(ways)], collapse = " * ")))

    plot <- suppressMessages(suppressWarnings(
      moplot_multpred(stats::lm(formula, data = data), data = data, style = "tree",
                      caption = FALSE, interaction_info = FALSE)
    ))

    formula_layer <- Filter(
      function(layer) inherits(layer$geom, "GeomMoplotFitText") && "Band_low" %in% names(layer$data),
      plot[[1L]]$layers
    )[[1L]]

    nodes <- formula_layer$data

    expect_true(all(nodes$Band_low <= nodes$Formula_y + 1e-12), info = ways)
    expect_true(all(nodes$Formula_y <= nodes$Band_high + 1e-12), info = ways)

    size <- attr(plot, "moplot_export_size")

    scales <- fitted_scales(plot, size[["width"]], size[["height"]])

    expect_equal(scales$formulas, 1, info = ways)
    expect_equal(scales$labels, 1, info = ways)

    # On a small device the formulas shrink, but not below 60%.
    small <- fitted_scales(plot, 6, 5)

    expect_lt(small$formulas, 1)
    expect_gte(small$formulas, 0.6)
  }
})
