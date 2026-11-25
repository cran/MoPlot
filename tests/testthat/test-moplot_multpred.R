test_that("moplot_multpred() returns a patchwork object for a two-way interaction", {

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2"))
  )

  dat <- dat[
    rep(seq_len(nrow(dat)), each = 20),
  ]

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

  plot <- moplot_multpred(
    model = mod,
    data = dat
  )

  expect_s3_class(
    plot,
    "patchwork"
  )
})


test_that("moplot_multpred() returns only the upper panel when show_coef = FALSE", {

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2"))
  )

  dat <- dat[
    rep(seq_len(nrow(dat)), each = 20),
  ]

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

  plot <- moplot_multpred(
    model = mod,
    data = dat,
    show_coef = FALSE
  )

  expect_s3_class(
    plot,
    "ggplot"
  )
})


test_that("multi-predictor coefficient positions match upper-panel block midpoints", {

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2"))
  )

  dat <- dat[
    rep(seq_len(nrow(dat)), each = 20),
  ]

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

  contrast_output <- contrasts_multpred(
    model = mod,
    data = dat,
    configuration = "all"
  )

  prep <- .prepare_moplot_multpred_data(
    contrast_output
  )

  expected_midpoints <- aggregate(
    x ~ coefficient,
    data = prep$point_data,
    FUN = function(z) mean(range(z))
  )

  observed_midpoints <- prep$coefficient_data[
    ,
    c("coefficient", "x")
  ]

  expected_midpoints <- expected_midpoints[
    match(
      observed_midpoints$coefficient,
      expected_midpoints$coefficient
    ),
    ,
    drop = FALSE
  ]

  expect_equal(
    observed_midpoints$x,
    expected_midpoints$x
  )
})


test_that("moplot_multpred() respects top_split", {

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2"))
  )

  dat <- dat[
    rep(seq_len(nrow(dat)), each = 20),
  ]

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

  contrast_output <- contrasts_multpred(
    model = mod,
    data = dat,
    top_split = "A",
    configuration = "all"
  )

  expect_equal(
    unique(contrast_output$split_var),
    "A"
  )

  plot <- moplot_multpred(
    model = mod,
    data = dat,
    top_split = "A"
  )

  expect_s3_class(
    plot,
    "patchwork"
  )
})

test_that("moplot_multpred() coloured styles assign one colour per contrast", {

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
    rnorm(
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

  out <- contrasts_multpred(
    model = mod,
    data = dat,
    configuration = "all",
    methodological_note = FALSE
  )

  prep <- .prepare_moplot_multpred_data(
    out
  )

  for (current_style in c("tree")) {

    upper <- .plot_moplot_multpred_upper(
      prep = prep,
      style = current_style
    )

    build <- ggplot2::ggplot_build(
      upper
    )

    segment_layer_index <- which(
      vapply(
        build$data,
        function(layer) {

          has_segment_columns <- all(
            c(
              "x",
              "y",
              "xend",
              "yend"
            ) %in% names(layer)
          )

          is_dashed <- has_segment_columns &&
            "linetype" %in% names(layer) &&
            all(
              layer$linetype == "dashed"
            )

          is_dashed
        },
        logical(1)
      )
    )

    expect_length(
      segment_layer_index,
      1L
    )

    segment_layer <- build$data[[segment_layer_index]]

    expect_equal(
      length(
        unique(
          segment_layer$colour
        )
      ),
      length(
        unique(
          prep$segment_data$contrast_id
        )
      )
    )
  }
})

test_that("moplot_multpred() coloured styles keep upper and lower contrast colours aligned", {

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
    rnorm(
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

  out <- contrasts_multpred(
    model = mod,
    data = dat,
    configuration = "all",
    methodological_note = FALSE
  )

  prep <- .prepare_moplot_multpred_data(
    out
  )

  parameters <- .prepare_moplot_multpred_parameters(
    model = mod,
    coefficient_data = prep$coefficient_data
  )

  for (current_style in c("tree")) {

    upper <- .plot_moplot_multpred_upper(
      prep = prep,
      style = current_style
    )

    lower <- .plot_moplot_multpred_lower(
      prep = prep,
      parameters = parameters,
      style = current_style
    )

    upper_build <- ggplot2::ggplot_build(
      upper
    )

    lower_build <- ggplot2::ggplot_build(
      lower
    )

    # Elementary dashed contrast segments in the upper panel.
    upper_segment_index <- which(
      vapply(
        upper_build$data,
        function(layer) {

          has_segment_columns <- all(
            c(
              "x",
              "y",
              "xend",
              "yend"
            ) %in% names(layer)
          )

          is_dashed <- has_segment_columns &&
            "linetype" %in% names(layer) &&
            all(
              layer$linetype == "dashed"
            )

          is_dashed
        },
        logical(1)
      )
    )

    expect_length(
      upper_segment_index,
      1L
    )

    upper_colours <- unique(
      upper_build$data[[upper_segment_index]]$colour
    )

    # X markers representing model coefficients in the lower panel.
    lower_point_indices <- which(
      vapply(
        lower_build$data,
        function(layer) {
          all(
            c(
              "x",
              "y",
              "shape",
              "colour"
            ) %in% names(layer)
          ) &&
            !all(
              c(
                "xmin",
                "xmax"
              ) %in% names(layer)
            )
        },
        logical(1)
      )
    )

    expect_true(
      length(lower_point_indices) >= 1L
    )

    lower_point_index <- lower_point_indices[
      which.max(
        vapply(
          lower_build$data[
            lower_point_indices
          ],
          function(layer) {
            length(
              unique(
                layer$colour
              )
            )
          },
          integer(1)
        )
      )
    ]

    lower_colours <- unique(
      lower_build$data[[lower_point_index]]$colour
    )

    expect_setequal(
      upper_colours,
      lower_colours
    )
  }
})

test_that("moplot_multpred() respects coefficient selection", {

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
    rnorm(nrow(dat), 0, 1)

  contrasts(dat$A) <- reverse_helmert(3)
  contrasts(dat$B) <- MASS::contr.sdif(2)

  mod <- stats::lm(
    y ~ A * B,
    data = dat
  )

  all_out <- contrasts_multpred(
    model = mod,
    data = dat,
    configuration = "all"
  )

  selected_coefficient <- unique(
    all_out$coefficient
  )[1L]

  selected_out <- contrasts_multpred(
    model = mod,
    data = dat,
    coefficients = selected_coefficient,
    configuration = "all"
  )

  expect_equal(
    unique(selected_out$coefficient),
    selected_coefficient
  )

  expect_equal(
    length(unique(selected_out$contrast_id)),
    1L
  )

  plot <- moplot_multpred(
    model = mod,
    data = dat,
    coefficients = selected_coefficient
  )

  expect_s3_class(
    plot,
    "patchwork"
  )
})

test_that("moplot_multpred() respects sig_level", {

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2")),
    rep = seq_len(30)
  )

  dat$y <- with(
    dat,
    10 +
      2 * (A == "A2") +
      4 * (A == "A3") +
      2 * (B == "B2") +
      1.5 * (A == "A2" & B == "B2")
  ) +
    rnorm(nrow(dat), 0, 2)

  contrasts(dat$A) <- reverse_helmert(3)
  contrasts(dat$B) <- MASS::contr.sdif(2)

  mod <- stats::lm(
    y ~ A * B,
    data = dat
  )

  out_05 <- contrasts_multpred(
    model = mod,
    data = dat,
    sig_level = 0.05,
    configuration = "all"
  )

  out_strict <- contrasts_multpred(
    model = mod,
    data = dat,
    sig_level = 1e-10,
    configuration = "all"
  )

  sig_05 <- unique(
    out_05[
      ,
      c("coefficient", "p_value", "significance")
    ]
  )

  sig_strict <- unique(
    out_strict[
      ,
      c("coefficient", "p_value", "significance")
    ]
  )

  expect_equal(
    sig_05$p_value,
    sig_strict$p_value
  )

  expect_true(
    sum(sig_strict$significance == "Sig") <=
      sum(sig_05$significance == "Sig")
  )
})

test_that("moplot_multpred() propagates confidence level", {

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2")),
    rep = seq_len(20)
  )

  dat$y <- rnorm(
    nrow(dat),
    10,
    2
  )

  contrasts(dat$A) <- reverse_helmert(3)
  contrasts(dat$B) <- MASS::contr.sdif(2)

  mod <- stats::lm(
    y ~ A * B,
    data = dat
  )

  out_95 <- contrasts_multpred(
    model = mod,
    data = dat,
    conf_level = 0.95,
    configuration = "all"
  )

  out_80 <- contrasts_multpred(
    model = mod,
    data = dat,
    conf_level = 0.80,
    configuration = "all"
  )

  width_95 <- out_95$ci_high - out_95$ci_low
  width_80 <- out_80$ci_high - out_80$ci_low

  expect_true(
    all(
      width_80 <= width_95,
      na.rm = TRUE
    )
  )
})

test_that("moplot_multpred() supports interpreted coefficient scale", {

  skip_if_not_installed(
    "lme4"
  )

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2")),
    subject = factor(seq_len(40))
  )

  eta <- with(
    dat,
    -1 +
      0.5 * (A == "A2") +
      1.0 * (A == "A3") +
      0.4 * (B == "B2") +
      0.8 * (A == "A2" & B == "B2")
  )

  dat$y <- stats::rbinom(
    nrow(dat),
    1,
    stats::plogis(eta)
  )

  contrasts(dat$A) <- reverse_helmert(3)
  contrasts(dat$B) <- MASS::contr.sdif(2)

  mod <- lme4::glmer(
    y ~ A * B + (1 | subject),
    data = dat,
    family = stats::binomial()
  )

  out <- contrasts_multpred(
    model = mod,
    data = dat,
    configuration = "all",
    estimate_type = "predicted"
  )

  prep <- .prepare_moplot_multpred_data(
    out
  )

  interpreted <- .prepare_moplot_multpred_parameters(
    model = mod,
    coefficient_data = prep$coefficient_data,
    coefficient_scale = "interpreted"
  )

  expect_true(
    all(
      interpreted$parameter_data$Display_type == "Odds ratio"
    )
  )
})

test_that("moplot_multpred() accepts custom panel heights", {

  set.seed(123)

  dat <- expand.grid(
    A = factor(c("A1", "A2", "A3")),
    B = factor(c("B1", "B2")),
    rep = seq_len(20)
  )

  dat$y <- rnorm(
    nrow(dat)
  )

  contrasts(dat$A) <- reverse_helmert(3)
  contrasts(dat$B) <- MASS::contr.sdif(2)

  mod <- stats::lm(
    y ~ A * B,
    data = dat
  )

  plot <- moplot_multpred(
    model = mod,
    data = dat,
    panel_heights = c(3, 1)
  )

  expect_s3_class(
    plot,
    "patchwork"
  )
})

test_that("moplot_multpred() rejects removed delta style", {

  expect_error(
    moplot_multpred(
      model = mod3,
      interaction = "A:B:C",
      style = "delta"
    )
  )
})


test_that("moplot_multpred() handles variable and level names with underscores", {

  set.seed(51)

  dat <- expand.grid(
    time_point = factor(c("t_1", "t_2", "t_3")),
    drug_type = factor(c("d_a", "d_b")),
    site_id = factor(c("s_1", "s_2")),
    replicate = 1:10
  )
  dat$y <- stats::rnorm(nrow(dat), (dat$time_point == "t_3") * (dat$drug_type == "d_b"))

  contrasts(dat$time_point) <- stats::contr.helmert(3)

  two_way <- stats::lm(y ~ time_point * drug_type, data = dat)
  three_way <- stats::lm(y ~ time_point * drug_type * site_id, data = dat)

  for (style in c("minimal", "tree")) {

    plot <- suppressMessages(
      moplot_multpred(two_way, data = dat, style = style, caption = FALSE, interaction_info = FALSE)
    )
    expect_s3_class(plot, "patchwork")

    plot <- suppressMessages(
      moplot_multpred(three_way, data = dat, style = style, caption = FALSE, interaction_info = FALSE)
    )
    expect_s3_class(plot, "patchwork")
  }

  decomposition <- attr(
    suppressMessages(contrasts_multpred(two_way, data = dat, methodological_note = FALSE)),
    "decomposition"
  )

  expect_identical(unique(decomposition$variable_path), "drug_type > time_point")
})


test_that("moplot_multpred() handles transformed covariates and factor outcomes", {

  set.seed(92)

  dat <- expand.grid(A = factor(c("a1", "a2", "a3")), B = factor(c("b1", "b2")), replicate = 1:15)
  dat$x <- stats::rexp(nrow(dat)) + 0.1
  dat$y <- stats::rnorm(nrow(dat), (dat$A == "a3") * (dat$B == "b2") + log(dat$x))
  dat$answer <- factor(ifelse(stats::runif(nrow(dat)) < 0.4, "yes", "no"))

  plot <- suppressMessages(
    moplot_multpred(
      stats::lm(y ~ A * B + log(x), data = dat),
      estimate_type = "predicted",
      caption = FALSE,
      interaction_info = FALSE
    )
  )
  expect_s3_class(plot, "patchwork")

  plot <- suppressMessages(
    moplot_multpred(
      stats::glm(answer ~ A * B, family = stats::binomial, data = dat),
      caption = FALSE,
      interaction_info = FALSE
    )
  )
  expect_s3_class(plot, "patchwork")
})


test_that("the final tree marker keeps clear of the x axis", {

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

  for (formula in list(y ~ A * B, y ~ A * B * C)) {

    model <- stats::lm(formula, data = data)

    plot <- suppressMessages(suppressWarnings(
      moplot_multpred(model, data = data, style = "tree", caption = FALSE, interaction_info = FALSE)
    ))

    upper <- plot[[1L]]

    y_limits <- upper$coordinates$limits$y

    final_layer <- Filter(
      function(layer) inherits(layer$geom, "GeomPoint") && "Final_x_y" %in% names(layer$data),
      upper$layers
    )[[1L]]

    # The marker centre lies at least 7% of the panel height above the axis,
    # whatever the size of the plot.
    clearance <- (min(final_layer$data$Final_x_y) - y_limits[1L]) / diff(y_limits)

    expect_gte(clearance, 0.07 - 1e-8)
  }
})


test_that("four-way trees keep clear of the x axis and label zero exactly", {

  set.seed(4)
  data <- expand.grid(
    A = factor(c("a1", "a2")),
    B = factor(c("b1", "b2")),
    C = factor(c("c1", "c2")),
    D = factor(c("d1", "d2")),
    replicate = 1:12
  )
  data$y <- stats::rnorm(
    nrow(data),
    0.8 * (data$A == "a2") * (data$B == "b2") * (data$C == "c2") * (data$D == "d2") +
      0.5 * (data$A == "a2")
  )
  for (variable in c("A", "B", "C", "D")) {
    stats::contrasts(data[[variable]]) <- stats::contr.sum(2)
  }

  model <- stats::lm(y ~ A * B * C * D, data = data)

  plot <- suppressMessages(suppressWarnings(
    moplot_multpred(model, data = data, style = "tree", caption = FALSE, interaction_info = FALSE)
  ))

  upper <- plot[[1L]]

  y_limits <- upper$coordinates$limits$y

  final_layer <- Filter(
    function(layer) inherits(layer$geom, "GeomPoint") && "Final_x_y" %in% names(layer$data),
    upper$layers
  )[[1L]]

  expect_gte(
    (min(final_layer$data$Final_x_y) - y_limits[1L]) / diff(y_limits),
    0.07 - 1e-8
  )

  # seq() left -5.55e-17 in place of zero, and the axis switched to
  # scientific notation.
  y_breaks <- upper$scales$get_scales("y")$breaks

  expect_true(0 %in% y_breaks)
  expect_false(any(grepl("e", format(y_breaks, trim = TRUE))))
})
