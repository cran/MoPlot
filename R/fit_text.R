# ============================================================
# TEXT LABELS THAT SHRINK TO AVOID OVERLAPS
# ============================================================
#
# The hierarchical group labels and the tree formulas of moplot_multpred()
# have a fixed size in points, whereas the distances between blocks, rows,
# and tree branches depend on the size of the plot. Their size is therefore
# fitted when the plot is drawn, and again whenever the device is resized.


#' Common scale that keeps neighbouring labels apart
#'
#' Labels in the same row (same `y`) are checked pairwise in order of `x`,
#' taking their horizontal justification into account. Vertically, a label
#' with a band (`ymin`, `ymax`) must stay inside it; labels without a band
#' sit on their baseline and each row must fit below the next one. The
#' scale is the largest value not above one at which every constraint holds
#' with at least `padding` to spare, bounded below by `min_scale`. A single
#' scale is returned, so that all labels keep one text size.
#'
#' @param x Anchor positions, in inches.
#' @param y Anchor positions, in inches; labels with the same value form a
#'   row.
#' @param widths,heights Label widths and heights at full size, in inches.
#' @param hjust,vjust Justification of each label.
#' @param ymin,ymax Optional vertical band of each label, in inches; `NA`
#'   for labels without one.
#' @param padding Minimum gap between neighbouring labels, in inches.
#' @param min_scale Smallest scale returned.
#'
#' @return A single number between `min_scale` and one.
#'
#' @keywords internal
#' @noRd
.fit_label_scale <- function(x,
                             y,
                             widths,
                             heights = rep(0, length(x)),
                             hjust = rep(0.5, length(x)),
                             vjust = rep(0, length(x)),
                             ymin = rep(NA_real_, length(x)),
                             ymax = rep(NA_real_, length(x)),
                             padding = 0.06,
                             min_scale = 0.6) {

  scale <- 1

  has_band <- !is.na(ymin) & !is.na(ymax)

  # Vertical, labels with a band: the part below the anchor must stay above
  # ymin and the part above it below ymax.
  below <- has_band & vjust > 0
  above <- has_band & vjust < 1

  band_scales <- c(
    (y[below] - ymin[below] - padding / 2) / (vjust[below] * heights[below]),
    (ymax[above] - y[above] - padding / 2) / ((1 - vjust[above]) * heights[above])
  )

  scale <- min(scale, band_scales[is.finite(band_scales)])

  # Vertical, labels without a band: the tallest label of a row must fit
  # below the next row.
  rows <- sort(unique(y[!has_band]))

  if (length(rows) > 1L) {

    row_heights <- vapply(rows, function(row) max(heights[!has_band & y == row]), numeric(1))

    row_scales <- (diff(rows) - padding / 2) / utils::head(row_heights, -1L)

    scale <- min(scale, row_scales[is.finite(row_scales)])
  }

  # Horizontal: neighbours within a row must not touch. A label extends
  # hjust * width to the left of its anchor and (1 - hjust) * width to the
  # right.
  for (row in unique(y)) {

    in_row <- which(y == row)

    if (length(in_row) < 2L) {
      next
    }

    in_row <- in_row[order(x[in_row])]

    left <- utils::head(in_row, -1L)
    right <- utils::tail(in_row, -1L)

    pair_scales <- (x[right] - x[left] - padding) /
      ((1 - hjust[left]) * widths[left] + hjust[right] * widths[right])

    scale <- min(scale, pair_scales[is.finite(pair_scales) & pair_scales > 0])
  }

  max(min_scale, min(1, scale))
}


#' Geom for labels fitted at drawing time
#'
#' Behaves as `ggplot2::GeomText`, including its default aesthetics and
#' `parse`, but draws a grob whose text size is fitted by
#' `.fit_label_scale()` when it is rendered. The optional aesthetics `ymin`
#' and `ymax` give the vertical band available to each label.
#'
#' @keywords internal
#' @noRd
GeomMoplotFitText <- ggplot2::ggproto(
  "GeomMoplotFitText",
  ggplot2::GeomText,

  optional_aes = c("ymin", "ymax"),

  draw_panel = function(data,
                        panel_params,
                        coord,
                        parse = FALSE,
                        padding = 0.06,
                        min_scale = 0.6) {

    data <- coord$transform(data, panel_params)

    grid::gTree(
      label_data = data,
      parse = parse,
      padding = padding,
      min_scale = min_scale,
      cl = "moplot_fit_text"
    )
  }
)


#' Fit labels to the space available on the device
#'
#' @param x A `"moplot_fit_text"` grob.
#'
#' @return The grob with its text children.
#'
#' @exportS3Method grid::makeContent
#' @noRd
makeContent.moplot_fit_text <- function(x) {

  data <- x$label_data

  number_labels <- nrow(data)

  column_or <- function(name, default) {
    if (is.null(data[[name]])) rep(default, number_labels) else data[[name]]
  }

  labels <- if (isTRUE(x$parse)) {
    as.expression(lapply(as.character(data$label), function(text) parse(text = text, keep.source = FALSE)[[1L]]))
  } else {
    as.character(data$label)
  }

  font_sizes <- data$size * ggplot2::.pt

  font_faces <- column_or("fontface", 1)

  hjust <- column_or("hjust", 0.5)

  vjust <- column_or("vjust", 0.5)

  full_size_labels <- lapply(
    seq_len(number_labels),
    function(row) {
      grid::textGrob(
        labels[row],
        gp = grid::gpar(fontsize = font_sizes[row], fontface = font_faces[row])
      )
    }
  )

  label_widths <- vapply(
    full_size_labels,
    function(label) grid::convertWidth(grid::grobWidth(label), "inches", valueOnly = TRUE),
    numeric(1)
  )

  label_heights <- vapply(
    full_size_labels,
    function(label) grid::convertHeight(grid::grobHeight(label), "inches", valueOnly = TRUE),
    numeric(1)
  )

  to_inches <- function(values) {
    inches <- rep(NA_real_, number_labels)
    known <- !is.na(values)
    if (any(known)) {
      inches[known] <- grid::convertY(grid::unit(values[known], "npc"), "inches", valueOnly = TRUE)
    }
    inches
  }

  scale <- .fit_label_scale(
    x = grid::convertX(grid::unit(data$x, "npc"), "inches", valueOnly = TRUE),
    y = grid::convertY(grid::unit(data$y, "npc"), "inches", valueOnly = TRUE),
    widths = label_widths,
    heights = label_heights,
    hjust = hjust,
    vjust = vjust,
    ymin = to_inches(column_or("ymin", NA_real_)),
    ymax = to_inches(column_or("ymax", NA_real_)),
    padding = x$padding,
    min_scale = x$min_scale
  )

  text <- grid::textGrob(
    label = labels,
    x = grid::unit(data$x, "npc"),
    y = grid::unit(data$y, "npc"),
    hjust = hjust,
    vjust = vjust,
    gp = grid::gpar(
      col = ggplot2::alpha(data$colour, data$alpha),
      fontsize = font_sizes * scale,
      fontface = font_faces,
      fontfamily = column_or("family", "")
    )
  )

  grid::setChildren(x, grid::gList(text))
}
