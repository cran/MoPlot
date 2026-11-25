#' Reverse Helmert Contrasts
#'
#' Creates a Reverse Helmert contrast matrix by reversing both the rows and
#' columns of a standard Helmert contrast matrix.
#'
#' @param n A factor, whose levels are used; a vector of level names; or a
#'   whole number specifying the number of factor levels (at least 2).
#' @param contrasts Logical indicating whether to return contrasts (default `TRUE`).
#'
#' @return A numeric matrix with `n` rows and `n - 1` columns containing
#'   Reverse Helmert contrasts.
#'
#' @details
#' Reverse Helmert contrasts compare each factor level, except the last,
#' with the mean of the subsequent levels.
#'
#' The contrast matrix is constructed from [stats::contr.helmert()] by
#' reversing the row order and subsequently the column order. This preserves
#' the Helmert coding structure while reversing the direction in which the
#' sequential comparisons are defined.
#'
#' @examples
#' reverse_helmert(3)
#' reverse_helmert(c("Control", "Trt1", "Trt2"))
#' reverse_helmert(factor(c("b", "a", "b", "c")))
#'
#' dat <- data.frame(
#'   group = factor(rep(c("A", "B", "C"), each = 10)),
#'   y = rnorm(30)
#' )
#'
#' contrasts(dat$group) <- reverse_helmert(3)
#' model <- lm(y ~ group, data = dat)
#'
#' @seealso [stats::contr.helmert()]
#'
#' @export
reverse_helmert <- function(n, contrasts = TRUE) {
  # Handle input to accept a factor, a vector of levels, or a numeric scalar.
  # A factor contributes its levels, not its values.
  if (is.factor(n)) {
    n <- levels(n)
  }

  if (is.numeric(n) && length(n) == 1L) {
    if (is.na(n) || n != round(n)) {
      stop("'n' must be a whole number of levels.", call. = FALSE)
    }
    if (n < 2L) {
      stop("not enough degrees of freedom to define contrasts")
    }
    levs <- seq_len(n)
  } else {
    levs <- n
    n <- length(levs)
    if (n < 2L) {
      stop("not enough degrees of freedom to define contrasts")
    }
  }

  # Obtain standard Helmert matrix
  H <- stats::contr.helmert(levs, contrasts = contrasts)

  if (!contrasts) {
    return(H)
  }

  # Reverse rows and columns using direct matrix indexing
  res <- H[n:1L, ncol(H):1L, drop = FALSE]

  # Preserve row names (factor levels)
  rownames(res) <- levs

  res
}
