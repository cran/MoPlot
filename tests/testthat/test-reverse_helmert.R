test_that("reverse_helmert() returns the expected contrasts for three levels", {

  expected <- matrix(
    c(
      2, -1, -1,
      0,  1, -1
    ),
    nrow = 3,
    ncol = 2
  )

  rownames(expected) <- as.character(seq_len(3))

  expect_equal(
    reverse_helmert(3),
    expected
  )
})


test_that("reverse_helmert() returns the expected contrasts for four levels", {

  expected <- matrix(
    c(
      3, -1, -1, -1,
      0,  2, -1, -1,
      0,  0,  1, -1
    ),
    nrow = 4,
    ncol = 3
  )

  rownames(expected) <- as.character(seq_len(4))

  expect_equal(
    reverse_helmert(4),
    expected
  )
})


test_that("reverse_helmert() returns a contrast matrix of the expected dimensions", {

  for (n in 2:10) {

    result <- reverse_helmert(n)

    expect_true(is.matrix(result))
    expect_equal(nrow(result), n)
    expect_equal(ncol(result), n - 1L)
  }
})


test_that("reverse_helmert() uses the levels of a factor", {

  group <- factor(c("b", "a", "b", "c", "a", "c"))

  expect_equal(
    reverse_helmert(group),
    reverse_helmert(levels(group))
  )
  expect_identical(rownames(reverse_helmert(group)), c("a", "b", "c"))
})


test_that("reverse_helmert() rejects a non-integer number of levels", {

  expect_error(reverse_helmert(2.5), "whole number")
  expect_error(reverse_helmert(NA_real_), "whole number")
})


test_that("reverse_helmert() contrasts sum to zero", {

  for (n in 2:10) {

    result <- reverse_helmert(n)

    expect_equal(
      colSums(result),
      rep(0, n - 1L)
    )
  }
})
