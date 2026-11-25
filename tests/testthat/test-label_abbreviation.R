abbreviation_key <- function(labels) {
  key <- .abbreviate_moplot_labels(labels)$key
  if (is.null(key)) return(character(0))
  stats::setNames(key$abbreviated, key$full)
}


test_that("labels sharing a prefix keep the part that distinguishes them", {

  expect_identical(
    abbreviation_key(c("Control", "Treatment1", "Treatment2")),
    c(Control = "Cont.", Treatment1 = "Tre.1", Treatment2 = "Tre.2")
  )

  expect_identical(
    abbreviation_key(c("TreatmentHigh", "TreatmentLow")),
    c(TreatmentHigh = "Tre.High", TreatmentLow = "Tre.Low")
  )
})


test_that("abbreviations are always shorter than the labels", {

  labels <- c("Treatment1", "Treatment2", "Treatment10", "Group_A_long", "Group_B_long", "Placebo")

  key <- abbreviation_key(labels)

  expect_true(all(nchar(key) < nchar(names(key))))
  expect_false(anyDuplicated(key) > 0)
})


test_that("short labels are not abbreviated", {

  expect_length(abbreviation_key(c("alpha", "beta", "A", "B")), 0)
})


test_that("labels inside Mean() and W() groups are abbreviated", {

  result <- .abbreviate_moplot_labels(
    c("Mean(Treatment1, Treatment2)", "W(Control .20, Treatment2 .80)")
  )

  expect_identical(
    unname(result$labels),
    c("Mean(Tre.1,Tre.2)", "W(Cont. .20, Tre.2 .80)")
  )
})
