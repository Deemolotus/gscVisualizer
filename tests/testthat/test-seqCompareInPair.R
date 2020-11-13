context("seqCompareInPair calculation")
library(gscVisualizer)

test_that("running correctly", {

  seqInfo <- seqCompareInPair(
    "noncoding", "seq1", "ATCGTCGCC", "seq2", "TTCGTCGCG",
    "seq3", "ATCTTGGCC", "seq4", "ATGGTAGGG")

  expect_type(seqInfo, "list")
  expect_s3_class(seqInfo, "data.frame")
  expect_length(seqInfo, 3)
  expect_identical(seqInfo$difference[1], 0)
  expect_identical(seqInfo$difference[2], -2)
  expect_identical(seqInfo$difference[3], -0)
  expect_identical(seqInfo$difference[4], -5)
})

context("Checking for invalid user input for seqCompareInPair")
test_that("seqCompareInPair error upon invalid user input", {

  #type is invalid
  expect_error(seqInfo <- seqCompareInPair(
    "something", "seq1", "ATCGTCGCC", "seq2", "TTCGTCGCG",
    "seq3", "ATCTTGGCC", "seq4", "ATGGTAGGG"))

  #no item in list
  expect_error(seqInfo <- seqCompareInPair(
    "noncoding"))

  #item in list is invalid
  expect_error(seqInfo <- seqCompareInPair(
    "noncoding", 1, "ATCGTCGCC", "seq2", "TTCGTCGCG"))
})


