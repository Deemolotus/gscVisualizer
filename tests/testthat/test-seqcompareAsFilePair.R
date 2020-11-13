context("seqCompareAsFilePair calculation")
library(gscVisualizer)

test_that("running correctly", {

  filePath <- system.file("extdata", "test.fa", package = "gscVisualizer")
  seqInfo <- seqCompareAsFilePair("RNA", "noncoding", filePath)

  expect_type(seqInfo, "list")
  expect_s3_class(seqInfo, "data.frame")
  expect_length(seqInfo, 3)
  expect_length(seqInfo$id, 14)
  expect_identical(seqInfo$difference[1], 0)
  expect_identical(seqInfo$difference[2], -56)
  expect_identical(seqInfo$difference[3], -0)
  expect_identical(seqInfo$difference[4], -56)
})

context("Checking for invalid user input for seqCompareAsFilePair")
test_that("seqCompareAsFilePair error upon invalid user input", {

  filePath <- system.file("extdata", "test.fa", package = "gscVisualizer")

  #no such file
  expect_error(seqInfo <- seqCompareAsFilePair(
    "RNA", "noncoding", "nofilefind"))

  #seqType is invalid
  expect_error(seqInfo <- seqCompareAsFilePair(
    "tRNA", "noncoding", filePath))

  #type is invalid
  expect_error(seqInfo <- seqCompareAsFilePair(
    "RNA", "oding", filePath))

})
