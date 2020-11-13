context("seqCompareAsFile calculation")
library(gscVisualizer)

test_that("running correctly", {

  seqInfo <- seqCompareAsFile("RNA", "noncoding", "../../inst/extdata/test.fa")

  expect_type(seqInfo, "list")
  expect_s3_class(seqInfo, "data.frame")
  expect_length(seqInfo, 3)
  expect_length(seqInfo$id, 14)
  expect_identical(seqInfo$difference[1], 0)
  expect_identical(seqInfo$difference[2], -56)
  expect_identical(seqInfo$difference[3], -55)
  expect_identical(seqInfo$difference[4], -54)
})

context("Checking for invalid user input for seqCompareAsFile")
test_that("seqCompareAsFile error upon invalid user input", {

  #no such file
  expect_error(seqInfo <- seqCompareAsFile(
    "RNA", "noncoding", "nofilefind"))

  #seqType is invalid
  expect_error(seqInfo <- seqCompareAsFile(
    "tRNA", "noncoding", "../../inst/extdata/test.fa"))

  #type is invalid
  expect_error(seqInfo <- seqCompareAsFile(
    "RNA", "oding", "../../inst/extdata/test.fa"))

})
