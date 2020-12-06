context("plotter to plot all differences")
library(gscVisualizer)

test_that("running correctly", {

  result <- checkSeq("...(((...)))...")

  expect_type(result, "character")
  expect_equal(nchar(result), 0)
  expect_equal(result, "")
})

context("Checking for invalid user input for plotter")
test_that("checkSeq error upon invalid user input", {

  #when input is not a string
  expect_error(result <- checkSeq(46565465))

  #when input dot-bracket form is invalid
  result <- checkSeq("...(((...)...")

  expect_type(result, "character")
  expect_equal(nchar(result), 5)
  expect_equal(result, "error")

})
