context("plotter to plot all differences")
library(gscVisualizer)

test_that("running correctly", {

  result <- dotComp("...(((...)))...", ".(.(.(...)))...")

  expect_type(result, "double")
  expect_equal(result, 2)
})

test_that("running correctly", {

  result <- dotComp("....", ".(.).")

  expect_type(result, "double")
  expect_equal(result, 3)
})

context("Checking for invalid user input for plotter")
test_that("dotComp error upon invalid user input", {

  #when input is no string
  expect_error(dotComp(123,456))

  #when input string is invalid
  expect_error(dotComp("asdasdasd", "asdasd"))

  #when input sequence is not complete
  expect_error(dotComp("..(().", "....."))

})
