context("plotter to plot all differences")
library(gscVisualizer)

test_that("running correctly", {

  a <- list()
  x <- data.frame(id = c("1","2","3","4","5","6","7","8"),
                 sequences = c("a", "b", "c", "d", "e", "f", "g", "h"),
                 difference = c(1,5,3,4,5,6,6,32))
  a <- x
  seqInfo <- plotter(a)

  expect_type(seqInfo, "double")
  expect_length(seqInfo, 5)
})

test_that("running correctly", {

  filePath <- system.file("extdata", "test.fa", package = "gscVisualizer")
  a <- seqCompareAsFile("RNA", "noncoding", filePath)
  seqInfo <- plotter(a)

  expect_type(seqInfo, "double")
  expect_length(seqInfo, 8)
})


context("Checking for invalid user input for plotter")
test_that("plotter error upon invalid user input", {

  #no column difference
  a <- list()
  x <- data.frame(id = c("1","2","3","4","5","6","7","8"),
                  sequences = c("a", "b", "c", "d", "e", "f", "g", "h"),
                  somthing = c(1,5,3,4,5,6,6,32))
  a <- x
  expect_error(seqInfo <- plotter(a))

  #missing column
  a <- list()
  x <- data.frame(id = c("1","2","3","4","5","6","7","8"))
  a <- x
  expect_error(seqInfo <- plotter(a))

  #no right column at all
  a <- list()
  x <- data.frame(b = c("1","2","3","4","5","6","7","8"),
                  v = c("a", "b", "c", "d", "e", "f", "g", "h"),
                  somthing = c(1,5,3,4,5,6,6,32))
  a <- x
  expect_error(seqInfo <- plotter(a))


})
