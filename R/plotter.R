#' Plotting the differences
#'
#' A function that takes in a list of dataframes and convert the
#' numerical number to plot form for easier interpretation,
#' only use this function when there is a reference sequence,
#' if the function used is comparison in pair,
#' plotting them in this way is meaning less
#'
#' @param dataframes A list that contains all id, sequences, differences.
#'
#' @return Returns a histogram.
#'
#' @examples
#' # Examples 1:
#' # plot the difference
#' a <- list()
#' x <- data.frame(id = c("1","2","3","4","5","6","7","8"),
#'                 sequences = c("a", "b", "c", "d", "e", "f", "g", "h"),
#'                 difference = c(1,5,3,4,5,6,6,32))
#' a <- x
#' plotter(a)
#'
#' # Examples 2:
#' # connect this function to seqCompareAsInput.R
#'
#' library(gscVisualizer)
#' filePath <- system.file("extdata", "test.fa", package = "gscVisualizer")
#' a <- seqCompareAsFile("RNA", "noncoding", filePath)
#' plotter(a)
#'
#' @export
#' @importFrom graphics barplot

plotter <- function(dataframes){
  difPlot <- list()

  #check if the dataframes contains the right column
  if(! "id" %in% colnames(dataframes) || ! "sequences" %in% colnames(dataframes)){
    stop("missing column id or sequences")
  }

  if ("difference" %in% colnames(dataframes)) {
    arrangeData <- dataframes$difference[dataframes$difference != 0]
    difGraph <- table(arrangeData)
    difPlot <- graphics::barplot(difGraph, ylab = "Occurances", xlab = "Differences",
                    main = "difference frequence plot")
  } else {
    stop("there is some thing wrong with this dataframe")
  }
  return(difPlot)
}

# [END]
