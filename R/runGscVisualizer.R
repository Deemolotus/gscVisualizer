#' Launch Shiny App For Package gscVisualizer
#'
#' A function that launches the shiny app for this package.
#' The shiny app allows user to choose a input file, the app
#' will calculate the difference among those genes and generate
#' a plot shows the distribution for all differences.
#' The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but open up a shiny page.
#'
#' @examples
#' \dontrun{
#' gscVisualizer::runGscVisualizer()
#' }
#'
#' @author Zhiwen Tan, \email{leon.tan@mail.utoronto.ca}
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#' Steipe B., ABC project (.utility 4.07) A Bioinformatics Course: Applied Bioinformatics
#' \href{http://steipe.biochemistry.utoronto.ca/abc/index.php/Bioinformatics_Main_Page}{Link}.
#'
#' @export
#' @importFrom shiny runApp
runGscVisualizer <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "gscVisualizer")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}

# [END]
