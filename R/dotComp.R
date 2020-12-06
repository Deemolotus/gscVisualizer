#' Calculate the difference between RNA secondary dot bracket form sequences
#'
#' A function that calculates the difference between two dot bracket form of RNA sequences,
#' takes two sequences as function arguments.
#'
#' @param seq1 A string represent an RNA secondary structure in dot bracket form
#' @param seq2 A string represent an RNA secondary structure in dot bracket form
#'
#' @return Returns a number represent the different in two structure
#'
#' @examples
#' dotComp("...(((...)))...", ".(.(.(...)))...")
#'
#'
#' @export

dotComp <- function(seq1, seq2){
  count <- 0
  if (nchar(seq1) != nchar(seq2)) {
    stop("two sequenec should have the same length")
  }
  seq1Split <- strsplit(seq1, "")[[1]]
  seq2Split <- strsplit(seq2, "")[[1]]
  for (i in 1:length(seq1Split)){
    if(seq1Split[i] != seq2Split[i]){
      count <- count + 1
    }
  }
  return(count)
}
