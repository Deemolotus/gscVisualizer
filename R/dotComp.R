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
  possiableChar <- c(".", "(", ")")
  if (nchar(seq1) != nchar(seq2)) {
    stop("two sequenec should have the same length")
  }
  seq1Check <- checkSeq(seq1)
  seq2Check <- checkSeq(seq2)
  if (nchar(seq1Check) != 0 || nchar(seq2Check) != 0) {
    stop("dataErr0r, the number of ( not equal to number of )")
  }

  seq1Split <- strsplit(seq1, "")[[1]]
  seq2Split <- strsplit(seq2, "")[[1]]
  for (i in 1:length(seq1Split)){
    if (! seq1Split[i] %in% possiableChar || ! seq2Split[i] %in% possiableChar) {
      stop("wrong type of input, please use ?dotComp to see and example")
    }
    if(seq1Split[i] != seq2Split[i]){
      count <- count + 1
    }
  }
  return(count)
}

#' validate the dot-bracket form sequences
#'
#' A function that check if the sequence is validate, if the number of ( is not equal to
#' the number of ), then the motif sequence is not complete. When calling the dotComp
#' function will also run this funtion for validation.
#'
#' @param seqs A string represent an RNA secondary structure in dot bracket form
#'
#' @return Returns an empty string when the sequence is validate, return an string "error"
#'         if the sequence is not validate
#'
#' @examples
#' #this should return an empty string
#' checkSeq("...(((...)))...")
#'
#' @export
checkSeq <- function(seqs){
  lbCount <- 0
  rbCount <- 0
  seqSplit <- strsplit(seqs, "")[[1]]
  for (i in 1:length(seqSplit)){
    if (seqSplit[i] == "(") {
      lbCount <- lbCount + 1
    } else if (seqSplit[i] == ")") {
      rbCount <- rbCount + 1
    }
  }
  if (rbCount != lbCount) {
    return("error")
  }
  return("")
}

