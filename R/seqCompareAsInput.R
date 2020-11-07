seqCompareAsInput <- function(ID, ref, type, ...){
  #you need to make sure all sequences are coding genes,
  #or all sequence are non-coding genes, not both.
  seqList <- list()
  x <- data.frame(id = ID, sequences = ref, difference = 0)
  seqList <- x
  inputList <- unlist(list(...))
  n <- length(inputList)
  i <- 1
  dif <- 0
  while (i < n) {
    #calculate the difference between the reference sequence and mutated sequence
    #two situations, the first one is the input genes are coding for proteins,
    #the second one is the input genes are non-coding sequences
    if (type == "coding"){
      dif <- cCompareDif(ref, inputList[i + 1])
    } else if (type == "noncoding") {
      dif <- nCompareDif(ref, inputList[i + 1])
    } else {
      stop("There is no such type, please try again!")
    }
    #add this comparison into the seqList list
    y <- data.frame(id = inputList[i], sequences = inputList[i + 1], difference = dif)
    seqList <- rbind(seqList, y)
    i <- i + 2
    rm(y)
  }
  return(seqList)
}

#comparing poly-peptide chain, for coding gene comparison
cCompareDif <- function(ref, mut){
  dif <- 0
  convert <- seqinr::s2c(ref)
  normalAAs <- seqinr::translate(convert)
  mutated <- seqinr::s2c(mut)
  mutatedAAs <- seqinr::translate(mutated)
  nEnd <- checkEnd(normalAAs)
  mEnd <- checkEnd(mutatedAAs)
  difference <- mEnd-nEnd
  dif <- difference
  comp <- min(mEnd, nEnd)
  j <- 1
  while (j < comp + 1) {
    if (normalAAs[j] != mutatedAAs[j]) {
      if (difference < 0) { #normal longer than mutated
        dif <- dif - 1
      } else { #equal length or mutated longer than normal
        dif <- dif + 1
      }
      j <- j + 1
    } else {
      j <- j + 1
    }
  }
  return(dif)
}

checkEnd <- function(ref){
  stopPosition <- 0
  i <- 1
  while (i < length(ref)) {
    if (ref[i] == "*"){
      stopPosition <- i
      i <- length(ref) + 1
    } else {
      i <- i + 1
    }
  }
  if (stopPosition == 0) {
    stopPosition <- length(ref)
  }
  return(stopPosition)
}

#string alignment method to find the difference, for non-coding gene comparison
nCompareDif <- function(ref, mut){
  rule <- NameNeedle::defaultNeedleParams
  rule$MATCH <- 0
  dif <- NameNeedle::needles(ref, mut, params = rule)$score
  return(dif)
}

seqCompareAsInput(ID = "1", ref = "ATCTAAAAAAAA", type = "noncoding", "2", "ATC")











