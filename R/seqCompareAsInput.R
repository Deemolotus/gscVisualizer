#' Calculate the difference between gene sequences
#'
#' A function that calculates the difference between gene sequences, takes
#' sequences as function arguments, every sequences are compare to the reference sequence.
#' the function returns a list shows all differences. For coding genes,
#' the difference are calculated based on the amino acid sequence.
#' For noncoding gene, the difference are calculated based on Needleman–Wunsch algorithm.
#'
#' @param ID The sequence ID for refernce sequence.
#' @param ref A string for the reference sequence.
#' @param type A string indicating either the comparison is for coding gene,
#'        or noncoding gene. type = "coding" or type = "noncoding".
#' @param ... A list of mutated sequence string where the odd index indicate the mutated sequence ID,
#'        the even index indicate the mutated sequence.
#'
#' @return Returns an List object of class seqCompareAsInput with results.
#' \itemize{
#'   \item id - A string name for the correspond sequence.
#'   \item sequences - A string indicate the sequence.
#'   \item difference - A value of class "numeric" indicating the
#'         difference between this sequence and reference sequence.
#' }
#'
#' @examples
#' # Examples 1:
#' # Using random sequences for noncoding type
#' # Calculate the difference
#' seqCompareAsInput(ID = "hsa-let-7b", ref = "CGGGGUGAGGUAGUAGGUUGUGUGGUUUCAGGGCAGUGAUGUUGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG", type = "noncoding",
#'                   "hsa-let-7b_UCEC", "CGGGGUGAGGUAGUAGGUUGUAUGGUUUCAGGGCAGUGAUGUUGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",
#'                   "hsa-let-7b_LIRI", "CGGGGUGAGGUAGUAGGUUGUGUGGUCUCAGGGCAGUGAUGUUGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG")
#'
#' # Examples 2:
#' # Using random sequences for coding type
#' # Calculate the difference
#' seqCompareAsInput(ID = "NM_000546.6",
#'                   ref = "CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGC
#'                        TGGGAGCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGGCAGCCAGACTGCCTTCCGGGTCACTG
#'                        CCATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATG
#'                        GAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCC
#'                        CCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTG
#'                        CTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCT
#'                        GTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCT
#'                        GGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGA
#'                        CCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTA
#'                        CAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGAT
#'                        GGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAA
#'                        ACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCA
#'                        CTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACA
#'                        CTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGA
#'                        GAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAG
#'                        CACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAA
#'                        TATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAAC
#'                        TCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAA
#'                        GGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGACATTCT
#'                        CCACTTCTTGTTCCCCACTGACAGCCTCCCACCCCCATCTCTCCCTCCCCTGCCATTTTGGGTTTTGGGT
#'                        CTTTGAACCCTTGCTTGCAATAGGTGTGCGTCAGAAGCACCCAGGACTTCCATTTGCTTTGTCCCGGGGC
#'                        TCCACTGAACAAGTTGGCCTGCACTGGTGTTTTGTTGTGGGGAGGAGGATGGGGAGTAGGACATACCAGC
#'                        TTAGATTTTAAGGTTTTTACTGTGAGGGATGTTTGGGAGATGTAAGAAATGTTCTTGCAGTTAAGGGTTA
#'                        GTTTACAATCAGCCACATTCTAGGTAGGGGCCCACTTCACCGTACTAACCAGGGAAGCTGTCCCTCACTG
#'                        TTGAATTTTCTCTAACTTCAAGGCCCATATCTGTGAAATGCTGGCATTTGCACCTACCTCACAGAGTGCA
#'                        TTGTGAGGGTTAATGAAATAATGTACATCTGGCCTTGAAACCACCTTTTATTACATGGGGTCTAGAACTT
#'                        GACCCCCTTGAGGGTGCTTGTTCCCTCTCCCTGTTGGTCGGTGGGTTGGTAGTTTCTACAGTTGGGCAGC
#'                        TGGTTAGGTAGAGGGAGTTGTCAAGTCTCTGCTGGCCCAGCCAAACCCTGTCTGACAACCTCTTGGTGAA
#'                        CCTTAGTACCTAAAAGGAAATCTCACCCCATCCCACACCCTGGAGGATTTCATCTCTTGTATATGATGAT
#'                        CTGGATCCACCAAGACTTGTTTTATGCTCAGGGTCAATTTCTTTTTTCTTTTTTTTTTTTTTTTTTCTTT
#'                        TTCTTTGAGACTGGGTCTCGCTTTGTTGCCCAGGCTGGAGTGGAGTGGCGTGATCTTGGCTTACTGCAGC
#'                        CTTTGCCTCCCCGGCTCGAGCAGTCCTGCCTCAGCCTCCGGAGTAGCTGGGACCACAGGTTCATGCCACC
#'                        ATGGCCAGCCAACTTTTGCATGTTTTGTAGAGATGGGGTCTCACAGTGTTGCCCAGGCTGGTCTCAAACT
#'                        CCTGGGCTCAGGCGATCCACCTGTCTCAGCCTCCCAGAGTGCTGGGATTACAATTGTGAGCCACCACGTC
#'                        CAGCTGGAAGGGTCAACATCTTTTACATTCTGCAAGCACATCTGCATTTTCACCCCACCCTTCCCCTCCT
#'                        TCTCCCTTTTTATATCCCATTTTTATATCGATCTCTTATTTTACAATAAAACTTTGCTGCCA",
#'                  type = "coding",
#'                  "NM_001126113.2|:1-2651",
#'                  "GATGGGATTGGGGTTTTCCCCTCCCATGTGCTCAAGACTGGCGCTAAAAGTTTTGAGCTTCTCAAAAGTC
#'                  TAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTG
#'                  CTTTCCACGACGGTGACACGCTTCCCTGGATTGGCAGCCAGACTGCCTTCCGGGTCACTGCCATGGAGGA
#'                  GCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTT
#'                  CCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATA
#'                  TTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGT
#'                  GGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCT
#'                  GTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCA
#'                  AGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGT
#'                  GCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCA
#'                  CAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCC
#'                  CTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCG
#'                  ACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTAC
#'                  ATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACT
#'                  CCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCG
#'                  CACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGA
#'                  GCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCC
#'                  TTCAGATGCTACTTGACTTACGATGGTGTTACTTCCTGATAAACTCGTCGTAAGTTGAAAATATTATCCG
#'                  TGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGG
#'                  AAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCC
#'                  ATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGACATTCTCCACTTCTTGTTCCCCACTG
#'                  ACAGCCTCCCACCCCCATCTCTCCCTCCCCTGCCATTTTGGGTTTTGGGTCTTTGAACCCTTGCTTGCAA
#'                  TAGGTGTGCGTCAGAAGCACCCAGGACTTCCATTTGCTTTGTCCCGGGGCTCCACTGAACAAGTTGGCCT
#'                  GCACTGGTGTTTTGTTGTGGGGAGGAGGATGGGGAGTAGGACATACCAGCTTAGATTTTAAGGTTTTTAC
#'                  TGTGAGGGATGTTTGGGAGATGTAAGAAATGTTCTTGCAGTTAAGGGTTAGTTTACAATCAGCCACATTC
#'                  TAGGTAGGGGCCCACTTCACCGTACTAACCAGGGAAGCTGTCCCTCACTGTTGAATTTTCTCTAACTTCA
#'                  AGGCCCATATCTGTGAAATGCTGGCATTTGCACCTACCTCACAGAGTGCATTGTGAGGGTTAATGAAATA
#'                  ATGTACATCTGGCCTTGAAACCACCTTTTATTACATGGGGTCTAGAACTTGACCCCCTTGAGGGTGCTTG
#'                  TTCCCTCTCCCTGTTGGTCGGTGGGTTGGTAGTTTCTACAGTTGGGCAGCTGGTTAGGTAGAGGGAGTTG
#'                  TCAAGTCTCTGCTGGCCCAGCCAAACCCTGTCTGACAACCTCTTGGTGAACCTTAGTACCTAAAAGGAAA
#'                  TCTCACCCCATCCCACACCCTGGAGGATTTCATCTCTTGTATATGATGATCTGGATCCACCAAGACTTGT
#'                  TTTATGCTCAGGGTCAATTTCTTTTTTCTTTTTTTTTTTTTTTTTTCTTTTTCTTTGAGACTGGGTCTCG
#'                  CTTTGTTGCCCAGGCTGGAGTGGAGTGGCGTGATCTTGGCTTACTGCAGCCTTTGCCTCCCCGGCTCGAG
#'                  CAGTCCTGCCTCAGCCTCCGGAGTAGCTGGGACCACAGGTTCATGCCACCATGGCCAGCCAACTTTTGCA
#'                  TGTTTTGTAGAGATGGGGTCTCACAGTGTTGCCCAGGCTGGTCTCAAACTCCTGGGCTCAGGCGATCCAC
#'                  CTGTCTCAGCCTCCCAGAGTGCTGGGATTACAATTGTGAGCCACCACGTCCAGCTGGAAGGGTCAACATC
#'                  TTTTACATTCTGCAAGCACATCTGCATTTTCACCCCACCCTTCCCCTCCTTCTCCCTTTTTATATCCCAT
#'                  TTTTATATCGATCTCTTATTTTACAATAAAACTTTGCTGCCACCTGTGTGTCTGAGGGGTG"
#')
#'
#' @references
#'Charif D, Lobry J (2007). “SeqinR 1.0-2: a contributed package to the R project for statistical
#'computing devoted to biological sequences retrieval and analysis.” In Bastolla U, Porto M, Roman H,
#'Vendruscolo M (eds.), \emph{Structural approaches to sequence evolution: Molecules, networks, populations,
#'series Biological and Medical Physics, Biomedical Engineering}, 207-232.
#'Springer Verlag, New York. ISBN : 978-3-540-35305-8.
#'
#'Kevin R. Coombes (2020). NameNeedle: Using Needleman-Wunsch to Match Sample Names. R package version
#'1.2.6/r51. \href{https://R-Forge.R-project.org/projects/nameneedle/}{Link}
#'
#'
#' @export
#' @import seqinr
#' @import NameNeedle

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
      #calculate difference
      dif <- cCompareDif(ref, inputList[i + 1])
    } else if (type == "noncoding") {
      #calculate difference
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
  #use seqinr package to convert the sequence to amino acid sequence
  convert <- seqinr::s2c(ref)
  normalAAs <- seqinr::translate(convert)
  mutated <- seqinr::s2c(mut)
  mutatedAAs <- seqinr::translate(mutated)
  #check where is the stop codon
  nEnd <- checkEnd(normalAAs)
  mEnd <- checkEnd(mutatedAAs)
  #get the difference in length
  difference <- mEnd-nEnd
  dif <- difference
  comp <- min(mEnd, nEnd)
  j <- 1
  #calculate difference
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

#This helper function is use to calcuate the length of the amino acid sequence,
#and return the stop codon position
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
#using Needleman–Wunsch algorithm
nCompareDif <- function(ref, mut){
  #change to rule, make sure match get 0 point
  rule <- NameNeedle::defaultNeedleParams
  rule$MATCH <- 0
  #get the differece
  dif <- NameNeedle::needles(ref, mut, params = rule)$score
  return(dif)
}



