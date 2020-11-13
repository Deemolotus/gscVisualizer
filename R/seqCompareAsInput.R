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
#' seqCompareAsInput(ID = "hsa-let-7b",
#'                   ref = "CGGGGUGAGGUAGUAGGUUGUGUGGUUUCAGGGCAGUGAUGU
#'                          UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",
#'                   type = "noncoding",
#'                   "hsa-let-7b_UCEC",
#'                   "CGGGGUGAGGUAGUAGGUUGUAUGGUUUCAGGGCAGUGAUGU
#'                    UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",
#'                   "hsa-let-7b_LIRI",
#'                   "CGGGGUGAGGUAGUAGGUUGUGUGGUCUCAGGGCAGUGAUGU
#'                    UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG")
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
#' @importFrom seqinr s2c translate
#' @import NameNeedle

seqCompareAsInput <- function(ID, ref, type, ...){
  #you need to make sure all sequences are coding genes,
  #or all sequence are non-coding genes, not both.

  #check user input
  if (typeof(ID) != "character" || typeof(ref) != "character" ||
      typeof(type) != "character") {
    stop("Invalid input, ID, ref, and type should all be character")
  }

  seqList <- list()
  x <- data.frame(id = ID, sequences = ref, difference = 0)
  seqList <- x
  input <- list(...)
  inputList <- unlist(list(...))

  #check if user input is empty
  if (is.null(inputList)) {
    stop("no other input find")
  }

  #check if ... is valid
  for (item in seq_along(input)) {
    if (typeof(input[[item]]) != "character") {
      stop("Invalid input")
    }
  }

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


########################################################################################


#' Calculate the difference between each pair of gene sequences
#'
#' A function that calculates the difference between gene sequences, takes
#' sequences as function arguments, compare the difference in pairs of sequences.
#' the function returns a list shows all differences. For coding genes,
#' the difference are calculated based on the amino acid sequence.
#' For noncoding gene, the difference are calculated based on Needleman–Wunsch algorithm.
#'
#' @param type A string indicating either the comparison is for coding gene,
#'        or noncoding gene. type = "coding" or type = "noncoding".
#' @param ... A list of sequence string, comparing the difference between each pair of sequence.
#'        For example, first input compare with second sequence,
#'        third one compare to the forth one and so on.
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
#' seqCompareInPair("noncoding", "seq1", "ATCGTCGCC",
#'                  "seq2", "TTCGTCGCG", "seq3", "ATCTTGGCC", "seq4", "ATGGTAGGG")
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
#' @importFrom seqinr s2c translate
#' @import NameNeedle
#'
seqCompareInPair <- function(type, ...){

  #check user input
  if (typeof(type) != "character") {
    stop("Invalid input, ID, ref, and type should all be character")
  }

  seqList <- list()
  input <- list(...)
  inputList <- unlist(list(...))

  #check if user input is empty
  if (is.null(inputList)) {
    stop("no other input find")
  }

  #check if ... is valid
  for (item in seq_along(input)) {
    if (typeof(input[[item]]) != "character") {
      stop("Invalid input")
    }
  }

  n <- length(inputList)
  i <- 1
  dif <- 0
  while (i < n) {

    #calculate the difference between the reference sequence and mutated sequence
    #two situations, the first one is the input genes are coding for proteins,
    #the second one is the input genes are non-coding sequences
    if (type == "coding"){
      #calculate difference
      dif <- cCompareDif(inputList[i + 1], inputList[i + 3])
    } else if (type == "noncoding") {
      #calculate difference
      dif <- nCompareDif(inputList[i + 1], inputList[i + 3])
    } else {
      stop("There is no such type, please try again!")
    }

    #add this comparison into the seqList list
    x <- data.frame(id = inputList[i], sequences = inputList[i + 1], difference = 0)
    y <- data.frame(id = inputList[i + 2], sequences = inputList[i + 3], difference = dif)
    seqList <- rbind(seqList, x)
    seqList <- rbind(seqList, y)
    i <- i + 4
    rm(y)
  }
  return(seqList)
}

############################################################################################


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


############################################################################################


#' Calculate the difference between gene sequences
#'
#' A function that calculates the difference between gene sequences, takes
#' an .fa file as function arguments, the default reference is the first sequence in the .fa file
#' the function returns a list shows all differences. For coding genes,
#' the difference are calculated based on the amino acid sequence.
#' For noncoding gene, the difference are calculated based on Needleman–Wunsch algorithm.
#'
#' @param seqType The sequence ID for refernce sequence.
#' @param type A string indicating either the comparison is for coding gene,
#'        or noncoding gene. type = "coding" or type = "noncoding".
#' @param filePath A string shows the path for the specific .fa file
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
#' # Using a .fa file for miRNA, use this file to
#' # Calculate the difference
#' filePath <- system.file("extdata", "test.fa", package = "gscVisualizer")
#' seqCompareAsFile("RNA", "noncoding", filePath)
#'
#' # Examples 2:
#' # Using a .fa file for DNA sequence, use this file to
#' # Calculate the difference
#' filePath <- system.file("extdata", "seqCompareAsFileTest.fa", package = "gscVisualizer")
#' seqCompareAsFile("DNA", "coding", filePath)
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
#'Pagès H, Aboyoun P, Gentleman R, DebRoy S (2020). Biostrings: Efficient manipulation of biological strings.
#'R package version 2.58.0, \href{https://bioconductor.org/packages/Biostrings}{Link}.
#'
#' @export
#' @importFrom seqinr s2c translate
#' @importFrom  Biostrings readRNAStringSet
#' @import NameNeedle

seqCompareAsFile <- function(seqType, type, filePath){
  #This function takes a .fa file as input and compare every sequence
  #with the reference sequence (reference sequence is the first sequence
  #in the .fa file)

  #check user input
  if (typeof(seqType) != "character" || typeof(type) != "character" ||
      typeof(filePath) != "character") {
    stop("Invalid input, seqType, type and filePath should all be character")
  }

  dif <- 0
  #dataFormater function convert the .fa file to dataframes
  sequenceFA <- dataFormater(seqType, filePath)
  #loop through the sequenceFA dataframes,
  #and then write the difference in the dataframe
  for (geneSeq in seq_along(sequenceFA$id)) {
    #frist set the difference=0 for the first sequence in the dataframe.
    if (geneSeq == 1) {
      sequenceFA[["difference"]] <- 0
    } else {
      #after set the difference for the first sequence,
      #we can use this to compare with other sequences in the dataframes
      if (type == "coding") {
        dif <- cCompareDif(sequenceFA$sequences[1], sequenceFA$sequences[geneSeq])
      } else if (type == "noncoding") {
        dif <- nCompareDif(sequenceFA$sequences[1], sequenceFA$sequences[geneSeq])
      } else {
        stop("the input type is invalid, please enter either coding or noncoding")
      }
      sequenceFA$difference[geneSeq] <- dif
    }
  }
  return(sequenceFA)

}

if (FALSE){
  seqCompareAsFile("RNA", "noncoding", "inst/extdata/test.fa")
  seqCompareAsFile("DNA", "coding", "inst/extdata/seqCompareAsFileTest.fa")
}


##################################################################################


#' Calculate the difference between gene sequences
#'
#' A function that calculates the difference between gene sequences, takes
#' an .fa file as function arguments, the default reference is the first sequence in the .fa file
#' the function returns a list shows all differences. For coding genes,
#' the difference are calculated based on the amino acid sequence.
#' For noncoding gene, the difference are calculated based on Needleman–Wunsch algorithm.
#'
#' @param seqType The sequence ID for refernce sequence.
#' @param type A string indicating either the comparison is for coding gene,
#'        or noncoding gene. type = "coding" or type = "noncoding".
#' @param filePath A string shows the path for the specific .fa file
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
#' # Using a .fa file for miRNA, use this file to
#' # Calculate the difference
#' filePath <- system.file("extdata", "test.fa", package = "gscVisualizer")
#' seqCompareAsFilePair("RNA", "noncoding", filePath)
#'
#' # Examples 2:
#' # Using a .fa file for DNA sequence, use this file to
#' # Calculate the difference
#' filePath <- system.file("extdata", "seqCompareAsFileTest.fa", package = "gscVisualizer")
#' seqCompareAsFilePair("DNA", "coding", filePath)
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
#'Pagès H, Aboyoun P, Gentleman R, DebRoy S (2020). Biostrings: Efficient manipulation of biological strings.
#'R package version 2.58.0, \href{https://bioconductor.org/packages/Biostrings}{Link}.
#'
#' @export
#' @importFrom seqinr s2c translate
#' @importFrom  Biostrings readRNAStringSet
#' @import NameNeedle

seqCompareAsFilePair <- function(seqType, type, filePath){
  #This function takes a .fa file as input and compare the
  #difference between every pair of sequences

  #check user input
  if (typeof(seqType) != "character" || typeof(type) != "character" ||
      typeof(filePath) != "character") {
    stop("Invalid input, seqType, type and filePath should all be character")
  }

  dif <- 0
  sequenceFA <- dataFormater(seqType, filePath)
  sequenceFA[["difference"]] <- 0
  for (geneSeq in seq_along(sequenceFA$id)) {
    #make sure each reference sequence has a difference 0
    if (geneSeq %% 2 != 0) {
      sequenceFA$difference[geneSeq] <- 0
    } else {
      #make comparison between each pair
      if (type == "coding") {
        dif <- cCompareDif(sequenceFA$sequences[geneSeq - 1], sequenceFA$sequences[geneSeq])
      } else if (type == "noncoding") {
        dif <- nCompareDif(sequenceFA$sequences[geneSeq - 1], sequenceFA$sequences[geneSeq])
      } else {
        stop("the input type is invalid, please enter either coding or noncoding")
      }
      sequenceFA$difference[geneSeq] <- dif
    }
  }
  return(sequenceFA)

}

if (FALSE){
  seqCompareAsFilePair("RNA", "noncoding", "inst/extdata/test.fa")
  seqCompareAsFilePair("DNA", "coding", "inst/extdata/seqCompareAsFileTest.fa")
}


########################################################################################


#This function is use to convert the .fa file to dataframes.
#This can make the future analysis easier.
dataFormater <- function(seqType, filePath){
  myFA <- list()
  #check the sequence type in order to use the correct function,
  #if no such type, stop and return an error message.
  if (seqType == "RNA") {
    myFA <- Biostrings::readRNAStringSet(filePath)
  } else if (seqType == "DNA") {
    myFA <- Biostrings::readDNAStringSet(filePath)
  } else {
    stop("Please choose a valid seqType, either DNA or RNA. :(")
  }
  id = names(myFA)
  sequences = paste(myFA)
  sequenceFA <- data.frame(id, sequences)
  return(sequenceFA)
}

# [END]

