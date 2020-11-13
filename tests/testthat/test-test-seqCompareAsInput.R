context("seqCompareAsInput calculation")
library(gscVisualizer)

test_that("function running correctly", {

  seqInfo <- seqCompareAsInput(ID = "hsa-let-7b",
                               ref = "CGGGGUGAGGUAGUAGGUUGUGUGGUUUCAGGGCAGUGAUGU
                               UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",
                               type = "noncoding",
                               "hsa-let-7b_UCEC",
                               "CGGGGUGAGGUAGUAGGUUGUAUGGUUUCAGGGCAGUGAUGU
                                UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",
                               "hsa-let-7b_LIRI",
                               "CGGGGUGAGGUAGUAGGUUGUGUGGUCUCAGGGCAGUGAUGU
                                UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG")

  expect_type(seqInfo, "list")
  expect_s3_class(seqInfo, "data.frame")
  expect_length(seqInfo, 3)
  expect_identical(seqInfo$difference[1], 0)
  expect_identical(seqInfo$difference[2], -2)
  expect_identical(seqInfo$difference[3], -2)
})

context("Checking for invalid user input for seqCompareAsInput")
test_that("seqCompareAsInput error upon invalid user input", {

  #ID is not a string
  expect_error(seqInfo <- seqCompareAsInput(
    ID = 2,
    ref = "CGGGGUGAGGUAGUAGGUUGUGUGGUUUCAGGGCAGUGAUGU
           UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",
    type = "noncoding",
    "hsa-let-7b_UCEC",
    "CGGGGUGAGGUAGUAGGUUGUAUGGUUUCAGGGCAGUGAUGU
     UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",
    "hsa-let-7b_LIRI",
    "CGGGGUGAGGUAGUAGGUUGUGUGGUCUCAGGGCAGUGAUGU
     UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG"))

  #ref is not correct type
  expect_error(seqInfo <- seqCompareAsInput(
    ID = "hsa-let-7b",
    ref = 12,
    type = "noncoding",
    "hsa-let-7b_UCEC",
    "CGGGGUGAGGUAGUAGGUUGUAUGGUUUCAGGGCAGUGAUGU
     UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",
    "hsa-let-7b_LIRI",
    "CGGGGUGAGGUAGUAGGUUGUGUGGUCUCAGGGCAGUGAUGU
     UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG"))

  #type is invalid
  expect_error(seqInfo <- seqCompareAsInput(
    ID = "hsa-let-7b",
    ref = "CGGGGUGAGGUAGUAGGUUGUGUGGUUUCAGGGCAGUGAUGU
           UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",
    type = "something",
    "hsa-let-7b_UCEC",
    "CGGGGUGAGGUAGUAGGUUGUAUGGUUUCAGGGCAGUGAUGU
     UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",
    "hsa-let-7b_LIRI",
    "CGGGGUGAGGUAGUAGGUUGUGUGGUCUCAGGGCAGUGAUGU
     UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG"))

  #following list is invalid
  expect_error(seqInfo <- seqCompareAsInput(
    ID = "hsa-let-7b",
    ref = "CGGGGUGAGGUAGUAGGUUGUGUGGUUUCAGGGCAGUGAUGU
           UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",
    type = "something",
    "hsa-let-7b_UCEC",
    "CGGGGUGAGGUAGUAGGUUGUAUGGUUUCAGGGCAGUGAUGU
     UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",
    1,
    "CGGGGUGAGGUAGUAGGUUGUGUGGUCUCAGGGCAGUGAUGU
     UGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG"))

})
