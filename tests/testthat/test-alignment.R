library(testthat)
library(mutualinfobio)

test_that("FASTA sequences are read correctly", {
  aln <- Alignment(alignmentID = "test", alignment_path = system.file("extdata", "example.fasta", package = "mutualinfobio"))
  expect_true(nrow(aln@alignment_data) > 0)
  expect_equal(colnames(aln@alignment_data), c("lineID", "sequence"))
})

test_that("Sequence names and retrieval work", {
  aln <- Alignment("test", system.file("extdata", "example.fasta", package = "mutualinfobio"))
  names <- get_sequence_names(aln)
  expect_type(names, "character")
  expect_true("seq1" %in% names)

  seq <- get_sequence(aln, "seq1")
  expect_true(nchar(seq) > 0)
})

test_that("Alignment length is computed", {
  aln <- Alignment("test", system.file("extdata", "example.fasta", package = "mutualinfobio"))
  len <- get_alignment_length(aln)
  expect_gt(len, 0)
})

test_that("MI CSV can be attached using helper function", {
  aln <- Alignment("test", system.file("extdata", "example.fasta", package = "mutualinfobio"))
  #aln@mi_data <- readMICSV(system.file("extdata", "example_mi.csv", package = "mutualinfobio"))
  aln <- addMICSVData(aln, system.file("extdata", "example_mi.csv", package = "mutualinfobio"))

  expect_s3_class(aln@mi_data, "data.frame")
  expect_true(all(c("i", "j", "score") %in% names(aln@mi_data)))
})

test_that("getTopMIFor returns highest scoring pair for a given i", {
  fasta_path <- system.file("extdata", "example.fasta", package = "mutualinfobio")
  mi_path <- system.file("extdata", "example_mi.csv", package = "mutualinfobio")

  aln <- Alignment("test", fasta_path)
  aln@mi_data <- readMICSV(mi_path)

  top <- getTopMIFor(aln, 0)
  expect_type(top, "list")
  expect_equal(top$j, 5)
  expect_equal(top$score, 0.95)

  top2 <- getTopMIFor(aln, 99)
  expect_null(top2)
})

test_that("grepNames finds matching sequence names", {
  fasta_path <- system.file("extdata", "example.fasta", package = "mutualinfobio")
  aln <- Alignment("test", fasta_path)

  hits <- grepNames(aln, "^seq")
  expect_true(length(hits) == 3)
  expect_true(all(hits %in% c("seq1", "seq2", "seq3")))

  none <- grepNames(aln, "missing")
  expect_length(none, 0)
})


test_that("getScoredResidues uses getTopMIFor and skips gaps", {
  fasta_path <- system.file("extdata", "example.fasta", package = "mutualinfobio")
  mi_path <- system.file("extdata", "example_mi.csv", package = "mutualinfobio")

  aln <- Alignment("test", fasta_path)
  aln@mi_data <- readMICSV(mi_path)

  df <- getScoredResidues(aln, "seq1")
  #print(df)
  expect_s3_class(df, "data.frame")
  expect_named(df, c("residue_index", "alignment_position", "residue", "score"))
  expect_true(all(df$residue != "-"))
  expect_true(all(df$score >= 0))
})

test_that("readMask reads columns correctly", {
  mask_path <- system.file("extdata", "example_mask.txt", package = "mutualinfobio")
  mask <- readMask(mask_path)

  expect_type(mask, "integer")
  expect_true(all(mask >= 0))
  expect_true(length(mask) > 0)
})




