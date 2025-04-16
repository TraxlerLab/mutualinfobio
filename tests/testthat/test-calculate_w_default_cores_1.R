test_that("Loading and calculating from FASTA alignment without specifying number of cores (should default to 1) ", {
  #print("Calculating MI scores")
  input_file <- "to_test.fasta"
  ali <- Alignment(alignment_path = input_file)
  ali_scores<-getMICovarianceScore(ali)
  #print("Expecting change")
  expect_equal(ali_scores$i[18], 4)
  expect_equal(ali_scores$j[18], 3)

})
