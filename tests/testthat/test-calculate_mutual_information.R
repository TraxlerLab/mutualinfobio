test_that("Loading and calculating from FASTA alignment ", {
  #print("Calculating MI scores")
  input_file <- "to_test.fasta"
  ali <- Alignment(alignment_path = input_file)
  ali_scores<-getMICovarianceScore(ali, cores = 1)
  #print("Checking if positions are correct in matrix. This means alignment loaded, and calculated scores")
  expect_equal(ali_scores$i[18], 4)
  expect_equal(ali_scores$j[18], 3)
})
