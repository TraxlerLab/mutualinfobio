test_that("Loading and calculating from FASTA alignment ", {
  #print("Calculating MI scores")
  input_file <- "to_test.fasta"
  ali <- Alignment(alignment_path = input_file)
  ali_scores<-getMICovarianceScore(ali, cores = 1)
  ali@mi_data <- ali_scores
  #print mi scores
  print("Printing MI scores")
  print(ali@mi_data)
  # Entropy
  print("Printing Entropy Vector")
  entropy_vec <- getColumnEntropies(ali)
  print(entropy_vec)

  # NMI
  print("Printing Entropy Normalized MI Score")
  nmi_table <- getNormalizedMI(ali)
  print(nmi_table)
  #print("Checking if positions are correct in matrix. This means alignment loaded, and calculated scores")
  expect_equal(ali_scores$i[18], 4)
  expect_equal(ali_scores$j[18], 3)
})
