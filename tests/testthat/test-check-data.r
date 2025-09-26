test_that("check_data works as expected with matrices", {
  some_matrix = matrix(1:10, 5, 2)
  good_matrix = some_matrix
  rownames(good_matrix) = 1:5
  colnames(good_matrix) = c("MV1", "MV2")

  expect_identical(check_data(some_matrix), good_matrix)
})

test_that("check_data works as expected with data frames", {
  good_df = iris[c(1:5,51:55,101:105),]
  
  expect_identical(check_data(good_df), good_df)
})

test_that("check_data detects bad data", {
  bad_data = list(1:10)
  bad_matrix = matrix(letters[1:15], 5, 3)
  one_row = iris[1,]
  one_column = as.matrix(1:5)
  
  expect_error(check_data(1:10))
  expect_error(check_data("string"))
  expect_error(check_data(bad_data))
  expect_error(check_data(bad_matrix))
  expect_error(check_data(one_row))
  expect_error(check_data(one_column))
})
