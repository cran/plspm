test_that("check_manifest_scaling works for right scaling", {
  aux =  c(3, 4, 5, 1, 2, 5, 1)
  MV = iris[,aux]
  good_scaling = list(c("num", "num", "nom"),
                      c("raw", "raw"),
                      c("nom", "ord"))
  
  
#  expect_true(check_manifest_scaling(MV, good_scaling))
  expect_false(all(sapply(c("raw", "raw"),get_metric)))
  expect_false(all(sapply(c("num", "raw"),get_metric)))
})

test_that("check_manifest_scaling throws errors", {  
  aux =  c(3, 4, 5, 1, 2, 5, 1)
  MV = iris[,aux]
  bad_scaling = list(c("num", "num", "ord"),
                     c("raw", "raw"),
                     c("ord", "ord"))
  
#  expect_error(check_manifest_scaling(MV, bad_scaling))
})
