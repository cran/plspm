test_that("get_metric works", {  
  expect_true(get_metric(NULL))
  expect_false(all(sapply(c("ord", "ord"),get_metric)))
  expect_false(all(sapply(c("num", "raw"),get_metric)))
})
