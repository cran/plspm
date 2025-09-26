test_that("check_scheme detects centroid", {
  expect_identical(check_scheme("centroid"), "centroid")
  expect_identical(check_scheme("CENTROID"), "centroid")
  expect_identical(check_scheme("cen"), "centroid")
  expect_identical(check_scheme("CEN"), "centroid")
  expect_identical(check_scheme("c"), "centroid")
  expect_identical(check_scheme("C"), "centroid")
})


test_that("check_scheme detects factorial", {
  expect_identical(check_scheme("factorial"), "factorial")
  expect_identical(check_scheme("FACTORIAL"), "factorial")
  expect_identical(check_scheme("factor"), "factorial")
  expect_identical(check_scheme("FACTOR"), "factorial")
  expect_identical(check_scheme("f"), "factorial")
  expect_identical(check_scheme("F"), "factorial")
})

test_that("check_scheme detects path", {
  expect_identical(check_scheme("path"), "path")
  expect_identical(check_scheme("PATH"), "path")
  expect_identical(check_scheme("p"), "path")
  expect_identical(check_scheme("P"), "path")
})

test_that("check_scheme detects bad schemes", {
  expect_warning(check_scheme(1), 
                 "Invalid 'scheme'. Default 'scheme=centroid' is used.")
  expect_warning(check_scheme("inner"), 
                 "Invalid 'scheme'. Default 'scheme=centroid' is used.")
  expect_warning(check_scheme("centroide"), 
                 "Invalid 'scheme'. Default 'scheme=centroid' is used.")
})
