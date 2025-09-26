test_that("check_modes works as expected", {
  modA = c("A", "A", "A")
  modB = c("B", "B", "B")
  modnewA = c("NEWA", "NEWA", "NEWA")
  modAB = c("A", "B", "B")
  modAnewA = c("A", "newA", "A")
  modBnewA = c("newA", "B", "B")  
  modABnewA = c("A", "B", "newA")
  block3 = list(1:2, 3:4, 5:6)
  
  expect_identical(check_modes(NULL, block3), modA)
  expect_identical(check_modes(modA, block3), modA)
  expect_identical(check_modes(modB, block3), modB)
  expect_identical(check_modes(modnewA, block3), modnewA)
  expect_identical(check_modes(modAB, block3), modAB)
  expect_error(check_modes(modAnewA, block3), "Sorry. Can't work with both modes 'A' and 'newA'")
  expect_error(check_modes(modABnewA, block3), "Sorry. Can't work with both modes 'A' and 'newA'")
})

test_that("check_mdoes detects bad modes", {
  modA2 = c("A", "A")
  modAAS = c("A", "A", "S")  
  block3 = list(1:2, 3:4, 5:6)
  
  expect_warning(check_modes(modA2, block3))
  expect_error(check_modes(modAAS, block3), "Sorry. Unrecognized mode: 'S'")
})
