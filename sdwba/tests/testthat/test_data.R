context("Data sets")

test_that("All data sets exist", {
  expect_named(generic.krill.Conti2006)
  expect_is(generic.krill.Conti2006, "data.frame")
  expect_named(generic.krill.McGeehee1998)
  expect_is(generic.krill.McGeehee1998, "data.frame")
})
