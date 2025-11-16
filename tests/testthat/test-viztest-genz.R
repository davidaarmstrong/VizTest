test_that("gen_z returns appropriate values", {
  ests <- viztest_chick_estimates()
  z <- gen_z(coef(ests), vcov(ests))
  expect_length(z, 2)
  expect_named(z, c("ave_z", "all_z"))
  expect_equal(nrow(z$ave_z), 6)
  expect_equal(nrow(z$all_z), 15)
})
