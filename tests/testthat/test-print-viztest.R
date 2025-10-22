test_that("print.viztest default summary is stable", {
  v <- viztest_fit()
  expect_snapshot_output(print(v))   # stores text in _snaps/
})

test_that("print.viztest respects arguments (best / missed_tests / level)", {
  v <- viztest_fit()
  
  # only best level(s)
  expect_snapshot_output(print(v, best = TRUE,  missed_tests = FALSE))
  
  # show all levels
  expect_snapshot_output(print(v, best = FALSE, missed_tests = FALSE))
  
  # force a particular level (numeric or alias)
  lev <- v$tab$level[which.max(v$tab$easy)]  # pick one deterministically
  expect_snapshot_output(print(v, level = lev))
})
