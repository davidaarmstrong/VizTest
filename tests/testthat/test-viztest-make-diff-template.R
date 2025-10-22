test_that("make_diff_template returns well-formed object", {

  ests <- viztest_chick_estimates()
  pred <- ests$estimate
  names(pred) <- ests$feed
  diff_template <- make_diff_template(pred, include_zero=FALSE)
  expect_snapshot_output(print(diff_template))
  sig_diff <- c(0,0,1,1,1, 0,1,1,1,0,0,1,0,1,1)  
  
  v_diff <- viztest(ests, include_zero=FALSE, sig_diffs=sig_diff)
  expect_equal(round(max(v_diff$tab$psame), 2), 0.93)   
})
