test_that("viztest returns a well-formed viztest object", {
  v <- viztest_fit()
  
  expect_s3_class(v, "viztest")
  # expected top-level components per manual
  expect_true(all(c("tab","pw_test","ci_tests","combs","param_names","L","U","est") %in% names(v)))
  
  # tab has the documented columns
  expect_true(all(c("level","psame","pdiff","easy") %in% names(v$tab)))
  
  # est has the documented columns
  expect_true(all(c("vbl","est","se") %in% names(v$est)))
  
  # basic sanity: levels grid covers a range and psame/pdiff/easy are numeric
  expect_gt(nrow(v$tab), 1)
  expect_type(v$tab$psame, "double")
  expect_type(v$tab$pdiff, "double")
  expect_type(v$tab$easy,  "double")
  
  v2 <- viztest_fit(add_test_level =TRUE)
  expect_true(all(c("vbl", "est", "se", "lwr_add", "upr_add") %in% names(v2$est)))
  
  v3 <- viztest_fit(add_test_level = TRUE, adjust = "bonferroni")
  expect_equal(round(attr(v3, "test_level"), 2), .05)
  expect_true(attr(v3, "adjust") == "bonferroni")
  expect_snapshot_output(print(v3, level =.9))
  
  
  mod <- viztest_fixture_model()
  b <- coef(mod)
  v <- vcov(mod)
  v_md <- make_vt_data(b, v, type="est_var")
  v4 <- viztest(v_md, add_test_level = TRUE, include_zero=TRUE, include_intercept=TRUE)
  expect_true(all(c("(Intercept)", "zero") %in% v4$est$vbl))
  
  z <- gen_z(b, v)
  expect_snapshot_output(print(z$ave_z))
  
})
