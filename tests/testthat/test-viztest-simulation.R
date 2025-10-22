test_that("viztest simulation method works as expected", {

mod <- viztest_fixture_model()
b <- coef(mod)
v <- vcov(mod)
sim <- mvtnorm::rmvnorm(5000, b, v)
v_sim <- make_vt_data(sim, type="sim")
v5 <- VizTest:::viztest.vtsim(v_sim, include_intercept=FALSE, include_zero=FALSE, add_test_level=TRUE)
expect_true(!any(c("(Intercept)", "zero") %in% names(v5$est$vbl)))  
expect_true(all(c("lwr_add", "upr_add") %in% names(v5$est)))

v6 <- viztest(v_sim, include_intercept=FALSE, include_zero=TRUE, add_test_level=TRUE, cifun = "hdi")
expect_true(!any(c("(Intercept)", "zero") %in% names(v5$est$vbl)))  
expect_true(all(c("lwr_add", "upr_add") %in% names(v5$est)))


})
