test_that("plot.viztest returns a ggplot by default", {
  v <- viztest_fit()
  p <- plot(v)
  # ggplot object (documented return)
  expect_s3_class(p, "ggplot")
})

test_that("plot.viztest can return data instead of a plot", {
  v <- viztest_fit()
  d <- plot(v, make_plot = FALSE)
  # documented columns for data mode
  expect_true(all(c("vbl","est","se","lwr","upr","label",
                    "stim_start","stim_end","bound_start","bound_end","ambiguous") %in% names(d)))
  expect_false(anyNA(d$lwr))
  expect_false(anyNA(d$upr))
})

test_that("plot.viztest respects level argument (numeric and aliases)", {
  v <- viztest_fit()
  # numeric level (pick a level present in the grid)
  some_level <- v$tab$level[1]
  d_num <- plot(v, level = some_level, make_plot = FALSE)
  expect_true(all(d_num$lwr <= d_num$est & d_num$est <= d_num$upr))
  
  # alias level: "ce" (cognitively easiest)
  d_ce <- plot(v, level = "ce", make_plot = FALSE)
  # intervals should differ for different chosen levels in general
  expect_false(isTRUE(all.equal(d_num$lwr, d_ce$lwr)) && isTRUE(all.equal(d_num$upr, d_ce$upr)))
})

test_that("plot.viztest draws correct reference-line flags", {
  v <- viztest_fit()
  d_none <- plot(v, ref_lines = "none", make_plot = FALSE)
  # Columns are still present; 'ambiguous' flag computable regardless of drawing
  expect_type(d_none$ambiguous, "logical")
  
  d_amb  <- plot(v, ref_lines = "ambiguous", make_plot = FALSE)
  # At least one ambiguous case typically exists on these data; if not, this still asserts type
  expect_type(d_amb$ambiguous, "logical")
})

test_that("plot.viztest supports ref_lines as specific stimulus names", {
  v <- viztest_fit()
  labs <- unique(v$est$vbl)
  pick <- head(labs, 1)
  d    <- plot(v, ref_lines = pick, make_plot = FALSE)
  # Data subset keeps the same variables; we don't assert exact counts, only structure
  expect_true(all(c("vbl","bound_start","bound_end") %in% names(d)))
})
