viztest_fixture_model <- function() {
  dat <- mtcars
  dat$cyl <- as.factor(dat$cyl)
  dat$hp  <- scale(dat$hp)
  dat$wt  <- scale(dat$wt)
  stats::lm(qsec ~ hp + wt + cyl, data = dat)
}

viztest_fit <- function(...) {
  VizTest::viztest(viztest_fixture_model(), ...)
}

viztest_chick_estimates <- function() {
  data("chickwts", package = "datasets")
  m <- stats::lm(weight ~ feed, data = chickwts)
  marginaleffects::predictions(m, variables="feed", by="feed")
}
