test_that("hurdle_flex with threshold=0 produces valid output", {
  set.seed(42)
  n <- 500
  x <- rnorm(n)
  y <- rpois(n, exp(0.5 + 0.3 * x))
  df <- data.frame(y = y, x = x)

  m <- hurdle_flex(y ~ x, data = df, threshold = 0)
  expect_s3_class(m, "hurdle_flex")
  expect_s3_class(m, "ecoflex")
  expect_true(!is.null(m$coefficients))
  expect_true(is.finite(m$logLik))
  expect_true(is.finite(m$AIC))
  expect_equal(m$threshold, 0)
})

test_that("hurdle_flex with threshold > 0 works", {
  set.seed(42)
  n <- 1000
  x <- rnorm(n)
  y <- ifelse(runif(n) < 0.4, sample(0:5, n, replace = TRUE),
              rpois(n, exp(1 + 0.5 * x)) + 5)
  df <- data.frame(y = y, x = x)
  m <- hurdle_flex(y ~ x, data = df, threshold = 5)

  expect_s3_class(m, "hurdle_flex")
  expect_equal(m$threshold, 5)
  expect_true(!is.null(m$coefficients))
})

test_that("hurdle_flex negbin works", {
  set.seed(123)
  df <- data.frame(y = rnbinom(500, mu = 5, size = 2), x = rnorm(500))
  m <- hurdle_flex(y ~ x, data = df, threshold = 3, dist = "negbin")

  expect_true("log_alpha" %in% names(coef(m)))
  expect_true(is.finite(m$AIC))
})

test_that("multi-part formula works in hurdle_flex", {
  set.seed(42)
  df <- data.frame(y = rpois(500, 3), x1 = rnorm(500), x2 = rnorm(500), z = rnorm(500))
  m <- hurdle_flex(y ~ x1 + x2 | x1 + z, data = df, threshold = 2)
  expect_true(length(coef(m)) > 0)
})

test_that("summary and print work for hurdle_flex", {
  set.seed(42)
  df <- data.frame(y = rpois(200, 3), x = rnorm(200))
  m <- hurdle_flex(y ~ x, data = df)
  s <- summary(m)
  expect_s3_class(s, "summary.ecoflex")
  expect_true(!is.null(s$coef_table))
  expect_output(print(s))
})

test_that("hurdle_flex two-step works", {
  set.seed(42)
  df <- data.frame(y = rpois(300, 2), x = rnorm(300))
  m <- hurdle_flex(y ~ x, data = df, method = "two-step")
  expect_s3_class(m, "hurdle_flex")
  expect_true(is.finite(m$logLik))
})

test_that("coef, logLik, nobs work", {
  set.seed(42)
  df <- data.frame(y = rpois(200, 3), x = rnorm(200))
  m <- hurdle_flex(y ~ x, data = df)
  expect_true(length(coef(m)) > 0)
  expect_true(is.finite(logLik(m)))
  expect_equal(nobs(m), 200)
})

test_that("predict works for hurdle_flex", {
  set.seed(42)
  df <- data.frame(y = rpois(200, 3), x = rnorm(200))
  m <- hurdle_flex(y ~ x, data = df)
  p <- predict(m)
  expect_true(length(p) > 0)
})
