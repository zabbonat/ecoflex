test_that("hurdle_flex returns correct class and structure", {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  y <- rpois(n, exp(0.5 + 0.3 * x))
  df <- data.frame(y = y, x = x)

  m <- hurdle_flex(y ~ x, data = df, threshold = 0)

  expect_s3_class(m, "hurdle_flex")
  expect_s3_class(m, "ecoflex")
  expect_true(is.numeric(coef(m)))
  expect_true(length(coef(m)) > 0)
  expect_true(is.matrix(vcov(m)))
  expect_equal(nobs(m), n)
})

test_that("hurdle_flex summary works", {
  set.seed(1)
  df <- data.frame(y = rpois(150, 2), x = rnorm(150))
  m  <- hurdle_flex(y ~ x, data = df)
  s  <- summary(m)
  expect_true(!is.null(s))
})

test_that("hurdle_flex predict returns numeric vector", {
  set.seed(2)
  df <- data.frame(y = rpois(150, 2), x = rnorm(150))
  m  <- hurdle_flex(y ~ x, data = df)
  p  <- predict(m, type = "response")
  expect_true(is.numeric(p))
  expect_equal(length(p), 150)
})

test_that("hurdle_flex threshold argument is validated", {
  df <- data.frame(y = rpois(50, 2), x = rnorm(50))
  expect_error(hurdle_flex(y ~ x, data = df, threshold = -1))
  expect_error(hurdle_flex(y ~ x, data = df, threshold = 0.5))
})

test_that("hurdle_flex works with negbin distribution", {
  set.seed(3)
  df <- data.frame(y = rpois(200, 3), x = rnorm(200))
  m  <- hurdle_flex(y ~ x, data = df, dist = "negbin")
  expect_s3_class(m, "hurdle_flex")
})

test_that("hurdle_flex two-step estimation works", {
  set.seed(4)
  df <- data.frame(y = rpois(200, 2), x = rnorm(200))
  m  <- hurdle_flex(y ~ x, data = df, method = "two-step")
  expect_s3_class(m, "hurdle_flex")
  expect_true(is.numeric(coef(m)))
})
