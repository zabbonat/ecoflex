test_that("tobit_flex type 1 normal works", {
  set.seed(42)
  n <- 500; x <- rnorm(n)
  y_star <- 2 + 1.5 * x + rnorm(n)
  y <- pmax(y_star, 0)
  df <- data.frame(y = y, x = x)

  m <- tobit_flex(y ~ x, data = df, tobit_type = 1, left = 0)
  expect_s3_class(m, "tobit_flex")
  expect_s3_class(m, "ecoflex")
  expect_true(m$n_censored_left > 0)
  expect_true(is.finite(m$logLik))
})

test_that("tobit_flex type 1 t-distribution works", {
  set.seed(42)
  n <- 300; x <- rnorm(n)
  y <- pmax(1 + x + rnorm(n), 0)
  df <- data.frame(y = y, x = x)

  m <- tobit_flex(y ~ x, data = df, tobit_type = 1, left = 0, error_dist = "t")
  expect_s3_class(m, "tobit_flex")
  expect_true(is.finite(m$AIC))
})

test_that("tobit_flex type 1 double censoring works", {
  set.seed(42)
  n <- 500; x <- rnorm(n)
  y_star <- 5 + 2*x + rnorm(n, sd = 2)
  y <- pmax(pmin(y_star, 10), 0)
  df <- data.frame(y = y, x = x)

  m <- tobit_flex(y ~ x, data = df, tobit_type = 1, left = 0, right = 10)
  expect_true(m$n_censored_left > 0 || m$n_censored_right > 0)
})

test_that("predict works for tobit_flex", {
  set.seed(42)
  n <- 200; x <- rnorm(n)
  y <- pmax(1 + x + rnorm(n), 0)
  df <- data.frame(y = y, x = x)

  m <- tobit_flex(y ~ x, data = df, tobit_type = 1, left = 0)
  p <- predict(m)
  expect_equal(length(p), n)
  expect_true(all(p >= 0))
})
