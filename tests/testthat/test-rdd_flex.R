test_that("rdd_flex sharp returns correct class and structure", {
  set.seed(123)
  n  <- 400
  x  <- runif(n, -1, 1)
  y  <- 1 + 0.5 * x + 2 * (x >= 0) + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x)

  m <- rdd_flex(y ~ x, data = df, cutoff = 0, type = "sharp")

  expect_s3_class(m, "rdd_flex")
  expect_s3_class(m, "ecoflex")
  expect_true(is.numeric(coef(m)))
})

test_that("rdd_flex summary works", {
  set.seed(1)
  n  <- 300
  x  <- runif(n, -1, 1)
  y  <- 2 * (x >= 0) + rnorm(n, sd = 0.3)
  df <- data.frame(y = y, x = x)
  m  <- rdd_flex(y ~ x, data = df, cutoff = 0)
  s  <- summary(m)
  expect_true(!is.null(s))
})

test_that("rdd_flex treatment effect is in plausible range", {
  set.seed(42)
  n    <- 600
  x    <- runif(n, -1, 1)
  true_effect <- 3
  y    <- true_effect * (x >= 0) + rnorm(n, sd = 0.4)
  df   <- data.frame(y = y, x = x)
  m    <- rdd_flex(y ~ x, data = df, cutoff = 0)
  est  <- coef(m)["RDD_effect"]
  expect_true(abs(est - true_effect) < 2.0)
})

test_that("rdd_manipulation_test returns a list", {
  set.seed(5)
  x    <- runif(300, -1, 1)
  res  <- rdd_manipulation_test(x, cutoff = 0)
  expect_true(is.list(res))
})
