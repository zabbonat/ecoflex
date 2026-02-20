test_that("heckman_flex twostep returns correct class and structure", {
  set.seed(42)
  n  <- 500
  x1 <- rnorm(n); z <- rnorm(n)
  s  <- as.integer(0.3 * x1 + 0.5 * z + rnorm(n) > 0)
  y  <- ifelse(s == 1, 1 + 2 * x1 + rnorm(n), NA)
  df <- data.frame(y = y, s = s, x1 = x1, z = z)

  m <- heckman_flex(y | s ~ x1 | x1 + z, data = df, method = "twostep")

  expect_s3_class(m, "heckman_flex")
  expect_s3_class(m, "ecoflex")
  expect_true(is.numeric(coef(m)))
  expect_true(length(coef(m)) > 0)
})

test_that("heckman_flex summary works", {
  set.seed(1)
  n  <- 400
  x1 <- rnorm(n); z <- rnorm(n)
  s  <- as.integer(0.4 * x1 + 0.4 * z + rnorm(n) > 0)
  y  <- ifelse(s == 1, 0.5 + 1.5 * x1 + rnorm(n), NA)
  df <- data.frame(y = y, s = s, x1 = x1, z = z)
  m <- heckman_flex(y | s ~ x1 | x1 + z, data = df, method = "twostep")
  s2 <- summary(m)
  expect_true(!is.null(s2))
})

test_that("heckman_flex predict returns numeric", {
  set.seed(2)
  n  <- 300
  x1 <- rnorm(n); z <- rnorm(n)
  s  <- as.integer(0.5 * x1 + 0.3 * z + rnorm(n) > 0)
  y  <- ifelse(s == 1, 1 + x1 + rnorm(n), NA)
  df <- data.frame(y = y, s = s, x1 = x1, z = z)
  m  <- heckman_flex(y | s ~ x1 | x1 + z, data = df, method = "twostep")
  p  <- predict(m, type = "response")
  expect_true(is.numeric(p))
})
