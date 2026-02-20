test_that("ivgmm_flex 2SLS returns correct class and structure", {
  set.seed(42)
  n  <- 300
  z  <- rnorm(n)
  x  <- 0.5 * z + rnorm(n)
  y  <- 1 + 2 * x + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  m <- ivgmm_flex(y ~ x | z, data = df, method = "2sls")

  expect_s3_class(m, "ivgmm_flex")
  expect_s3_class(m, "ecoflex")
  expect_true(is.numeric(coef(m)))
  expect_true(length(coef(m)) > 0)
  expect_true(is.matrix(vcov(m)))
  expect_equal(nobs(m), n)
})

test_that("ivgmm_flex summary works", {
  set.seed(1)
  n  <- 200; z <- rnorm(n); x <- 0.5 * z + rnorm(n); y <- 2 * x + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)
  m  <- ivgmm_flex(y ~ x | z, data = df, method = "2sls")
  s  <- summary(m)
  expect_true(!is.null(s))
})

test_that("ivgmm_flex coefficient close to true value", {
  set.seed(99)
  n     <- 1000
  z     <- rnorm(n)
  x     <- 0.8 * z + rnorm(n, sd = 0.5)
  y     <- 1 + 3 * x + rnorm(n)
  df    <- data.frame(y = y, x = x, z = z)
  m     <- ivgmm_flex(y ~ x | z, data = df, method = "2sls")
  est_x <- coef(m)["x"]
  expect_true(abs(est_x - 3) < 0.5)
})

test_that("ivgmm_flex errors on under-identified model", {
  set.seed(2)
  n  <- 100; y <- rnorm(n); x1 <- rnorm(n); x2 <- rnorm(n); z <- rnorm(n)
  df <- data.frame(y = y, x1 = x1, x2 = x2, z = z)
  expect_error(ivgmm_flex(y ~ x1 + x2 | z, data = df))
})

test_that("ivgmm_flex diagnostics are present", {
  set.seed(3)
  n  <- 200; z <- rnorm(n); x <- 0.6 * z + rnorm(n); y <- 2 * x + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)
  m  <- ivgmm_flex(y ~ x | z, data = df, method = "2sls")
  expect_true(!is.null(m$diagnostics))
  expect_true(!is.null(m$diagnostics$weak_iv_test))
})
