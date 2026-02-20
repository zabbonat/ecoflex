test_that("ivgmm_flex 2SLS produces valid output", {
  set.seed(42)
  n <- 500
  z1 <- rnorm(n); z2 <- rnorm(n)
  x <- 0.5 * z1 + 0.3 * z2 + rnorm(n)
  y <- 1 + 2 * x + rnorm(n)
  df <- data.frame(y = y, x = x, z1 = z1, z2 = z2, const = 1)

  m <- ivgmm_flex(y ~ x + const | z1 + z2 + const, data = df, method = "2sls")
  expect_s3_class(m, "ivgmm_flex")
  expect_s3_class(m, "ecoflex")
  expect_true(length(coef(m)) > 0)
})

test_that("2SLS estimate is close to true value", {
  set.seed(42)
  n <- 1000
  z <- rnorm(n); x <- 0.8 * z + rnorm(n)
  y <- 1 + 2 * x + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  m <- ivgmm_flex(y ~ x | z, data = df, method = "2sls")
  expect_true(abs(coef(m)[2] - 2) < 1)
})

test_that("GMM iterative converges", {
  set.seed(42)
  n <- 1000
  z1 <- rnorm(n); z2 <- rnorm(n)
  x <- 0.5 * z1 + 0.3 * z2 + rnorm(n)
  y <- 1 + 2 * x + rnorm(n)
  df <- data.frame(y = y, x = x, z1 = z1, z2 = z2)

  m <- ivgmm_flex(y ~ x | z1 + z2, data = df, method = "gmm_iterative")
  expect_true(abs(coef(m)[2] - 2) < 0.5)
})

test_that("LIML produces valid output", {
  set.seed(42)
  n <- 500
  z <- rnorm(n); x <- 0.5 * z + rnorm(n)
  y <- 1 + 2 * x + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)

  m <- ivgmm_flex(y ~ x | z, data = df, method = "liml")
  expect_s3_class(m, "ivgmm_flex")
})

test_that("diagnostics are computed", {
  set.seed(42)
  n <- 500
  z1 <- rnorm(n); z2 <- rnorm(n)
  x <- 0.5 * z1 + 0.3 * z2 + rnorm(n)
  y <- 1 + x + rnorm(n)
  df <- data.frame(y = y, x = x, z1 = z1, z2 = z2)

  m <- ivgmm_flex(y ~ x | z1 + z2, data = df, method = "2sls")
  expect_true(!is.null(m$diagnostics$weak_iv_test))
  expect_true(!is.null(m$diagnostics$endogeneity_test))
  expect_true(!is.null(m$diagnostics$first_stage_F))
  expect_true(!is.null(m$diagnostics$overid_test))
})

test_that("predict works for ivgmm_flex", {
  set.seed(42)
  n <- 200
  z <- rnorm(n); x <- 0.5 * z + rnorm(n)
  y <- 1 + x + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)
  m <- ivgmm_flex(y ~ x | z, data = df, method = "2sls")
  p <- predict(m)
  expect_equal(length(p), n)
})
