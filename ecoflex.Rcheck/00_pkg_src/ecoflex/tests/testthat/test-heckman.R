test_that("heckman_flex twostep produces valid output", {
  set.seed(42)
  n <- 500
  z <- rnorm(n); x <- rnorm(n); kids <- rbinom(n, 2, 0.3)
  sel <- as.integer(0.5 + 0.3*x - 0.5*kids + rnorm(n) > 0)
  wage <- ifelse(sel == 1, 10 + 2*x + rnorm(n), NA)
  df <- data.frame(wage = ifelse(is.na(wage), 0, wage), sel = sel, x = x, kids = kids)

  m <- heckman_flex(wage ~ x | sel ~ x + kids, data = df, method = "twostep")
  expect_s3_class(m, "heckman_flex")
  expect_s3_class(m, "ecoflex")
  expect_true(length(coef(m)) > 0)
  expect_true(is.finite(m$logLik))
})

test_that("heckman_flex ML normal works", {
  set.seed(42)
  n <- 400
  x <- rnorm(n); kids <- rbinom(n, 2, 0.3)
  sel <- as.integer(0.5 + 0.3*x - 0.5*kids + rnorm(n) > 0)
  wage <- ifelse(sel == 1, 10 + 2*x + rnorm(n), 0)
  df <- data.frame(wage = wage, sel = sel, x = x, kids = kids)

  m <- heckman_flex(wage ~ x | sel ~ x + kids, data = df, method = "ml")
  expect_s3_class(m, "heckman_flex")
  expect_true(!is.null(m$sigma))
  expect_true(!is.null(m$rho))
})

test_that("heckman_flex semiparametric works", {
  set.seed(42)
  n <- 400
  x <- rnorm(n); z <- rnorm(n)
  sel <- as.integer(0.5*x + 0.3*z + rnorm(n) > 0)
  wage <- ifelse(sel == 1, 5 + x + rnorm(n), 0)
  df <- data.frame(wage = wage, sel = sel, x = x, z = z)

  m <- heckman_flex(wage ~ x | sel ~ x + z, data = df, method = "semipar")
  expect_s3_class(m, "heckman_flex")
  expect_true(!is.null(m$poly_order))
})

test_that("predict works for heckman_flex", {
  set.seed(42)
  n <- 300
  x <- rnorm(n); z <- rnorm(n)
  sel <- as.integer(0.5*x + 0.3*z + rnorm(n) > 0)
  wage <- ifelse(sel == 1, 5 + x + rnorm(n), 0)
  df <- data.frame(wage = wage, sel = sel, x = x, z = z)

  m <- heckman_flex(wage ~ x | sel ~ x + z, data = df, method = "twostep")
  p <- predict(m, type = "response")
  expect_true(length(p) == n)
})
