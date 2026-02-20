test_that("switching_flex exogenous works", {
  set.seed(42)
  n <- 300
  regime <- sample(0:1, n, replace = TRUE)
  x <- rnorm(n)
  y <- ifelse(regime == 1, 2 + 1.5*x + rnorm(n), -1 + 0.5*x + rnorm(n))
  df <- data.frame(y = y, x = x, regime = regime)

  m <- switching_flex(y ~ x, data = df, type = "exogenous", regime_var = "regime")
  expect_s3_class(m, "switching_flex")
  expect_s3_class(m, "ecoflex")
  expect_true(is.finite(m$logLik))
})

test_that("switching_flex endogenous works", {
  set.seed(42)
  n <- 400
  z <- rnorm(n); x <- rnorm(n)
  sel <- as.integer(0.5*z + rnorm(n) > 0)
  y <- ifelse(sel == 1, 3 + x + rnorm(n), 1 + 0.5*x + rnorm(n))
  df <- data.frame(y = y, x = x, sel = sel, z = z)

  m <- switching_flex(y ~ x | sel ~ z, data = df, type = "endogenous")
  expect_s3_class(m, "switching_flex")
  expect_true(length(coef(m)) > 0)
})

test_that("switching_flex Markov works", {
  set.seed(42)
  n <- 200
  x <- matrix(1, n, 1)  # intercept only
  # Generate from 2-state Markov switching
  state <- numeric(n); state[1] <- 1
  for (t in 2:n) state[t] <- sample(1:2, 1, prob = if (state[t-1]==1) c(0.9,0.1) else c(0.2,0.8))
  y <- ifelse(state == 1, 2 + rnorm(n, sd = 0.5), -1 + rnorm(n, sd = 1))
  df <- data.frame(y = y, const = 1)

  m <- switching_flex(y ~ const, data = df, type = "markov", n_regimes = 2)
  expect_s3_class(m, "switching_flex")
  expect_true(!is.null(m$filtered_probs))
  expect_equal(ncol(m$filtered_probs), 2)
})
