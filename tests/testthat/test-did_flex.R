test_that("did_flex twfe returns correct class and structure", {
  set.seed(42)
  panel        <- expand.grid(id = 1:40, time = 1:8)
  panel$treat  <- as.integer(panel$id > 20 & panel$time >= 5)
  panel$y      <- 1 + 0.5 * panel$treat + rnorm(nrow(panel))

  m <- did_flex(y ~ 1, data = panel,
                id_var = "id", time_var = "time", treat_var = "treat",
                estimator = "twfe")

  expect_s3_class(m, "did_flex")
  expect_s3_class(m, "ecoflex")
  expect_true(is.numeric(coef(m)))
})

test_that("did_flex summary works", {
  set.seed(1)
  panel       <- expand.grid(id = 1:30, time = 1:6)
  panel$treat <- as.integer(panel$id > 15 & panel$time >= 4)
  panel$y     <- 0.8 * panel$treat + rnorm(nrow(panel))
  m  <- did_flex(y ~ 1, data = panel,
                 id_var = "id", time_var = "time", treat_var = "treat",
                 estimator = "twfe")
  s  <- summary(m)
  expect_true(!is.null(s))
})

test_that("did_flex nobs equals panel rows", {
  set.seed(2)
  panel       <- expand.grid(id = 1:20, time = 1:5)
  panel$treat <- as.integer(panel$id > 10 & panel$time >= 3)
  panel$y     <- panel$treat + rnorm(nrow(panel))
  m <- did_flex(y ~ 1, data = panel,
                id_var = "id", time_var = "time", treat_var = "treat",
                estimator = "twfe")
  expect_equal(nobs(m), nrow(panel))
})
