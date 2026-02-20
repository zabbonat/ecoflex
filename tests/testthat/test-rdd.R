test_that("rdd_flex sharp produces valid output", {
  set.seed(42)
  n <- 1000
  x <- runif(n, -1, 1)
  y <- 1 + 0.5 * x + 2 * (x >= 0) + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x)

  m <- rdd_flex(y ~ x, data = df, cutoff = 0, type = "sharp")
  expect_s3_class(m, "rdd_flex")
  expect_s3_class(m, "ecoflex")
  expect_true(!is.na(m$tau))
  expect_true(abs(m$tau - 2) < 1.5)
})

test_that("rdd_manipulation_test works", {
  set.seed(42)
  x <- rnorm(500)
  result <- rdd_manipulation_test(x, cutoff = 0)
  expect_true(!is.null(result$statistic))
  expect_true(!is.null(result$p.value))
})

test_that("rdd_placebo_cutoffs works", {
  set.seed(42)
  n <- 500
  x <- runif(n, -2, 2)
  y <- 1 + x + 1.5 * (x >= 0) + rnorm(n)
  df <- data.frame(y = y, x = x)
  result <- rdd_placebo_cutoffs(y ~ x, data = df, cutoff = 0,
                                 placebo_cutoffs = c(-1, 1))
  expect_true(nrow(result) == 2)
})
