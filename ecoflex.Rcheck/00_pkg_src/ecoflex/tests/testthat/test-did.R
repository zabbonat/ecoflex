test_that("did_flex TWFE works", {
  set.seed(42)
  n_units <- 50; n_periods <- 10
  panel <- expand.grid(id = 1:n_units, time = 1:n_periods)
  panel$treated <- as.integer(panel$id > 25 & panel$time > 5)
  panel$y <- 1 + 0.5 * panel$time + 2 * panel$treated + rnorm(nrow(panel))

  m <- did_flex(y ~ 1, data = panel, id_var = "id", time_var = "time",
                treat_var = "treated", estimator = "twfe")
  expect_s3_class(m, "did_flex")
  expect_true(!is.na(m$att))
})

test_that("did_flex Callaway-Sant'Anna works", {
  set.seed(42)
  n_units <- 60; n_periods <- 8
  panel <- expand.grid(id = 1:n_units, time = 1:n_periods)
  panel$cohort <- ifelse(panel$id <= 20, 5,
                          ifelse(panel$id <= 40, 7, Inf))
  panel$treated <- as.integer(panel$time >= panel$cohort)
  panel$y <- 1 + 0.3 * panel$time + 1.5 * panel$treated + rnorm(nrow(panel))

  m <- did_flex(y ~ 1, data = panel, id_var = "id", time_var = "time",
                treat_var = "treated", cohort_var = "cohort",
                estimator = "callaway_santanna")
  expect_s3_class(m, "did_flex")
})

test_that("did_compare works", {
  set.seed(42)
  n_units <- 40; n_periods <- 6
  panel <- expand.grid(id = 1:n_units, time = 1:n_periods)
  panel$treated <- as.integer(panel$id > 20 & panel$time > 3)
  panel$y <- 1 + panel$time + 2 * panel$treated + rnorm(nrow(panel))

  m1 <- did_flex(y ~ 1, data = panel, id_var = "id", time_var = "time",
                 treat_var = "treated", estimator = "twfe")
  m2 <- did_flex(y ~ 1, data = panel, id_var = "id", time_var = "time",
                 treat_var = "treated", estimator = "stacked")

  comp <- did_compare(m1, m2, labels = c("TWFE", "Stacked"))
  expect_s3_class(comp, "did_comparison")
  expect_equal(nrow(comp$comparison), 2)
})
