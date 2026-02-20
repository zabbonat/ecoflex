test_that("panel_nl Poisson FE works", {
  set.seed(42)
  n_units <- 30; n_periods <- 5
  panel <- expand.grid(id = 1:n_units, time = 1:n_periods)
  panel$x <- rnorm(nrow(panel))
  panel$y <- rpois(nrow(panel), exp(0.5 + 0.3 * panel$x))

  m <- panel_nl(y ~ x, data = panel, id = "id", time = "time",
                family = "poisson", model = "fe")
  expect_s3_class(m, "panel_nl")
  expect_s3_class(m, "ecoflex")
  expect_true(length(coef(m)) > 0)
})

test_that("panel_nl correlated RE works", {
  set.seed(42)
  n_units <- 30; n_periods <- 5
  panel <- expand.grid(id = 1:n_units, time = 1:n_periods)
  panel$x <- rnorm(nrow(panel))
  panel$y <- rpois(nrow(panel), exp(0.5 + 0.3 * panel$x))

  m <- panel_nl(y ~ x, data = panel, id = "id", time = "time",
                family = "poisson", model = "correlated_re")
  expect_s3_class(m, "panel_nl")
})

test_that("panel_nl logit FE works", {
  set.seed(42)
  n_units <- 30; n_periods <- 5
  panel <- expand.grid(id = 1:n_units, time = 1:n_periods)
  panel$x <- rnorm(nrow(panel))
  panel$y <- rbinom(nrow(panel), 1, plogis(0.5 + 0.3 * panel$x))

  m <- panel_nl(y ~ x, data = panel, id = "id", time = "time",
                family = "logit", model = "fe")
  expect_s3_class(m, "panel_nl")
})
