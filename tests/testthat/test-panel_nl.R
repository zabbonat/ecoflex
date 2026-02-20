test_that("panel_nl Poisson FE returns correct class and structure", {
  set.seed(42)
  panel     <- expand.grid(id = 1:20, time = 1:5)
  panel$x   <- rnorm(nrow(panel))
  panel$y   <- rpois(nrow(panel), exp(0.5 + 0.3 * panel$x))

  m <- panel_nl(y ~ x, data = panel, id = "id", time = "time",
                family = "poisson", model = "fe")

  expect_s3_class(m, "panel_nl")
  expect_s3_class(m, "ecoflex")
  expect_true(is.numeric(coef(m)))
  expect_true(length(coef(m)) > 0)
})

test_that("panel_nl summary works", {
  set.seed(1)
  panel   <- expand.grid(id = 1:15, time = 1:4)
  panel$x <- rnorm(nrow(panel))
  panel$y <- rpois(nrow(panel), exp(0.3 * panel$x))
  m  <- panel_nl(y ~ x, data = panel, id = "id", time = "time",
                 family = "poisson", model = "fe")
  s  <- summary(m)
  expect_true(!is.null(s))
})

test_that("panel_nl predict returns numeric", {
  set.seed(2)
  panel   <- expand.grid(id = 1:20, time = 1:5)
  panel$x <- rnorm(nrow(panel))
  panel$y <- rpois(nrow(panel), exp(0.4 * panel$x))
  m  <- panel_nl(y ~ x, data = panel, id = "id", time = "time",
                 family = "poisson", model = "fe")
  p  <- predict(m, type = "response")
  expect_true(is.numeric(p))
  expect_equal(length(p), nrow(panel))
})

test_that("panel_nl logit FE runs without error", {
  set.seed(3)
  panel     <- expand.grid(id = 1:20, time = 1:5)
  panel$x   <- rnorm(nrow(panel))
  panel$y_b <- as.integer(panel$x + rnorm(nrow(panel)) > 0)
  expect_no_error(
    panel_nl(y_b ~ x, data = panel, id = "id", time = "time",
             family = "logit", model = "fe")
  )
})
