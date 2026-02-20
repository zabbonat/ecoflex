test_that("tobit_flex type I returns correct class and structure", {
  m <- tobit_flex(mpg ~ hp + wt, data = mtcars, tobit_type = 1, left = 15)

  expect_s3_class(m, "tobit_flex")
  expect_s3_class(m, "ecoflex")
  expect_true(is.numeric(coef(m)))
  expect_true(length(coef(m)) > 0)
  expect_true(is.matrix(vcov(m)))
  expect_equal(nobs(m), nrow(mtcars))
})

test_that("tobit_flex summary works", {
  m <- tobit_flex(mpg ~ hp + wt, data = mtcars, tobit_type = 1, left = 15)
  s <- summary(m)
  expect_true(!is.null(s))
})

test_that("tobit_flex predict returns numeric vector", {
  m <- tobit_flex(mpg ~ hp + wt, data = mtcars, tobit_type = 1, left = 15)
  p <- predict(m, type = "response")
  expect_true(is.numeric(p))
  expect_equal(length(p), nrow(mtcars))
})

test_that("tobit_flex validates tobit_type argument", {
  expect_error(tobit_flex(mpg ~ hp, data = mtcars, tobit_type = 0))
  expect_error(tobit_flex(mpg ~ hp, data = mtcars, tobit_type = 6))
})

test_that("tobit_flex logLik is finite", {
  m  <- tobit_flex(mpg ~ hp + wt, data = mtcars, tobit_type = 1, left = 15)
  ll <- logLik(m)
  expect_true(is.finite(as.numeric(ll)))
})
