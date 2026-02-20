test_that("switching_flex exogenous returns correct class and structure", {
  m <- switching_flex(mpg ~ hp + wt, data = mtcars,
                      type = "exogenous", regime_var = "am")

  expect_s3_class(m, "switching_flex")
  expect_s3_class(m, "ecoflex")
  expect_true(is.numeric(coef(m)))
  expect_true(length(coef(m)) > 0)
})

test_that("switching_flex summary works", {
  m <- switching_flex(mpg ~ hp + wt, data = mtcars,
                      type = "exogenous", regime_var = "am")
  s <- summary(m)
  expect_true(!is.null(s))
})

test_that("switching_flex predict returns numeric", {
  m <- switching_flex(mpg ~ hp + wt, data = mtcars,
                      type = "exogenous", regime_var = "am")
  p <- predict(m, type = "response")
  expect_true(is.numeric(p))
  expect_equal(length(p), nrow(mtcars))
})

test_that("switching_flex requires regime_var for exogenous type", {
  expect_error(
    switching_flex(mpg ~ hp, data = mtcars, type = "exogenous")
  )
})
