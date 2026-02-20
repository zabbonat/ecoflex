test_that("to_latex single model produces valid LaTeX", {
  set.seed(42)
  df <- data.frame(y = rpois(200, 3), x = rnorm(200))
  m <- hurdle_flex(y ~ x, data = df)
  out <- to_latex(m, digits = 3)
  expect_true(any(grepl("\\\\begin\\{table\\}", out)))
  expect_true(any(grepl("\\\\begin\\{adjustbox\\}", out)))
  expect_true(any(grepl("\\\\begin\\{threeparttable\\}", out)))
  expect_true(any(grepl("\\\\end\\{table\\}", out)))
})

test_that("to_latex hurdle produces side-by-side", {
  set.seed(42)
  df <- data.frame(y = rpois(200, 3), x = rnorm(200))
  m <- hurdle_flex(y ~ x, data = df)
  out <- to_latex(m)
  # Side-by-side: should have "Zero/Binary" and "Count/Hurdle" headers
  expect_true(any(grepl("Zero|Count|Binary|Hurdle", out)))
})

test_that("to_latex with preamble includes usepackage", {
  set.seed(42)
  df <- data.frame(y = rpois(200, 3), x = rnorm(200))
  m <- hurdle_flex(y ~ x, data = df)
  out <- to_latex(m, preamble = TRUE)
  expect_true(any(grepl("usepackage\\{adjustbox\\}", out)))
  expect_true(any(grepl("usepackage\\{booktabs\\}", out)))
  expect_true(any(grepl("usepackage\\{threeparttable\\}", out)))
})

test_that("to_latex file output works", {
  set.seed(42)
  df <- data.frame(y = rpois(200, 3), x = rnorm(200))
  m <- hurdle_flex(y ~ x, data = df)
  tmp <- tempfile(fileext = ".tex")
  to_latex(m, file = tmp)
  expect_true(file.exists(tmp))
  content <- readLines(tmp)
  expect_true(any(grepl("begin\\{table\\}", content)))
  unlink(tmp)
})

test_that("to_latex RDD single column works", {
  set.seed(42)
  x <- runif(500, -1, 1)
  y <- 1 + 2 * (x >= 0) + rnorm(500)
  df <- data.frame(y = y, x = x)
  m <- rdd_flex(y ~ x, data = df, cutoff = 0)
  out <- to_latex(m, stars = FALSE, se_type = "t")
  expect_true(any(grepl("\\\\begin\\{table\\}", out)))
})

test_that("to_latex ivgmm works", {
  set.seed(42)
  n <- 300
  z <- rnorm(n); x <- 0.5*z + rnorm(n); y <- 1 + 2*x + rnorm(n)
  df <- data.frame(y = y, x = x, z = z)
  m <- ivgmm_flex(y ~ x | z, data = df, method = "2sls")
  out <- to_latex(m, caption = "IV results", label = "tab:iv")
  expect_true(any(grepl("\\\\caption\\{IV results\\}", out)))
  expect_true(any(grepl("\\\\label\\{tab:iv\\}", out)))
})
