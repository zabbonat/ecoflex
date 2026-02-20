## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)
library(ecoflex)

## ----hurdle-example, eval = FALSE---------------------------------------------
# set.seed(42)
# n <- 1000
# x <- rnorm(n)
# y <- rpois(n, exp(0.5 + 0.3 * x))
# 
# df <- data.frame(y = y, x = x)
# 
# # Joint ML hurdle with threshold = 0
# m <- hurdle_flex(y ~ x, data = df, threshold = 0)
# summary(m)

## ----rdd-example, eval = FALSE------------------------------------------------
# x <- runif(1000, -1, 1)
# y <- 1 + 0.5 * x + 2 * (x >= 0) + rnorm(1000, sd = 0.5)
# df <- data.frame(y = y, x = x)
# 
# m <- rdd_flex(y ~ x, data = df, cutoff = 0, type = "sharp")
# summary(m)
# 
# # Diagnostics
# rdd_manipulation_test(df$x, cutoff = 0)

## ----iv-example, eval = FALSE-------------------------------------------------
# n <- 500
# z <- rnorm(n); x <- 0.5 * z + rnorm(n)
# y <- 1 + 2 * x + rnorm(n)
# df <- data.frame(y = y, x = x, z = z)
# 
# m <- ivgmm_flex(y ~ x | z, data = df, method = "2sls")
# summary(m)
# m$diagnostics$weak_iv_test

## ----vcov-example, eval = FALSE-----------------------------------------------
# # Cluster-robust SE
# V_cl <- flex_vcov(m, type = "cluster", cluster = df$cluster_id)
# 
# # Bootstrap SE
# V_boot <- flex_vcov(m, type = "bootstrap", R = 1000)

