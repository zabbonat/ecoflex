# Quick test — write output from within R
pkg_dir <- "c:/Users/Dilet/Desktop/Milano/Novelty/Causality/aRRpackage/R"
out_file <- "c:/Users/Dilet/Desktop/Milano/Novelty/Causality/aRRpackage/test_results.txt"
sink(out_file, split = TRUE)

suppressPackageStartupMessages({
  library(stats); library(Formula); library(sandwich)
  library(numDeriv); library(MASS); library(Matrix)
})

for (f in c("utils.R","generics.R","vcov_methods.R","latex_output.R",
            "hurdle_flex.R","heckman_flex.R","rdd_flex.R","ivgmm_flex.R",
            "did_flex.R","tobit_flex.R","panel_nl.R","switching_flex.R")) {
  source(file.path(pkg_dir, f))
}
cat("All sourced OK\n\n")

results <- list()

# ---- TEST 1: Tobit on mtcars ----
cat("=== TEST 1: Tobit Type I (mtcars, mpg censored at 15) ===\n")
tryCatch({
  m <- tobit_flex(mpg ~ hp + wt, data = mtcars, tobit_type = 1, left = 15)
  cat("Coefficients:\n"); print(coef(m))
  s <- summary(m)
  cat("\nSummary:\n"); print(s)
  results$tobit <- "PASS"
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); results$tobit <<- "FAIL" })

# ---- TEST 2: IV/2SLS ----
cat("\n=== TEST 2: IV/2SLS (simulated) ===\n")
tryCatch({
  set.seed(42); n <- 500
  z1 <- rnorm(n); z2 <- rnorm(n)
  x <- 0.5*z1 + 0.3*z2 + rnorm(n)
  y <- 1 + 2*x + rnorm(n)
  df <- data.frame(y=y, x=x, z1=z1, z2=z2)
  m <- ivgmm_flex(y ~ x | z1 + z2, data = df, method = "2sls")
  cat("Coefficients:\n"); print(coef(m))
  cat("Weak IV F:", m$diagnostics$weak_iv_test$first_stage_F, "\n")
  cat("Sargan p:", m$diagnostics$overid_test$p.value, "\n")
  s <- summary(m)
  cat("\nSummary:\n"); print(s)
  results$iv <- "PASS"
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); results$iv <<- "FAIL" })

# ---- TEST 3: RDD ----
cat("\n=== TEST 3: Sharp RDD (simulated) ===\n")
tryCatch({
  set.seed(123)
  x <- runif(1000, -1, 1)
  y <- 1 + 0.5*x + 3*(x >= 0) + rnorm(1000, sd = 0.5)
  df <- data.frame(y=y, x=x)
  m <- rdd_flex(y ~ x, data = df, cutoff = 0, type = "sharp")
  cat("tau =", m$tau, "\nse =", m$se_tau, "\n")
  s <- summary(m)
  cat("\nSummary:\n"); print(s)
  results$rdd <- "PASS"
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); results$rdd <<- "FAIL" })

# ---- TEST 4: Hurdle ----
cat("\n=== TEST 4: Hurdle (warpbreaks) ===\n")
tryCatch({
  data(warpbreaks)
  m <- hurdle_flex(breaks ~ wool + tension, data = warpbreaks,
                   threshold = 0, dist = "poisson")
  cat("Coefficients:\n"); print(coef(m))
  s <- summary(m)
  cat("\nSummary:\n"); print(s)
  results$hurdle <- "PASS"
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); results$hurdle <<- "FAIL" })

# ---- TEST 5: Heckman ----
cat("\n=== TEST 5: Heckman (twostep, simulated) ===\n")
tryCatch({
  set.seed(42); n <- 1000
  x1 <- rnorm(n); x2 <- rnorm(n); z <- rnorm(n)
  s <- as.integer(0.5 + 0.3*x1 + 0.5*z + rnorm(n) > 0)
  y <- ifelse(s==1, 1 + 2*x1 + 0.5*x2 + rnorm(n), 0)
  df <- data.frame(y=y, s=s, x1=x1, x2=x2, z=z)
  m <- heckman_flex(y | s ~ x1 + x2 | x1 + z, data = df, method = "twostep")
  cat("Coefficients:\n"); print(coef(m))
  s <- summary(m)
  cat("\nSummary:\n"); print(s)
  results$heckman <- "PASS"
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); results$heckman <<- "FAIL" })

# ---- TEST 6: Switching ----
cat("\n=== TEST 6: Switching (exogenous, mtcars) ===\n")
tryCatch({
  m <- switching_flex(mpg ~ hp + wt, data = mtcars,
                      type = "exogenous", regime_var = "am")
  cat("Coefficients:\n"); print(coef(m))
  s <- summary(m)
  cat("\nSummary:\n"); print(s)
  results$switching <- "PASS"
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); results$switching <<- "FAIL" })

# ---- TEST 7: summary(m, latex=TRUE) ----
cat("\n=== TEST 7: summary(m, latex=TRUE) ===\n")
tryCatch({
  summary(m, latex = TRUE)
  results$latex_summary <- "PASS"
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); results$latex_summary <<- "FAIL" })

# ---- RESULTS ----
cat("\n\n========== RESULTS ==========\n")
for (nm in names(results)) {
  cat(sprintf("  %-20s : %s\n", nm, results[[nm]]))
}
cat("=============================\n")

sink()
cat("Output written to", out_file, "\n")
