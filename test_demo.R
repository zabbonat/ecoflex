# =============================================================================
# ecoflex — Demo tests on real R datasets
# =============================================================================

cat("=== Loading dependencies ===\n")
suppressPackageStartupMessages({
  library(stats)
  library(Formula)
  library(sandwich)
  library(numDeriv)
  library(MASS)
  library(Matrix)
})

cat("=== Sourcing ecoflex modules ===\n")
pkg_dir <- "c:/Users/Dilet/Desktop/Milano/Novelty/Causality/aRRpackage/R"
source(file.path(pkg_dir, "utils.R"))
source(file.path(pkg_dir, "generics.R"))
source(file.path(pkg_dir, "vcov_methods.R"))
source(file.path(pkg_dir, "latex_output.R"))
source(file.path(pkg_dir, "hurdle_flex.R"))
source(file.path(pkg_dir, "heckman_flex.R"))
source(file.path(pkg_dir, "rdd_flex.R"))
source(file.path(pkg_dir, "ivgmm_flex.R"))
source(file.path(pkg_dir, "did_flex.R"))
source(file.path(pkg_dir, "tobit_flex.R"))
source(file.path(pkg_dir, "panel_nl.R"))
source(file.path(pkg_dir, "switching_flex.R"))

# =============================================================================
# TEST 1: Hurdle Model — warpbreaks dataset (count data)
# =============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("TEST 1: Hurdle Model on warpbreaks\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

data(warpbreaks)
m_hurdle <- hurdle_flex(breaks ~ wool + tension, data = warpbreaks,
                        threshold = 0, dist = "poisson")
cat("Summary:\n")
print(summary(m_hurdle))
cat("\nPredictions (first 6):\n")
print(head(predict(m_hurdle)))
cat("\nLaTeX output:\n")
summary(m_hurdle, latex = TRUE)

# =============================================================================
# TEST 2: Hurdle Model — mtcars (treat carb as count, threshold = 2)
# =============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("TEST 2: Hurdle on mtcars (carb, threshold=2)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

m_hurdle2 <- hurdle_flex(carb ~ hp + wt + mpg, data = mtcars,
                         threshold = 2, dist = "poisson")
print(summary(m_hurdle2))

# =============================================================================
# TEST 3: Tobit Model — mtcars (mpg censored at 15)
# =============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("TEST 3: Tobit Type I on mtcars (mpg censored at 15)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

m_tobit <- tobit_flex(mpg ~ hp + wt + disp, data = mtcars,
                      tobit_type = 1, left = 15)
print(summary(m_tobit))
cat("\nLaTeX output:\n")
summary(m_tobit, latex = TRUE, caption = "Tobit Model -- mtcars", label = "tab:tobit")

# =============================================================================
# TEST 4: Tobit with t-distribution
# =============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("TEST 4: Tobit Type I with t-errors on mtcars\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

m_tobit_t <- tobit_flex(mpg ~ hp + wt, data = mtcars,
                        tobit_type = 1, left = 15, error_dist = "t")
print(summary(m_tobit_t))

# =============================================================================
# TEST 5: IV/GMM — simulated data (endogeneity)
# =============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("TEST 5: IV/2SLS on simulated data\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

set.seed(42)
n <- 500
z1 <- rnorm(n); z2 <- rnorm(n)
x_endog <- 0.5 * z1 + 0.3 * z2 + rnorm(n)
y_iv <- 1 + 2 * x_endog + rnorm(n)

df_iv <- data.frame(y = y_iv, x = x_endog, z1 = z1, z2 = z2)
m_2sls <- ivgmm_flex(y ~ x | z1 + z2, data = df_iv, method = "2sls")
print(summary(m_2sls))
cat("\nDiagnostics:\n")
cat("  Weak IV F-stat:", m_2sls$diagnostics$weak_iv_test$first_stage_F, "\n")
cat("  Sargan-Hansen J:", m_2sls$diagnostics$overid_test$statistic,
    " p =", m_2sls$diagnostics$overid_test$p.value, "\n")
cat("\nLaTeX output:\n")
summary(m_2sls, latex = TRUE)

# =============================================================================
# TEST 6: RDD — simulated sharp design
# =============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("TEST 6: Sharp RDD on simulated data\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

set.seed(123)
x_rdd <- runif(1000, -1, 1)
y_rdd <- 1 + 0.5 * x_rdd + 3 * (x_rdd >= 0) + rnorm(1000, sd = 0.5)
df_rdd <- data.frame(y = y_rdd, x = x_rdd)

m_rdd <- rdd_flex(y ~ x, data = df_rdd, cutoff = 0, type = "sharp")
print(summary(m_rdd))
cat("\nManipulation test:\n")
print(rdd_manipulation_test(df_rdd$x, cutoff = 0))
cat("\nLaTeX output:\n")
summary(m_rdd, latex = TRUE)

# =============================================================================
# TEST 7: Heckman Selection — simulated
# =============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("TEST 7: Heckman Selection (two-step)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

set.seed(42)
n <- 1000
x1 <- rnorm(n); x2 <- rnorm(n); z <- rnorm(n)
s_star <- 0.5 + 0.3 * x1 + 0.5 * z + rnorm(n)
s <- as.integer(s_star > 0)
y_star <- 1 + 2 * x1 + 0.5 * x2 + rnorm(n)
y <- ifelse(s == 1, y_star, 0)

df_heck <- data.frame(y = y, s = s, x1 = x1, x2 = x2, z = z)
# Formula: y | s ~ outcome_vars | selection_vars
m_heck <- heckman_flex(y | s ~ x1 + x2 | x1 + z, data = df_heck, method = "twostep")
print(summary(m_heck))
cat("\nLaTeX (side-by-side Selection/Outcome):\n")
summary(m_heck, latex = TRUE)

# =============================================================================
# TEST 8: Panel Poisson FE — simulated panel
# =============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("TEST 8: Panel Poisson FE\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

set.seed(42)
N <- 50; T_periods <- 10
panel_df <- expand.grid(id = 1:N, time = 1:T_periods)
panel_df$x1 <- rnorm(nrow(panel_df))
panel_df$x2 <- rnorm(nrow(panel_df))
alpha_i <- rep(rnorm(N, sd = 0.3), each = T_periods)
panel_df$y <- rpois(nrow(panel_df), exp(0.5 + 0.3 * panel_df$x1 - 0.2 * panel_df$x2 + alpha_i))

m_panel <- panel_nl(y ~ x1 + x2, data = panel_df, id = "id", time = "time",
                    family = "poisson", model = "fe", effect = "individual")
print(summary(m_panel))

# =============================================================================
# TEST 9: Switching Regression (exogenous) — mtcars by transmission
# =============================================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("TEST 9: Exogenous Switching on mtcars (by am)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

m_switch <- switching_flex(mpg ~ hp + wt, data = mtcars,
                           type = "exogenous", regime_var = "am")
print(summary(m_switch))
cat("\nLaTeX (side-by-side Regime 0 / Regime 1):\n")
summary(m_switch, latex = TRUE)

# =============================================================================
cat("\n\n", paste(rep("*", 60), collapse = ""), "\n")
cat("ALL TESTS COMPLETED SUCCESSFULLY!\n")
cat(paste(rep("*", 60), collapse = ""), "\n")
