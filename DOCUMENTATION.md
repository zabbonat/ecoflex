# ecoflex — Full Documentation & Examples

This guide provides **complete, runnable examples** for every function in
`ecoflex`. Each example generates its own data, so you can copy-paste any
block into R and run it immediately.

```r
library(ecoflex)
```

---

## Table of Contents

1. [Hurdle Models — `hurdle_flex()`](#1-hurdle-models)
2. [Heckman Selection — `heckman_flex()`](#2-heckman-selection-models)
3. [Regression Discontinuity — `rdd_flex()`](#3-regression-discontinuity)
4. [Instrumental Variables / GMM — `ivgmm_flex()`](#4-instrumental-variables--gmm)
5. [Difference-in-Differences — `did_flex()`](#5-difference-in-differences)
6. [Tobit Models — `tobit_flex()`](#6-tobit-models)
7. [Nonlinear Panel Data — `panel_nl()`](#7-nonlinear-panel-data)
8. [Switching Regression — `switching_flex()`](#8-switching-regression)
9. [Flexible Standard Errors — `flex_vcov()`](#9-flexible-standard-errors)
10. [LaTeX Output — `to_latex()`](#10-latex-output)

---

## 1. Hurdle Models

`hurdle_flex()` estimates count data models where the generating process
differs *above* vs *at-or-below* the hurdle. Unlike standard zero-inflated
models, the threshold can be any non-negative integer.

### Full Signature

```r
hurdle_flex(formula, data,
  threshold = 0,                       # hurdle threshold (any integer >= 0)
  dist      = "poisson",               # "poisson" | "negbin" | "geometric"
  zero_dist = "binomial",              # "binomial" (logit) | "probit"
  method    = "ml",                    # "ml" (joint MLE) | "two-step"
  start     = NULL,                    # optional numeric vector of starting values
  latex     = FALSE                    # if TRUE, prints LaTeX table
)
```

### Example: Basic Hurdle (threshold = 0)

```r
set.seed(42)
n <- 1000
x <- rnorm(n)
y <- rpois(n, exp(0.5 + 0.3 * x))

df <- data.frame(y = y, x = x)

m <- hurdle_flex(y ~ x, data = df, threshold = 0, dist = "poisson")
summary(m)
coef(m)         # named coefficient vector
logLik(m)       # log-likelihood
m$AIC           # AIC
```

### Example: Higher Threshold + Negative Binomial

```r
m2 <- hurdle_flex(y ~ x, data = df, threshold = 2, dist = "negbin")
summary(m2)
```

### Example: Two-Step Estimation with Probit

```r
m3 <- hurdle_flex(y ~ x, data = df, method = "two-step", zero_dist = "probit")
summary(m3)
```

---

## 2. Heckman Selection Models

`heckman_flex()` corrects for sample selection bias using the classic
Heckman two-equation framework. Supports ML, Heckman two-step, Lee's
robust two-step, and semiparametric correction.

### Formula Syntax

```
outcome | selection ~ outcome_covariates | selection_covariates
```

The `|` separates the two equations on both the LHS (dependent variables)
and RHS (covariates).

### Full Signature

```r
heckman_flex(formula, data,
  method          = "ml",              # "ml" | "twostep" | "twostep_robust" | "semipar"
  selection_type  = "binary",          # "binary" | "ordered" | "multinomial"
  error_dist      = "normal",          # "normal" | "t" | "copula"
  copula_family   = "gaussian",        # "gaussian" | "frank" | "clayton" | "gumbel"
  heteroscedastic = NULL,              # formula for heteroscedastic errors
  vcov            = "standard",        # "standard" | "robust" | "bootstrap"
  latex           = FALSE
)
```

### Example: Two-Step Heckman

```r
set.seed(42)
n <- 500
x <- rnorm(n)
kids <- rbinom(n, 2, 0.3)
sel <- as.integer(0.5 + 0.3 * x - 0.5 * kids + rnorm(n) > 0)
wage <- ifelse(sel == 1, 10 + 2 * x + rnorm(n), 0)

df <- data.frame(wage = wage, sel = sel, x = x, kids = kids)

m <- heckman_flex(wage | sel ~ x | x + kids, data = df, method = "twostep")
summary(m)

# Key outputs:
m$sigma    # residual SD
m$rho      # correlation between errors (selection bias indicator)
```

### Example: ML with t-distributed Errors

```r
m2 <- heckman_flex(wage | sel ~ x | x + kids, data = df,
                   method = "ml", error_dist = "t")
summary(m2)
```

### Prediction

```r
# Predicted values (E[y | observed])
p <- predict(m, type = "response")
head(p)
```

---

## 3. Regression Discontinuity

`rdd_flex()` implements sharp and fuzzy RDD with local polynomial
regression, MSE-optimal bandwidth selection, and various kernel functions.

### Full Signature

```r
rdd_flex(formula, data,
  cutoff           = 0,               # cutoff value for the running variable
  type             = "sharp",         # "sharp" | "fuzzy"
  fuzzy_treatment  = NULL,            # column name of actual treatment (fuzzy only)
  bandwidth        = "mserd",         # "mserd" (auto) or numeric value
  kernel           = "triangular",    # "triangular" | "epanechnikov" | "uniform" | "gaussian"
  polynomial       = 1,              # polynomial order (1 = local linear)
  polynomial_right = NULL,           # separate order for right side of cutoff
  discrete_running = FALSE,          # TRUE if running var is discrete
  cluster          = NULL,           # cluster variable name
  covs_adjust      = "none",         # "none" | "linear" | "interacted"
  masspoints       = "adjust",       # "adjust" | "check" | "off"
  all_bandwidths   = FALSE,          # return results for multiple bandwidths
  latex            = FALSE
)
```

### Example: Sharp RDD

```r
set.seed(42)
n <- 1000
x <- runif(n, -1, 1)
y <- 1 + 0.5 * x + 2 * (x >= 0) + rnorm(n, sd = 0.5)
df <- data.frame(y = y, x = x)

m <- rdd_flex(y ~ x, data = df, cutoff = 0, type = "sharp")
summary(m)

# Treatment effect at the cutoff:
m$tau       # point estimate
m$se_tau    # standard error
```

### Example: Different Kernels

```r
m_tri <- rdd_flex(y ~ x, data = df, cutoff = 0, kernel = "triangular")
m_epa <- rdd_flex(y ~ x, data = df, cutoff = 0, kernel = "epanechnikov")
m_uni <- rdd_flex(y ~ x, data = df, cutoff = 0, kernel = "uniform")
```

### RDD Diagnostics

```r
# McCrary/Cattaneo manipulation test (is the running variable manipulated?)
rdd_manipulation_test(df$x, cutoff = 0, method = "cattaneo")

# Placebo cutoffs (should NOT find effects at non-true cutoffs)
rdd_placebo_cutoffs(y ~ x, data = df, cutoff = 0)

# Balance test (covariates should be smooth at the cutoff)
df$z <- rnorm(n)  # a covariate
rdd_balance_test(y ~ x, data = df, cutoff = 0, covariates = "z")
```

### Plotting

```r
# Requires ggplot2
plot(m)  # scatterplot with loess fits on both sides of cutoff
```

---

## 4. Instrumental Variables / GMM

`ivgmm_flex()` estimates models with endogenous regressors using 6
methods, with automatic diagnostic tests.

### Formula Syntax

```
outcome ~ endogenous_regressors + exogenous_regressors | all_instruments + exogenous_regressors
```

The `|` separates regressors (LHS) from instruments (RHS).

### Full Signature

```r
ivgmm_flex(formula, data,
  method         = "2sls",            # "2sls" | "liml" | "fuller" | "gmm_twostep"
                                      # | "gmm_iterative" | "cue"
  vcov           = "classical",       # "classical" | "robust" | "HAC" | "cluster" | "bootstrap"
  weights_matrix = "robust",          # GMM weight matrix
  HAC_kernel     = "bartlett",        # "bartlett" | "parzen" | "qs"
  HAC_bw         = "auto",            # HAC bandwidth (numeric or "auto")
  fuller_alpha   = 1,                 # Fuller modification parameter
  cluster        = NULL,              # cluster variable name (for cluster vcov)
  maxiter        = 100,               # max iterations (iterative GMM)
  tol            = 1e-8,              # convergence tolerance
  latex          = FALSE
)
```

### Example: 2SLS

```r
set.seed(42)
n <- 500
z1 <- rnorm(n)
z2 <- rnorm(n)
x_endog <- 0.5 * z1 + 0.3 * z2 + rnorm(n)   # endogenous regressor
y <- 1 + 2 * x_endog + rnorm(n)
df <- data.frame(y = y, x = x_endog, z1 = z1, z2 = z2)

m <- ivgmm_flex(y ~ x | z1 + z2, data = df, method = "2sls")
summary(m)
```

### Example: LIML with Robust Standard Errors

```r
m2 <- ivgmm_flex(y ~ x | z1 + z2, data = df,
                 method = "liml", vcov = "robust")
summary(m2)
```

### Diagnostics (automatic)

```r
m$diagnostics$weak_iv_test        # Weak instrument F-statistics
m$diagnostics$overid_test         # Sargan-Hansen over-identification test
m$diagnostics$endogeneity_test    # Durbin-Wu-Hausman test
m$diagnostics$first_stage_F       # First-stage F-statistics per regressor
```

---

## 5. Difference-in-Differences

`did_flex()` implements TWFE plus 6 modern heterogeneity-robust DiD estimators
for staggered adoption designs.

### Full Signature

```r
did_flex(formula, data, id_var, time_var, treat_var,
  cohort_var    = NULL,                # auto-computed if NULL
  estimator     = "callaway_santanna", # see below for all 7 options
  control_group = "nevertreated",      # "nevertreated" | "notyettreated"
  anticipation  = 0,                   # number of anticipation periods
  aggregation   = "dynamic",           # "dynamic" | "simple" | "group" | "calendar"
  cluster       = NULL,                # cluster variable name
  boot_R        = 999,                 # bootstrap replications
  latex         = FALSE
)
```

**Available estimators:**
| Estimator | Reference |
|:----------|:----------|
| `"twfe"` | Classic two-way fixed effects |
| `"callaway_santanna"` | Callaway & Sant'Anna (2021) |
| `"sun_abraham"` | Sun & Abraham (2021) |
| `"borusyak_jaravel_spiess"` | Borusyak, Jaravel & Spiess (2024) |
| `"dechaisemartin_dhaultfoeuille"` | de Chaisemartin & d'Haultfœuille (2020) |
| `"gardner"` | Gardner (2022) |
| `"stacked"` | Stacked regression |

### Example: Callaway & Sant'Anna

```r
set.seed(42)
n_units <- 50; n_time <- 10
panel <- expand.grid(id = 1:n_units, year = 2001:2010)
panel$treated <- as.integer(panel$id > 30 & panel$year >= 2006)
panel$y <- 1 + 0.5 * panel$treated + rnorm(nrow(panel))

m <- did_flex(y ~ 1, data = panel, id_var = "id", time_var = "year",
             treat_var = "treated", estimator = "callaway_santanna")
summary(m)
m$att      # average treatment effect on the treated
m$att_se   # standard error
```

### Example: TWFE

```r
m_twfe <- did_flex(y ~ 1, data = panel, id_var = "id", time_var = "year",
                   treat_var = "treated", estimator = "twfe")
summary(m_twfe)
```

### Compare Estimators

```r
did_compare(
  did_flex(y ~ 1, data = panel, id_var = "id", time_var = "year",
           treat_var = "treated", estimator = "twfe"),
  did_flex(y ~ 1, data = panel, id_var = "id", time_var = "year",
           treat_var = "treated", estimator = "callaway_santanna"),
  labels = c("TWFE", "CS")
)
```

### Event Study Plot

```r
# Requires ggplot2
plot(m)  # plots event study coefficients with confidence intervals
```

---

## 6. Tobit Models

`tobit_flex()` estimates censored regression models of Types I through V,
with flexible error distributions and Powell's censored LAD estimator.

### Formula Syntax

- **Type I**: `y ~ x1 + x2` (standard censored)
- **Type II**: `y | selection ~ x1 | z1` (Heckman-like)
- **Types III–V**: Multi-part formulas with `|` separators

### Full Signature

```r
tobit_flex(formula, data,
  tobit_type      = 1,                 # integer 1-5
  left            = 0,                 # left censoring point
  right           = Inf,               # right censoring point (Inf = no right censoring)
  error_dist      = "normal",          # "normal" | "t" | "logistic"
  heteroscedastic = NULL,              # formula for heteroscedastic errors
  method          = "ml",             # "ml" | "twostep" | "powell" (censored LAD)
  latex           = FALSE
)
```

### Example: Type I Tobit (left-censored at 0)

```r
set.seed(42)
n <- 500
x <- rnorm(n)
y_star <- 2 + 1.5 * x + rnorm(n)
y <- pmax(y_star, 0)   # censored at 0
df <- data.frame(y = y, x = x)

m <- tobit_flex(y ~ x, data = df, tobit_type = 1, left = 0)
summary(m)
```

### Example: Type I with t-distributed Errors

```r
m2 <- tobit_flex(y ~ x, data = df, tobit_type = 1, left = 0, error_dist = "t")
summary(m2)
```

### Example: Powell's Censored LAD (robust to non-normality)

```r
m3 <- tobit_flex(y ~ x, data = df, tobit_type = 1, left = 0, method = "powell")
summary(m3)
```

### Example: Double Censoring (left = 0, right = 5)

```r
y_double <- pmin(pmax(y_star, 0), 5)
df2 <- data.frame(y = y_double, x = x)
m4 <- tobit_flex(y ~ x, data = df2, tobit_type = 1, left = 0, right = 5)
summary(m4)
```

---

## 7. Nonlinear Panel Data

`panel_nl()` estimates nonlinear panel models with fixed effects, random effects,
or Chamberlain's correlated random effects.

### Full Signature

```r
panel_nl(formula, data, id, time,
  family       = "poisson",           # "poisson" | "negbin" | "logit" | "probit" | "gamma"
  effect       = "individual",        # "individual" | "time" | "twoways"
  model        = "fe",                # "fe" | "re" | "correlated_re"
  vcov         = "standard",          # "standard" | "robust" | "cluster"
                                      # | "driscoll_kraay" | "bootstrap"
  cluster      = NULL,                # cluster variable name
  DK_bandwidth = NULL,                # Driscoll-Kraay bandwidth
  latex        = FALSE
)
```

### Example: Poisson Fixed Effects

```r
set.seed(42)
n_firms <- 50; n_years <- 8
panel <- expand.grid(firm = 1:n_firms, year = 1:n_years)
panel$x <- rnorm(nrow(panel))
panel$firm_fe <- rep(rnorm(n_firms, 0, 0.5), each = n_years)
panel$y <- rpois(nrow(panel), exp(0.5 + 0.3 * panel$x + panel$firm_fe))

m <- panel_nl(y ~ x, data = panel, id = "firm", time = "year",
              family = "poisson", model = "fe")
summary(m)
```

### Example: Logit Random Effects with Cluster SE

```r
panel$y_bin <- as.integer(panel$y > 1)
m2 <- panel_nl(y_bin ~ x, data = panel, id = "firm", time = "year",
               family = "logit", model = "re",
               vcov = "cluster", cluster = "firm")
summary(m2)
```

---

## 8. Switching Regression

`switching_flex()` estimates models where observations come from different
regimes. Supports exogenous switching (known regimes), endogenous switching
(Maddala 1983, à la Heckman), and Markov switching (Hamilton 1989).

### Full Signature

```r
switching_flex(formula, data,
  type                  = "endogenous",   # "endogenous" | "exogenous" | "markov"
  n_regimes             = 2,              # number of regimes (markov)
  regime_var            = NULL,           # column name with regime labels (exogenous)
  error_dist            = "normal",       # "normal" | "t"
  transition            = "constant",     # "constant" | "time_varying" (markov)
  transition_covariates = NULL,           # covariates for transition probs
  latex                 = FALSE
)
```

### Example: Exogenous Switching

```r
set.seed(42)
n <- 300
regime <- sample(1:2, n, replace = TRUE)
x <- rnorm(n)
y <- ifelse(regime == 1, 3 + x + rnorm(n), 1 + 0.5*x + rnorm(n))
df <- data.frame(y = y, x = x, regime = regime)

m <- switching_flex(y ~ x, data = df, type = "exogenous", regime_var = "regime")
summary(m)
```

### Example: Endogenous Switching (Maddala)

Uses a selection equation to model regime assignment jointly:

```r
set.seed(42)
n <- 500
x <- rnorm(n); z <- rnorm(n)
sel <- as.integer(0.3 * x + 0.5 * z + rnorm(n) > 0)
y <- ifelse(sel == 1, 3 + x + rnorm(n), 1 + 0.5*x + rnorm(n))
df <- data.frame(y = y, x = x, sel = sel, z = z)

# Formula: outcome | selection ~ outcome_covs | selection_covs
m <- switching_flex(y | sel ~ x | z, data = df, type = "endogenous")
summary(m)
```

### Example: Markov Switching (Hamilton)

```r
set.seed(42)
n <- 300
state <- numeric(n); state[1] <- 1
for (t in 2:n) {
  state[t] <- sample(1:2, 1,
    prob = if (state[t-1] == 1) c(0.9, 0.1) else c(0.2, 0.8))
}
y <- ifelse(state == 1, 2 + rnorm(n, sd = 0.5), -1 + rnorm(n, sd = 1))
df <- data.frame(y = y, const = 1)

m <- switching_flex(y ~ const, data = df, type = "markov", n_regimes = 2)
summary(m)

# Filtered regime probabilities (n x n_regimes matrix):
head(m$filtered_probs)
```

---

## 9. Flexible Standard Errors

`flex_vcov()` recomputes the variance-covariance matrix for any ecoflex
model *after* estimation.

```r
flex_vcov(object,
  type     = "sandwich",               # "hessian" | "opg" | "sandwich" | "cluster" | "bootstrap"
  cluster  = NULL,                     # vector of cluster IDs
  R        = 500,                      # bootstrap replications
  parallel = FALSE,                    # parallel bootstrap
  cores    = 2                         # number of cores
)
```

### Examples

```r
# From any model:
m <- hurdle_flex(y ~ x, data = df)

# Sandwich (robust) standard errors
V_sandwich <- flex_vcov(m, type = "sandwich")

# Cluster-robust
V_cluster <- flex_vcov(m, type = "cluster", cluster = df$group_id)

# Bootstrap
V_boot <- flex_vcov(m, type = "bootstrap", R = 1000)

# Use the new vcov to get corrected SEs:
se_robust <- sqrt(diag(V_sandwich))
```

---

## 10. LaTeX Output

`to_latex()` generates publication-ready LaTeX tables. For two-stage
models (Heckman, Tobit II+, Switching), it displays the selection and
outcome equations side by side.

```r
# Print to console
to_latex(m)

# Write to file
to_latex(m, file = "results_table.tex")

# With custom caption
to_latex(m, caption = "Hurdle Model Estimates")

# Or generate LaTeX directly when estimating:
m <- hurdle_flex(y ~ x, data = df, latex = TRUE)
```

---

## Common Methods (All Models)

Every ecoflex model object supports these standard R generics:

```r
summary(m)       # Coefficient table: estimate, SE, z-value, p-value
coef(m)          # Named coefficient vector
vcov(m)          # Variance-covariance matrix
logLik(m)        # Log-likelihood value
nobs(m)          # Number of observations
predict(m)       # Fitted values
predict(m, newdata = new_df)   # Out-of-sample prediction

# For DiD and RDD models:
plot(m)          # Visual diagnostics (event study / RD plot)
```
