# ecoflex <img src="man/figures/logo.png" align="right" height="139" />

[![R-CMD-check](https://github.com/zabbonat/ecoflex/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main1)](https://github.com/zabbonat/ecoflex/actions/workflows/R-CMD-check.yaml)

**ecoflex** provides a unified, flexible R interface for advanced
econometric models commonly used in applied microeconomics and causal
inference research. Every model supports **flexible standard errors**,
**LaTeX output**, and a consistent `summary()` / `coef()` / `predict()` API.

## Installation

```r
# install.packages("remotes")
remotes::install_github("zabbonat/ecoflex", ref = "main1")
```

> 📖 **[Full Documentation & Running Examples →](DOCUMENTATION.md)**
>
> Complete guide with runnable code for every function.

---

## Model Families

### 1. Hurdle Models — `hurdle_flex()`

Count data models with an arbitrary hurdle threshold.

```r
hurdle_flex(formula, data,
  threshold = 0,                             # any non-negative integer
  dist      = c("poisson", "negbin", "geometric"),
  zero_dist = c("binomial", "probit"),       # binary component link
  method    = c("ml", "two-step"),
  start     = NULL,                          # optional starting values
  latex     = FALSE                          # print LaTeX table
)
```

**Example:**
```r
m <- hurdle_flex(y ~ x1 + x2, data = df, threshold = 2, dist = "negbin")
summary(m)
to_latex(m)
```

---

### 2. Heckman Selection Models — `heckman_flex()`

Sample selection correction with multiple estimation methods.

```r
heckman_flex(outcome | selection ~ outcome_covs | selection_covs, data,
  method         = c("ml", "twostep", "twostep_robust", "semipar"),
  selection_type = c("binary", "ordered", "multinomial"),
  error_dist     = c("normal", "t", "copula"),
  copula_family  = c("gaussian", "frank", "clayton", "gumbel"),
  heteroscedastic = NULL,                    # formula for heteroscedasticity
  vcov           = c("standard", "robust", "bootstrap"),
  latex          = FALSE
)
```

**Example:**
```r
m <- heckman_flex(wage | employed ~ education + experience | education + kids,
                  data = df, method = "twostep")
summary(m)
```

---

### 3. Regression Discontinuity — `rdd_flex()`

Sharp and fuzzy RDD with MSE-optimal bandwidth selection.

```r
rdd_flex(formula, data,
  cutoff          = 0,
  type            = c("sharp", "fuzzy"),
  fuzzy_treatment = NULL,                   # variable name for fuzzy treatment
  bandwidth       = "mserd",                # or numeric value
  kernel          = c("triangular", "epanechnikov", "uniform", "gaussian"),
  polynomial      = 1,                      # polynomial order
  polynomial_right = NULL,                  # separate order for right of cutoff
  discrete_running = FALSE,
  cluster         = NULL,
  covs_adjust     = c("none", "linear", "interacted"),
  masspoints      = c("adjust", "check", "off"),
  all_bandwidths  = FALSE,
  latex           = FALSE
)
```

**Diagnostic functions:**
```r
rdd_manipulation_test(x, cutoff = 0, method = c("cattaneo", "mccrary"))
rdd_balance_test(formula, data, cutoff = 0, covariates = c("x1", "x2"))
rdd_placebo_cutoffs(formula, data, cutoff = 0, placebo_cutoffs = NULL)
```

**Example:**
```r
m <- rdd_flex(y ~ running_var, data = df, cutoff = 0, kernel = "epanechnikov")
summary(m)
plot(m)

# Diagnostics
rdd_manipulation_test(df$running_var, cutoff = 0)
```

---

### 4. Instrumental Variables / GMM — `ivgmm_flex()`

IV estimation with 6 methods and full diagnostics.

```r
ivgmm_flex(y ~ endogenous_regressors | instruments, data,
  method        = c("2sls", "liml", "fuller", "gmm_twostep",
                    "gmm_iterative", "cue"),
  vcov          = c("classical", "robust", "HAC", "cluster", "bootstrap"),
  weights_matrix = "robust",
  HAC_kernel    = c("bartlett", "parzen", "qs"),
  HAC_bw        = "auto",
  fuller_alpha  = 1,                        # Fuller modification parameter
  cluster       = NULL,
  maxiter       = 100,                      # for iterative GMM
  tol           = 1e-8,
  latex         = FALSE
)
```

**Built-in diagnostics** (automatic):
- Sargan-Hansen over-identification test
- Weak instrument F-test
- Durbin-Wu-Hausman endogeneity test
- First-stage diagnostics

**Example:**
```r
m <- ivgmm_flex(y ~ x_endog + x_exog | z1 + z2 + x_exog, data = df,
                method = "2sls", vcov = "robust")
summary(m)
m$diagnostics  # all test results
```

---

### 5. Difference-in-Differences — `did_flex()`

TWFE plus 6 modern heterogeneity-robust estimators.

```r
did_flex(formula, data, id_var, time_var, treat_var,
  cohort_var    = NULL,                     # auto-computed if NULL
  estimator     = c("callaway_santanna", "sun_abraham",
                    "borusyak_jaravel_spiess",
                    "dechaisemartin_dhaultfoeuille",
                    "gardner", "stacked", "twfe"),
  control_group = c("nevertreated", "notyettreated"),
  anticipation  = 0,
  aggregation   = c("dynamic", "simple", "group", "calendar"),
  cluster       = NULL,
  boot_R        = 999,
  latex         = FALSE
)
```

**Compare estimators:**
```r
did_compare(
  did_flex(..., estimator = "twfe"),
  did_flex(..., estimator = "callaway_santanna"),
  labels = c("TWFE", "CS")
)
```

**Example:**
```r
m <- did_flex(y ~ 1, data = panel, id_var = "id", time_var = "year",
             treat_var = "treated", estimator = "callaway_santanna")
summary(m)
plot(m)  # event study plot
```

---

### 6. Tobit Models (Types I–V) — `tobit_flex()`

Censored regression with flexible error distributions.

```r
tobit_flex(formula, data,
  tobit_type = 1,                           # integer 1-5
  left       = 0,                           # left censoring point
  right      = Inf,                         # right censoring point
  error_dist = c("normal", "t", "logistic"),
  heteroscedastic = NULL,                   # formula for heteroscedasticity
  method     = c("ml", "twostep", "powell"),# Powell = censored LAD
  latex      = FALSE
)
```

For types II–V, use multi-part formulas with `|`:
```r
# Type II (Heckman-like): outcome | selection ~ covs_outcome | covs_selection
# Type V: y1 | y2 | selection ~ x1 | x2 | z
```

**Example:**
```r
m <- tobit_flex(y ~ x1 + x2, data = df, tobit_type = 1,
                left = 0, error_dist = "t")
summary(m)
```

---

### 7. Nonlinear Panel Data — `panel_nl()`

Fixed, random, and correlated random effects for nonlinear models.

```r
panel_nl(formula, data, id, time,
  family       = c("poisson", "negbin", "logit", "probit", "gamma"),
  effect       = c("individual", "time", "twoways"),
  model        = c("fe", "re", "correlated_re"),
  vcov         = c("standard", "robust", "cluster", "driscoll_kraay",
                   "bootstrap"),
  cluster      = NULL,
  DK_bandwidth = NULL,                      # for Driscoll-Kraay
  latex        = FALSE
)
```

**Example:**
```r
m <- panel_nl(y ~ x1 + x2, data = panel_df, id = "firm_id", time = "year",
              family = "poisson", model = "fe", vcov = "cluster",
              cluster = "firm_id")
summary(m)
```

---

### 8. Switching Regression — `switching_flex()`

Exogenous, endogenous (Maddala), and Markov regime-switching models.

```r
switching_flex(formula, data,
  type                  = c("endogenous", "exogenous", "markov"),
  n_regimes             = 2,
  regime_var            = NULL,              # for exogenous switching
  error_dist            = c("normal", "t"),
  transition            = c("constant", "time_varying"),
  transition_covariates = NULL,              # covariates for transition probs
  latex                 = FALSE
)
```

**Formulas by type:**
```r
# Exogenous:  y ~ x (uses regime_var to split)
# Endogenous: y | selection ~ x | z  (à la Maddala)
# Markov:     y ~ x (hidden regimes estimated via EM)
```

**Example:**
```r
m <- switching_flex(y | sel ~ x | z, data = df, type = "endogenous")
summary(m)

m <- switching_flex(y ~ x, data = df, type = "markov", n_regimes = 2)
m$filtered_probs  # regime probabilities
```

---

## Cross-Cutting Features

### Flexible Standard Errors — `flex_vcov()`

Post-estimation variance-covariance recomputation for any model.

```r
flex_vcov(object,
  type     = c("hessian", "opg", "sandwich", "cluster", "bootstrap"),
  cluster  = NULL,
  R        = 500,          # bootstrap replications
  parallel = FALSE,
  cores    = 2
)
```

```r
V <- flex_vcov(m, type = "cluster", cluster = df$firm_id)
```

### LaTeX Output — `to_latex()`

Publication-ready tables for any model. Two-stage models display
selection and outcome equations side by side.

```r
to_latex(m)                          # print to console
to_latex(m, file = "table.tex")      # write to file
to_latex(m, caption = "My results")
```

Or pass `latex = TRUE` directly when estimating:
```r
m <- hurdle_flex(y ~ x, data = df, latex = TRUE)  # prints table on creation
```

### Common Methods

All ecoflex objects support:

| Method | Description |
|:---|:---|
| `summary(m)` | Coefficient table with standard errors, z-values, p-values |
| `coef(m)` | Extract coefficients |
| `vcov(m)` | Variance-covariance matrix |
| `logLik(m)` | Log-likelihood |
| `nobs(m)` | Number of observations |
| `predict(m)` | Fitted values / predictions |
| `plot(m)` | Visual diagnostics (for DiD and RDD) |

---

## License

MIT © ecoflex authors
