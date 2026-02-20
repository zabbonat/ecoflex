# ecoflex <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/PLACEHOLDER/ecoflex/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PLACEHOLDER/ecoflex/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**ecoflex** provides a unified, flexible R interface for advanced
econometric models commonly used in applied microeconomics and causal
inference research.

## Installation

```r
# Install from GitHub
# devtools::install_github("PLACEHOLDER/ecoflex")
```

## Model Families

| Module | Function | Key Features |
|:---|:---|:---|
| **Hurdle** | `hurdle_flex()` | Arbitrary threshold, ML/two-step, Poisson/NegBin/Geometric |
| **Heckman** | `heckman_flex()` | ML/two-step/Lee/semiparametric, normal/t errors |
| **RDD** | `rdd_flex()` | Sharp/fuzzy, local polynomial, MSE-optimal bandwidth |
| **IV/GMM** | `ivgmm_flex()` | 2SLS/LIML/Fuller/GMM/CUE |
| **DiD** | `did_flex()` | TWFE + 6 modern heterogeneity-robust estimators |
| **Tobit** | `tobit_flex()` | Types I–V, normal/t/logistic, censored LAD |
| **Panel** | `panel_nl()` | FE/RE/CRE, Poisson/NegBin/Logit/Probit/Gamma |
| **Switching** | `switching_flex()` | Exogenous/endogenous/Markov switching |

## Quick Start

```r
library(ecoflex)

# Hurdle model
m <- hurdle_flex(y ~ x, data = df, threshold = 0)
summary(m)

# Sharp RDD
m <- rdd_flex(y ~ running_var, data = df, cutoff = 0)
summary(m)

# IV/GMM
m <- ivgmm_flex(y ~ endog | instruments, data = df, method = "2sls")
summary(m)
```

## Flexible Standard Errors

All models support `flex_vcov()` for robust/cluster/bootstrap SE:

```r
V <- flex_vcov(m, type = "cluster", cluster = df$firm_id)
```

## License

MIT © ecoflex authors
