# ecoflex 0.1.0

* Initial CRAN release.
* Eight econometric model families:
  - `hurdle_flex()`: Hurdle models with arbitrary thresholds (Poisson, NegBin, Geometric).
  - `heckman_flex()`: Extended Heckman selection (ML, two-step, Lee, semiparametric).
  - `rdd_flex()`: Sharp and fuzzy RDD with local polynomial estimation.
  - `ivgmm_flex()`: Unified IV/GMM (2SLS, LIML, Fuller, GMM, CUE).
  - `did_flex()`: Difference-in-differences with modern heterogeneity-robust estimators.
  - `tobit_flex()`: Tobit types I-V with flexible error distributions.
  - `panel_nl()`: Nonlinear panel models (Poisson, NegBin, Logit, Probit, Gamma).
  - `switching_flex()`: Exogenous, endogenous, and Markov switching regressions.
* Unified S3 interface: `summary()`, `coef()`, `vcov()`, `predict()`, `logLik()`, `nobs()`.
* LaTeX table output via `to_latex()` and `summary(m, latex = TRUE)`.
* Post-estimation variance: robust, cluster, bootstrap, HAC via `flex_vcov()`.
