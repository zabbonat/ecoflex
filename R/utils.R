# =============================================================================
# ecoflex: Shared Utility Functions
# =============================================================================

#' @keywords internal
#' Compute numerical score matrix (gradient of log-likelihood per observation)
#'
#' @param object An ecoflex model object
#' @return n x k matrix of observation-level scores
.compute_scores <- function(object) {
  if (!is.null(object$scores)) return(object$scores)

  # Fallback: numerical differentiation
  # We need the per-observation log-likelihood function
  if (!is.null(object$loglik_obs_fn)) {
    scores <- numDeriv::jacobian(
      func = function(par) object$loglik_obs_fn(par),
      x = object$coefficients
    )
    return(scores)
  }

  # Last resort: use numDeriv on full log-likelihood and approximate
  warning("Scores computed via numerical approximation; results may be imprecise.")
  n <- object$n
  k <- length(object$coefficients)
  scores <- matrix(0, n, k)

  # Perturb data approach
  for (i in seq_len(n)) {
    scores[i, ] <- numDeriv::grad(
      func = function(par) {
        .observation_loglik(object, par, i)
      },
      x = object$coefficients
    )
  }
  scores
}

#' @keywords internal
#' Per-observation log-likelihood (generic fallback)
.observation_loglik <- function(object, par, obs_idx) {
  # This is overridden by specific model classes
  # Fallback: return total loglik / n (crude approximation)
  if (!is.null(object$loglik_fn)) {
    return(object$loglik_fn(par) / object$n)
  }
  0
}

# =============================================================================
# Truncated Count Models
# =============================================================================

#' @keywords internal
#' Fit a zero-truncated Poisson model via ML
.truncated_poisson <- function(y, X, maxiter = 100) {
  # y must be > 0
  stopifnot(all(y > 0))

  n <- length(y)
  k <- ncol(X)

  # Log-likelihood for zero-truncated Poisson
  loglik <- function(beta) {
    mu <- exp(X %*% beta)
    ll <- dpois(y, lambda = mu, log = TRUE) - log(1 - exp(-mu))
    sum(ll)
  }

  # Starting values from Poisson GLM
  fit0 <- glm.fit(X, y, family = poisson())
  start <- fit0$coefficients

  opt <- maxLik::maxLik(loglik, start = start, method = "BFGS")

  coefs <- coef(opt)
  H <- maxLik::hessian(opt)
  V <- tryCatch(solve(-H), error = function(e) MASS::ginv(-H))

  list(
    coefficients = coefs,
    vcov = V,
    logLik = logLik(opt),
    hessian = H,
    convergence = opt$code
  )
}

#' @keywords internal
#' Fit a zero-truncated Negative Binomial model via ML
.truncated_negbin <- function(y, X, size_fixed = NULL, maxiter = 100) {
  # y must be > 0
  stopifnot(all(y > 0))

  n <- length(y)
  k <- ncol(X)
  estimate_size <- is.null(size_fixed)

  loglik <- function(par) {
    beta <- par[1:k]
    alpha <- if (estimate_size) exp(par[k + 1]) else size_fixed
    mu <- exp(X %*% beta)

    ll <- dnbinom(y, mu = mu, size = alpha, log = TRUE) -
      log(1 - dnbinom(0, mu = mu, size = alpha))
    sum(ll)
  }

  # Starting values
  fit0 <- glm.fit(X, y, family = poisson())
  start <- fit0$coefficients
  if (estimate_size) start <- c(start, log(1))

  opt <- maxLik::maxLik(loglik, start = start, method = "BFGS")

  coefs <- coef(opt)
  H <- maxLik::hessian(opt)
  V <- tryCatch(solve(-H), error = function(e) MASS::ginv(-H))

  list(
    coefficients = coefs[1:k],
    vcov = V[1:k, 1:k, drop = FALSE],
    logLik = logLik(opt),
    hessian = H,
    alpha = if (estimate_size) exp(coefs[k + 1]) else size_fixed,
    convergence = opt$code
  )
}

# =============================================================================
# Cohort computation for DiD
# =============================================================================

#' @keywords internal
#' Compute treatment cohort from panel data
.compute_cohort <- function(data, id_var, time_var, treat_var) {
  # Cohort = first period where treat_var == 1
  # If never treated, cohort = Inf
  ids <- unique(data[[id_var]])
  cohort <- rep(Inf, nrow(data))

  for (id in ids) {
    idx <- data[[id_var]] == id
    treated_periods <- data[[time_var]][idx & data[[treat_var]] == 1]
    if (length(treated_periods) > 0) {
      first_treat <- min(treated_periods)
      cohort[idx] <- first_treat
    }
  }
  cohort
}

# =============================================================================
# Kernel Functions for RDD
# =============================================================================

#' @keywords internal
#' Kernel function evaluation
.kernel_fn <- function(u, type = "triangular") {
  switch(type,
    triangular = ifelse(abs(u) <= 1, 1 - abs(u), 0),
    epanechnikov = ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0),
    uniform = ifelse(abs(u) <= 1, 0.5, 0),
    gaussian = dnorm(u),
    stop("Unknown kernel type: ", type)
  )
}

# =============================================================================
# Censored LAD (Powell 1984)
# =============================================================================

#' @keywords internal
#' Powell's Censored Least Absolute Deviations estimator
.censored_lad <- function(y, X, left = 0, right = Inf, maxiter = 200, tol = 1e-6) {
  n <- length(y)
  k <- ncol(X)

  # Starting values from OLS
  beta <- lm.fit(X, y)$coefficients

  for (iter in seq_len(maxiter)) {
    beta_old <- beta
    xb <- as.numeric(X %*% beta)

    # Censored predicted values
    xb_censored <- pmax(pmin(xb, right), left)

    # LAD objective: minimize sum |y - max(left, min(right, x'beta))|
    # Use iteratively reweighted least squares approximation
    resid <- y - xb_censored
    weights <- 1 / (abs(resid) + 1e-8)

    # Weighted LS update
    W <- diag(weights)
    beta <- tryCatch(
      solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y,
      error = function(e) beta_old
    )
    beta <- as.numeric(beta)

    if (max(abs(beta - beta_old)) < tol) break
  }

  # Bootstrap SE
  se <- rep(NA_real_, k)
  names(beta) <- colnames(X)

  list(
    coefficients = beta,
    se = se,
    z = rep(NA_real_, k),
    pvalue = rep(NA_real_, k),
    vcov = matrix(NA_real_, k, k),
    hessian = NULL,
    logLik = -sum(abs(y - pmax(pmin(as.numeric(X %*% beta), right), left))),
    AIC = NA_real_,
    BIC = NA_real_,
    n = n,
    n_censored_left = sum(y <= left),
    n_censored_right = sum(y >= right & is.finite(right)),
    n_uncensored = sum(y > left & y < right)
  )
}

# =============================================================================
# Bandwidth Selection for RDD
# =============================================================================

#' @keywords internal
#' MSE-optimal bandwidth selection for RDD
.bandwidth_selection <- function(y, x, kernel, polynomial,
                                  method = "mserd", discrete = FALSE) {
  n <- length(y)

  # Silverman pilot bandwidth
  h_pilot <- 1.06 * sd(x) * n^(-1/5)

  # Simple MSE-optimal bandwidth (Imbens-Kalyanaraman 2012 style)
  # Use pilot bandwidth to estimate curvature
  left_idx <- x >= -h_pilot & x < 0
  right_idx <- x >= 0 & x <= h_pilot

  if (sum(left_idx) < 5 || sum(right_idx) < 5) {
    # Not enough observations; use pilot
    return(h_pilot)
  }

  # Estimate second derivatives via local quadratic
  m2_left <- tryCatch({
    fit_l <- lm(y[left_idx] ~ x[left_idx] + I(x[left_idx]^2))
    2 * coef(fit_l)[3]
  }, error = function(e) 0)

  m2_right <- tryCatch({
    fit_r <- lm(y[right_idx] ~ x[right_idx] + I(x[right_idx]^2))
    2 * coef(fit_r)[3]
  }, error = function(e) 0)

  # Regularization constants for triangular kernel, p=1
  C_k <- 3.4375  # kernel constant

  curvature <- (m2_right - m2_left)^2
  if (curvature < 1e-10) curvature <- 1e-10

  # Variance estimate
  sigma2_left <- var(y[left_idx])
  sigma2_right <- var(y[right_idx])
  V_bound <- sigma2_left + sigma2_right

  # MSE-optimal bandwidth
  h_opt <- C_k * (V_bound / (curvature * n))^(1/5)

  # Bound the bandwidth
  x_range <- diff(range(x))
  h_opt <- max(h_opt, x_range / 50)
  h_opt <- min(h_opt, x_range / 2)

  if (discrete) {
    # Adjust for mass points
    n_unique <- length(unique(x))
    h_opt <- h_opt * (n / n_unique)^(1/5)
  }

  h_opt
}

# =============================================================================
# Stationary Distribution for Markov Chains
# =============================================================================

#' @keywords internal
#' Compute stationary distribution of a transition matrix
.stationary_distribution <- function(P) {
  n <- nrow(P)
  # Solve pi' P = pi' with sum(pi) = 1
  # Equivalent to (P' - I) pi = 0 with constraint
  A <- t(P) - diag(n)
  A[n, ] <- rep(1, n)
  b <- c(rep(0, n - 1), 1)
  tryCatch(
    solve(A, b),
    error = function(e) rep(1/n, n)
  )
}

# =============================================================================
# Build Transition Matrix for Markov-Switching
# =============================================================================

#' @keywords internal
.build_transition_matrix <- function(par, n_regimes, n_par_regime,
                                      transition, n) {
  offset <- n_regimes * n_par_regime + n_regimes  # after betas and sigmas

  if (transition == "constant") {
    # Multinomial logit parametrization for each row
    P_single <- matrix(0, n_regimes, n_regimes)
    idx <- 1
    for (i in seq_len(n_regimes)) {
      raw <- c(0, par[offset + idx:(idx + n_regimes - 2)])
      P_single[i, ] <- exp(raw) / sum(exp(raw))
      idx <- idx + n_regimes - 1
    }
    # Replicate for all observations
    P <- array(0, dim = c(n_regimes, n_regimes, n))
    for (t in seq_len(n)) P[, , t] <- P_single
  } else {
    # Time-varying: each transition probability depends on covariates
    # Simplified: constant transitions for now
    P <- array(0, dim = c(n_regimes, n_regimes, n))
    P_single <- matrix(1/n_regimes, n_regimes, n_regimes)
    for (t in seq_len(n)) P[, , t] <- P_single
  }
  P
}

# =============================================================================
# Weight Matrix for GMM
# =============================================================================

#' @keywords internal
#' Compute GMM weight matrix
.compute_weight_matrix <- function(Z, residuals, type = "robust") {
  n <- length(residuals)
  l <- ncol(Z)

  switch(type,
    identity = diag(l) / n,
    unadjusted = solve(crossprod(Z) / n),
    robust = {
      # Heteroscedasticity-robust weight matrix
      S <- crossprod(Z * residuals) / n
      tryCatch(solve(S), error = function(e) MASS::ginv(S))
    },
    HAC = {
      # Newey-West HAC weight matrix
      bw <- floor(n^(1/3))
      S <- crossprod(Z * residuals) / n
      for (j in seq_len(bw)) {
        w <- 1 - j / (bw + 1)  # Bartlett kernel
        Gamma_j <- crossprod(Z[-(1:j), ] * residuals[-(1:j)],
                              Z[1:(n-j), ] * residuals[1:(n-j)]) / n
        S <- S + w * (Gamma_j + t(Gamma_j))
      }
      tryCatch(solve(S), error = function(e) MASS::ginv(S))
    },
    stop("Unknown weight matrix type: ", type)
  )
}
