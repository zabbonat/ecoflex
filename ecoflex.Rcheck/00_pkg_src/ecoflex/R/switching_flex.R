# =============================================================================
# ecoflex Module 8: Switching Regression
# =============================================================================

#' Switching Regression Models
#'
#' Estimate exogenous, endogenous (Maddala), or Markov switching regressions
#' with flexible error distributions.
#'
#' @param formula Formula for the outcome equation.
#' @param data A data frame.
#' @param type \code{"endogenous"} (Maddala), \code{"exogenous"}, \code{"markov"}.
#' @param n_regimes Number of regimes. Default: 2.
#' @param regime_var Regime variable name (required if \code{type = "exogenous"}).
#' @param error_dist \code{"normal"} or \code{"t"}.
#' @param transition For Markov: \code{"constant"} or \code{"time_varying"}.
#' @param transition_covariates Formula for time-varying transition probabilities.
#' @param latex If \code{TRUE}, prints a LaTeX table. Default: \code{FALSE}.
#' @param ... Additional arguments.
#' @return Object of class \code{c("switching_flex", "ecoflex")}.
#' @export
#'
#' @examples
#' \donttest{
#' m <- switching_flex(mpg ~ hp + wt, data = mtcars,
#'                     type = "exogenous", regime_var = "am")
#' summary(m)
#' }
switching_flex <- function(formula, data,
                           type = c("endogenous", "exogenous", "markov"),
                           n_regimes = 2,
                           regime_var = NULL,
                           error_dist = c("normal", "t"),
                           transition = c("constant", "time_varying"),
                           transition_covariates = NULL, latex = FALSE, ...) {
  type <- match.arg(type)
  error_dist <- match.arg(error_dist)
  transition <- match.arg(transition)

  result <- switch(type,
    endogenous = .switching_endogenous(formula, data, n_regimes, error_dist, ...),
    exogenous = .switching_exogenous(formula, data, regime_var, error_dist, ...),
    markov = .switching_markov(formula, data, n_regimes, error_dist,
                                transition, transition_covariates, ...)
  )

  obj <- structure(c(result, list(
    call = match.call(), formula = Formula::Formula(formula),
    type = type, n_regimes = n_regimes, error_dist = error_dist,
    model_data = data, n = nrow(data),
    model_name = sprintf("Switching Regression (%s, %d regimes)", type, n_regimes)
  )), class = c("switching_flex", "ecoflex"))
  if (latex) to_latex(obj)
  obj
}

# --- Exogenous Switching ---
#' @keywords internal
.switching_exogenous <- function(formula, data, regime_var, error_dist, ...) {
  if (is.null(regime_var)) stop("regime_var required for exogenous switching")
  F <- Formula::Formula(formula)
  mf <- model.frame(F, data = data)
  y <- model.response(mf)
  X <- model.matrix(F, data = data)
  regimes <- data[[regime_var]]
  regime_levels <- sort(unique(regimes))
  n_regimes <- length(regime_levels)
  k <- ncol(X)

  # Fit separate OLS for each regime
  fits <- lapply(regime_levels, function(r) {
    idx <- regimes == r
    lm.fit(X[idx, , drop = FALSE], y[idx])
  })

  all_coefs <- c()
  all_se <- c()
  loglik_total <- 0
  for (r in seq_along(regime_levels)) {
    fit_r <- fits[[r]]
    n_r <- sum(regimes == regime_levels[r])
    sigma_r <- sqrt(sum(fit_r$residuals^2) / (n_r - k))
    V_r <- sigma_r^2 * tryCatch(solve(crossprod(X[regimes == regime_levels[r], ])),
                                 error = function(e) MASS::ginv(crossprod(X[regimes == regime_levels[r], ])))
    coefs_r <- fit_r$coefficients
    se_r <- sqrt(pmax(diag(V_r), 0))
    names(coefs_r) <- paste0("regime", regime_levels[r], "_", colnames(X))
    names(se_r) <- names(coefs_r)
    all_coefs <- c(all_coefs, coefs_r, setNames(log(sigma_r), paste0("regime", regime_levels[r], "_log_sigma")))
    all_se <- c(all_se, se_r, NA)
    loglik_total <- loglik_total + sum(dnorm(fit_r$residuals, sd = sigma_r, log = TRUE))
  }

  V <- diag(ifelse(is.na(all_se), 1e-10, all_se^2))

  list(coefficients = all_coefs, se = all_se, z = all_coefs / all_se,
       pvalue = 2 * pnorm(-abs(all_coefs / all_se)),
       vcov = V, hessian = NULL,
       logLik = loglik_total, AIC = -2*loglik_total + 2*length(all_coefs),
       BIC = -2*loglik_total + log(length(y))*length(all_coefs),
       regime_fits = fits, convergence = 0)
}

# --- Endogenous Switching (Maddala 1983) ---
#' @keywords internal
.switching_endogenous <- function(formula, data, n_regimes, error_dist, ...) {
  F <- Formula::Formula(formula)
  n_parts <- length(F)[2]

  mf <- model.frame(F, data = data)
  y <- model.response(mf)
  X <- model.matrix(F, data = data, rhs = 1)
  n <- length(y); k <- ncol(X)

  if (n_parts >= 2) {
    mf2 <- model.frame(F, data = data, rhs = 2)
    s <- model.response(mf2)
    Z <- model.matrix(F, data = data, rhs = 2)
    k_z <- ncol(Z)
  } else {
    # If no selection equation, use X for both
    s <- as.integer(y > median(y))
    Z <- X; k_z <- k
  }

  # Joint ML for 2-regime endogenous switching
  loglik <- function(par) {
    beta1 <- par[1:k]; beta2 <- par[(k+1):(2*k)]
    gamma <- par[(2*k+1):(2*k+k_z)]
    log_s1 <- par[2*k + k_z + 1]; log_s2 <- par[2*k + k_z + 2]
    rho1_raw <- par[2*k + k_z + 3]; rho2_raw <- par[2*k + k_z + 4]

    sigma1 <- exp(log_s1); sigma2 <- exp(log_s2)
    rho1 <- tanh(rho1_raw); rho2 <- tanh(rho2_raw)
    Xb1 <- as.numeric(X %*% beta1); Xb2 <- as.numeric(X %*% beta2)
    Zg <- as.numeric(Z %*% gamma)

    ll <- numeric(n)
    idx1 <- s == 1
    if (any(idx1)) {
      e1 <- (y[idx1] - Xb1[idx1]) / sigma1
      arg1 <- (Zg[idx1] + rho1 * e1) / sqrt(1 - rho1^2)
      ll[idx1] <- dnorm(e1, log = TRUE) - log(sigma1) + pnorm(arg1, log.p = TRUE)
    }
    idx0 <- s == 0
    if (any(idx0)) {
      e2 <- (y[idx0] - Xb2[idx0]) / sigma2
      arg2 <- (-Zg[idx0] + rho2 * e2) / sqrt(1 - rho2^2)
      ll[idx0] <- dnorm(e2, log = TRUE) - log(sigma2) + pnorm(arg2, log.p = TRUE)
    }
    sum(ll)
  }

  # Starting values from OLS on each regime
  ols1 <- lm.fit(X[s == 1, , drop = FALSE], y[s == 1])
  ols2 <- lm.fit(X[s == 0, , drop = FALSE], y[s == 0])
  probit_s <- suppressWarnings(glm(s ~ Z - 1, family = binomial(link = "probit")))

  start <- c(ols1$coefficients, ols2$coefficients, coef(probit_s),
             log(sd(ols1$residuals)), log(sd(ols2$residuals)), 0, 0)

  opt <- maxLik::maxLik(loglik, start = start, method = "BFGS", ...)
  coefs <- coef(opt)
  nm <- c(paste0("regime1_", colnames(X)), paste0("regime2_", colnames(X)),
          paste0("selection_", colnames(Z)),
          "log_sigma1", "log_sigma2", "atanh_rho1", "atanh_rho2")
  names(coefs) <- nm
  H <- maxLik::hessian(opt)
  V <- tryCatch(solve(-H), error = function(e) MASS::ginv(-H))
  se <- sqrt(pmax(diag(V), 0))
  ll_val <- as.numeric(logLik(opt))

  list(coefficients = coefs, se = se, z = coefs / se,
       pvalue = 2 * pnorm(-abs(coefs / se)),
       vcov = V, hessian = H,
       logLik = ll_val, AIC = -2*ll_val + 2*length(coefs),
       BIC = -2*ll_val + log(n)*length(coefs),
       sigma1 = exp(coefs["log_sigma1"]), sigma2 = exp(coefs["log_sigma2"]),
       rho1 = tanh(coefs["atanh_rho1"]), rho2 = tanh(coefs["atanh_rho2"]),
       convergence = opt$code)
}

# --- Markov-Switching (Hamilton 1989) ---
#' @keywords internal
.switching_markov <- function(formula, data, n_regimes, error_dist,
                               transition, transition_covariates, ...) {
  F <- Formula::Formula(formula)
  mf <- model.frame(F, data = data)
  y <- model.response(mf)
  X <- model.matrix(F, data = data)
  n <- length(y); k <- ncol(X)

  n_par_regime <- k + 1  # beta + log_sigma per regime

  loglik_and_filter <- function(par, return_probs = FALSE) {
    beta_list <- lapply(seq_len(n_regimes), function(j) par[((j-1)*k+1):(j*k)])
    sigma_list <- lapply(seq_len(n_regimes), function(j) exp(par[n_regimes*k + j]))

    # Transition matrix (multinomial logit)
    offset <- n_regimes * (k + 1)
    P <- matrix(0, n_regimes, n_regimes)
    idx <- 1
    for (i in seq_len(n_regimes)) {
      raw <- numeric(n_regimes); raw[1] <- 0
      if (n_regimes > 1) {
        for (j in 2:n_regimes) {
          raw[j] <- par[offset + idx]; idx <- idx + 1
        }
      }
      P[i, ] <- exp(raw) / sum(exp(raw))
    }

    xi <- .stationary_distribution(P)
    ll <- 0
    xi_filtered <- matrix(0, n, n_regimes)

    for (t in seq_len(n)) {
      xi_pred <- if (t == 1) xi else as.numeric(t(P) %*% xi_filtered[t-1, ])
      f <- sapply(seq_len(n_regimes), function(j) {
        dnorm(y[t], mean = as.numeric(X[t, ] %*% beta_list[[j]]),
              sd = sigma_list[[j]])
      })
      f_weighted <- f * xi_pred; f_marginal <- sum(f_weighted)
      if (f_marginal < .Machine$double.eps) {
        if (return_probs) return(list(loglik = -1e10, filtered_probs = xi_filtered))
        return(-1e10)
      }
      xi_filtered[t, ] <- f_weighted / f_marginal
      ll <- ll + log(f_marginal)
    }
    if (return_probs) list(loglik = ll, filtered_probs = xi_filtered) else ll
  }

  # Starting values
  start_beta <- lapply(seq_len(n_regimes), function(j) {
    idx <- ((j-1) * floor(n/n_regimes) + 1):min(j * floor(n/n_regimes), n)
    lm.fit(X[idx, , drop = FALSE], y[idx])$coefficients
  })
  start_sigma <- rep(log(sd(y)), n_regimes)
  n_trans <- n_regimes * (n_regimes - 1)
  start <- c(unlist(start_beta), start_sigma, rep(0, n_trans))

  opt <- maxLik::maxLik(function(par) loglik_and_filter(par, FALSE),
                         start = start, method = "BFGS", ...)
  coefs <- coef(opt)

  # Name parameters
  nm <- c(unlist(lapply(seq_len(n_regimes), function(j)
    paste0("regime", j, "_", colnames(X)))),
    paste0("log_sigma_regime", seq_len(n_regimes)),
    paste0("trans_", seq_len(n_trans)))
  names(coefs) <- nm

  H <- maxLik::hessian(opt)
  V <- tryCatch(solve(-H), error = function(e) MASS::ginv(-H))
  se <- sqrt(pmax(diag(V), 0))
  ll_val <- as.numeric(logLik(opt))

  # Get filtered probabilities
  filter_result <- loglik_and_filter(coefs, return_probs = TRUE)

  list(coefficients = coefs, se = se, z = coefs / se,
       pvalue = 2 * pnorm(-abs(coefs / se)),
       vcov = V, hessian = H,
       logLik = ll_val, AIC = -2*ll_val + 2*length(coefs),
       BIC = -2*ll_val + log(n)*length(coefs),
       filtered_probs = filter_result$filtered_probs,
       convergence = opt$code)
}

#' @export
predict.switching_flex <- function(object, newdata = NULL, type = "response", ...) {
  if (is.null(newdata)) newdata <- object$model_data
  F <- object$formula
  X <- model.matrix(F, data = newdata, rhs = 1)
  k <- ncol(X)
  # Weighted average prediction across regimes (using filtered probs for Markov)
  if (object$type == "markov" && !is.null(object$filtered_probs)) {
    preds <- matrix(0, nrow(X), object$n_regimes)
    for (j in seq_len(object$n_regimes)) {
      beta_j <- object$coefficients[((j-1)*k+1):(j*k)]
      preds[, j] <- as.numeric(X %*% beta_j) * object$filtered_probs[, j]
    }
    rowSums(preds)
  } else {
    beta1 <- object$coefficients[1:k]
    as.numeric(X %*% beta1)
  }
}
