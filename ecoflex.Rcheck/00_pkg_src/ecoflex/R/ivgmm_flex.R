# =============================================================================
# ecoflex Module 4: Unified IV/GMM Estimation
# =============================================================================

#' Unified IV/GMM Estimation
#'
#' Estimate models with instrumental variables using 2SLS, LIML, Fuller,
#' or GMM with various weighting options.
#'
#' @param formula Formula: \code{outcome ~ endogenous | instruments}.
#' @param data A data frame.
#' @param method \code{"2sls"}, \code{"liml"}, \code{"fuller"},
#'   \code{"gmm_twostep"}, \code{"gmm_iterative"}, \code{"cue"}.
#' @param vcov \code{"classical"}, \code{"robust"}, \code{"HAC"},
#'   \code{"cluster"}, \code{"bootstrap"}.
#' @param weights_matrix For GMM: \code{"identity"}, \code{"unadjusted"},
#'   \code{"robust"}, \code{"HAC"}.
#' @param HAC_kernel \code{"bartlett"}, \code{"parzen"}, \code{"qs"}.
#' @param HAC_bw HAC bandwidth: numeric or \code{"auto"}.
#' @param fuller_alpha Fuller parameter (if \code{method = "fuller"}).
#' @param cluster Cluster variable name.
#' @param maxiter GMM max iterations.
#' @param tol GMM convergence tolerance.
#' @param latex If \code{TRUE}, prints a LaTeX table. Default: \code{FALSE}.
#' @param ... Additional arguments.
#' @return Object of class \code{c("ivgmm_flex", "ecoflex")} with
#'   diagnostic tests (weak IV, Sargan-Hansen).
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 300; z <- rnorm(n)
#' x <- 0.5 * z + rnorm(n)
#' y <- 1 + 2 * x + rnorm(n)
#' df <- data.frame(y = y, x = x, z = z)
#' m <- ivgmm_flex(y ~ x | z, data = df, method = "2sls")
#' summary(m)
#' m$diagnostics$weak_iv_test
#' }
ivgmm_flex <- function(formula, data,
                       method = c("2sls", "liml", "fuller", "gmm_twostep",
                                  "gmm_iterative", "cue"),
                       vcov = c("classical", "robust", "HAC", "cluster", "bootstrap"),
                       weights_matrix = "robust",
                       HAC_kernel = c("bartlett", "parzen", "qs"),
                       HAC_bw = "auto", fuller_alpha = 1,
                       cluster = NULL, maxiter = 100, tol = 1e-8, latex = FALSE, ...) {
  method <- match.arg(method)
  vcov_type <- match.arg(vcov)
  HAC_kernel <- match.arg(HAC_kernel)

  F <- Formula::Formula(formula)
  mf <- model.frame(F, data = data)
  y <- model.response(mf)
  X <- model.matrix(F, data = data, rhs = 1)
  Z <- model.matrix(F, data = data, rhs = 2)
  n <- length(y); k <- ncol(X); l <- ncol(Z)

  if (l < k) stop("Under-identified model: number of instruments < number of regressors")

  result <- switch(method,
    "2sls" = .estimate_2sls(y, X, Z),
    "liml" = .estimate_liml(y, X, Z),
    "fuller" = .estimate_fuller(y, X, Z, alpha = fuller_alpha),
    "gmm_twostep" = .estimate_gmm(y, X, Z, weights_matrix, iterative = FALSE),
    "gmm_iterative" = .estimate_gmm(y, X, Z, weights_matrix, iterative = TRUE,
                                      maxiter = maxiter, tol = tol),
    "cue" = .estimate_cue(y, X, Z, weights_matrix)
  )

  V <- switch(vcov_type,
    classical = result$vcov_classical,
    robust = .iv_robust_vcov(result, y, X, Z),
    HAC = .iv_hac_vcov(result, y, X, Z, kernel = HAC_kernel, bw = HAC_bw),
    cluster = .iv_cluster_vcov(result, y, X, Z, cluster = data[[cluster]]),
    bootstrap = .iv_bootstrap_vcov(result, formula, data, method, R = 500)
  )

  diagnostics <- list(
    overid_test = if (l > k) .sargan_hansen_test(result, y, X, Z) else NULL,
    weak_iv_test = .weak_iv_test(y, X, Z),
    endogeneity_test = .durbin_wu_hausman(y, X, Z),
    first_stage_F = .first_stage_diagnostics(y, X, Z)
  )

  coefs <- result$coefficients
  names(coefs) <- colnames(X)
  se <- sqrt(pmax(diag(V), 0))

  obj <- structure(list(
    coefficients = coefs, se = se, z = coefs / se,
    pvalue = 2 * pnorm(-abs(coefs / se)),
    vcov = V, residuals = result$residuals, fitted.values = result$fitted,
    diagnostics = diagnostics, call = match.call(), formula = F,
    method = method, n = n, k = k, l = l,
    model_data = data,
    model_name = sprintf("IV/GMM Estimation (%s)", method),
    hessian = NULL,
    logLik = -n/2 * log(2 * pi * sum(result$residuals^2)/n) - n/2,
    AIC = n * log(sum(result$residuals^2)/n) + 2 * k,
    BIC = n * log(sum(result$residuals^2)/n) + k * log(n)
  ), class = c("ivgmm_flex", "ecoflex"))
  if (latex) to_latex(obj)
  obj
}

# --- Estimators ---

#' @keywords internal
.estimate_2sls <- function(y, X, Z) {
  n <- length(y); k <- ncol(X)
  Pz <- Z %*% solve(crossprod(Z)) %*% t(Z)
  beta <- solve(t(X) %*% Pz %*% X) %*% t(X) %*% Pz %*% y
  beta <- as.numeric(beta)
  resid <- as.numeric(y - X %*% beta)
  sigma2 <- sum(resid^2) / (n - k)
  V <- sigma2 * solve(t(X) %*% Pz %*% X)
  list(coefficients = beta, residuals = resid, fitted = as.numeric(X %*% beta),
       vcov_classical = V, Pz = Pz)
}

#' @keywords internal
.estimate_liml <- function(y, X, Z) {
  n <- length(y); k <- ncol(X)
  Pz <- Z %*% solve(crossprod(Z)) %*% t(Z)
  Mz <- diag(n) - Pz
  # LIML: find kappa as smallest eigenvalue of (Y'Mz Y)^{-1} (Y'Pz Y)
  Y <- cbind(y, X[, -1, drop = FALSE])  # treat first col as intercept
  A <- crossprod(Y, Pz %*% Y)
  B <- crossprod(Y, Mz %*% Y)
  eig <- eigen(solve(B) %*% A)
  kappa <- min(Re(eig$values))
  # LIML estimator
  W <- Pz - kappa * Mz
  beta <- solve(t(X) %*% (diag(n) + W) %*% X) %*% t(X) %*% (diag(n) + W) %*% y
  # Simplified: use 2SLS-like formula with kappa adjustment
  beta <- tryCatch({
    solve(t(X) %*% Pz %*% X - (kappa - 1) * t(X) %*% Mz %*% X) %*%
      (t(X) %*% Pz %*% y - (kappa - 1) * t(X) %*% Mz %*% y)
  }, error = function(e) .estimate_2sls(y, X, Z)$coefficients)
  beta <- as.numeric(beta)
  resid <- as.numeric(y - X %*% beta)
  sigma2 <- sum(resid^2) / (n - k)
  V <- sigma2 * tryCatch(solve(t(X) %*% Pz %*% X),
                          error = function(e) MASS::ginv(t(X) %*% Pz %*% X))
  list(coefficients = beta, residuals = resid, fitted = as.numeric(X %*% beta),
       vcov_classical = V, kappa = kappa, Pz = Pz)
}

#' @keywords internal
.estimate_fuller <- function(y, X, Z, alpha = 1) {
  n <- length(y); k <- ncol(X)
  liml_result <- .estimate_liml(y, X, Z)
  kappa <- liml_result$kappa
  kappa_fuller <- kappa - alpha / (n - ncol(Z))
  Pz <- liml_result$Pz; Mz <- diag(n) - Pz
  beta <- tryCatch({
    solve(t(X) %*% Pz %*% X - (kappa_fuller - 1) * t(X) %*% Mz %*% X) %*%
      (t(X) %*% Pz %*% y - (kappa_fuller - 1) * t(X) %*% Mz %*% y)
  }, error = function(e) liml_result$coefficients)
  beta <- as.numeric(beta)
  resid <- as.numeric(y - X %*% beta)
  sigma2 <- sum(resid^2) / (n - k)
  V <- sigma2 * tryCatch(solve(t(X) %*% Pz %*% X),
                          error = function(e) MASS::ginv(t(X) %*% Pz %*% X))
  list(coefficients = beta, residuals = resid, fitted = as.numeric(X %*% beta),
       vcov_classical = V, kappa_fuller = kappa_fuller, Pz = Pz)
}

#' @keywords internal
.estimate_gmm <- function(y, X, Z, weights_type, iterative = FALSE,
                           maxiter = 100, tol = 1e-8) {
  n <- length(y)
  Pz <- Z %*% solve(crossprod(Z)) %*% t(Z)
  beta <- solve(t(X) %*% Pz %*% X) %*% t(X) %*% Pz %*% y
  e <- as.numeric(y - X %*% beta)
  W <- .compute_weight_matrix(Z, e, type = weights_type)

  if (!iterative) {
    beta <- solve(t(X) %*% Z %*% W %*% t(Z) %*% X) %*%
            t(X) %*% Z %*% W %*% t(Z) %*% y
  } else {
    for (i in seq_len(maxiter)) {
      beta_old <- beta
      e <- as.numeric(y - X %*% beta)
      W <- .compute_weight_matrix(Z, e, type = weights_type)
      beta <- tryCatch(
        solve(t(X) %*% Z %*% W %*% t(Z) %*% X) %*%
          t(X) %*% Z %*% W %*% t(Z) %*% y,
        error = function(e) beta_old
      )
      if (max(abs(beta - beta_old)) < tol) break
    }
  }
  beta <- as.numeric(beta)
  e <- as.numeric(y - X %*% beta)
  V <- tryCatch(solve(t(X) %*% Z %*% W %*% t(Z) %*% X),
                error = function(er) MASS::ginv(t(X) %*% Z %*% W %*% t(Z) %*% X))
  list(coefficients = beta, residuals = e, fitted = as.numeric(X %*% beta),
       vcov_classical = V, W = W, Pz = Z %*% solve(crossprod(Z)) %*% t(Z))
}

#' @keywords internal
.estimate_cue <- function(y, X, Z, weights_type) {
  n <- length(y); k <- ncol(X)
  # CUE: minimize Q(beta) = g(beta)' W(beta)^{-1} g(beta) jointly
  g_fn <- function(beta) {
    e <- as.numeric(y - X %*% beta)
    as.numeric(t(Z) %*% e / n)
  }
  obj_fn <- function(beta) {
    e <- as.numeric(y - X %*% beta)
    g <- as.numeric(t(Z) %*% e / n)
    S <- crossprod(Z * e) / n
    S_inv <- tryCatch(solve(S), error = function(er) MASS::ginv(S))
    as.numeric(t(g) %*% S_inv %*% g) * n
  }
  start <- .estimate_2sls(y, X, Z)$coefficients
  opt <- optim(start, obj_fn, method = "BFGS")
  beta <- opt$par
  e <- as.numeric(y - X %*% beta)
  W <- .compute_weight_matrix(Z, e, type = weights_type)
  V <- tryCatch(solve(t(X) %*% Z %*% W %*% t(Z) %*% X),
                error = function(er) MASS::ginv(t(X) %*% Z %*% W %*% t(Z) %*% X))
  list(coefficients = beta, residuals = e, fitted = as.numeric(X %*% beta),
       vcov_classical = V, W = W, Pz = Z %*% solve(crossprod(Z)) %*% t(Z))
}

# --- Vcov helpers ---

#' @keywords internal
.iv_robust_vcov <- function(result, y, X, Z) {
  n <- length(y); k <- ncol(X); e <- result$residuals
  Pz <- result$Pz; XPzX_inv <- tryCatch(solve(t(X) %*% Pz %*% X),
                                          error = function(er) MASS::ginv(t(X) %*% Pz %*% X))
  meat <- t(X) %*% Pz %*% diag(e^2) %*% Pz %*% X
  n / (n - k) * XPzX_inv %*% meat %*% XPzX_inv
}

#' @keywords internal
.iv_hac_vcov <- function(result, y, X, Z, kernel = "bartlett", bw = "auto") {
  n <- length(y); e <- result$residuals
  if (is.character(bw) && bw == "auto") bw <- floor(n^(1/3))
  XZ <- t(X) %*% Z; ZX <- t(Z) %*% X
  Ze <- Z * e
  S <- crossprod(Ze) / n
  for (j in seq_len(bw)) {
    w <- switch(kernel, bartlett = 1 - j/(bw+1), parzen = {
      r <- j/(bw+1); if (r <= 0.5) 1-6*r^2+6*r^3 else 2*(1-r)^3
    }, qs = 25/(12*pi^2*(j/(bw+1))^2) * (sin(6*pi*(j/(bw+1))/5)/(6*pi*(j/(bw+1))/5) - cos(6*pi*(j/(bw+1))/5)))
    Gamma_j <- crossprod(Ze[-(1:j), , drop = FALSE], Ze[1:(n-j), , drop = FALSE]) / n
    S <- S + w * (Gamma_j + t(Gamma_j))
  }
  ZZ_inv <- solve(crossprod(Z))
  bread <- solve(XZ %*% ZZ_inv %*% ZX)
  bread %*% XZ %*% ZZ_inv %*% S %*% ZZ_inv %*% ZX %*% bread * n
}

#' @keywords internal
.iv_cluster_vcov <- function(result, y, X, Z, cluster) {
  if (is.null(cluster)) stop("Cluster variable required")
  n <- length(y); k <- ncol(X); e <- result$residuals
  Pz <- result$Pz
  XPzX_inv <- tryCatch(solve(t(X) %*% Pz %*% X), error = function(er) MASS::ginv(t(X) %*% Pz %*% X))
  cl_ids <- unique(cluster); G <- length(cl_ids)
  meat <- Reduce("+", lapply(cl_ids, function(g) {
    idx <- which(cluster == g)
    Xi_Pz_e <- t(X[idx, , drop = FALSE]) %*% Pz[idx, ] %*% (e * Pz[, idx]) # approximate
    # Simplified: use Z-projected residuals
    Zi_ei <- colSums(Z[idx, , drop = FALSE] * e[idx])
    XZ_inv_Ze <- as.numeric(XPzX_inv %*% t(X) %*% Z %*% solve(crossprod(Z)) %*% Zi_ei)
    tcrossprod(XZ_inv_Ze)
  }))
  G / (G - 1) * (n - 1) / (n - k) * meat
}

#' @keywords internal
.iv_bootstrap_vcov <- function(result, formula, data, method, R = 500) {
  n <- nrow(data); k <- length(result$coefficients)
  boot_coefs <- matrix(NA_real_, R, k)
  for (i in seq_len(R)) {
    idx <- sample(n, replace = TRUE)
    d <- data[idx, ]
    fit <- tryCatch(ivgmm_flex(formula, data = d, method = method), error = function(e) NULL)
    if (!is.null(fit)) boot_coefs[i, ] <- coef(fit)
  }
  boot_coefs <- boot_coefs[complete.cases(boot_coefs), , drop = FALSE]
  if (nrow(boot_coefs) < 10) warning("Few successful bootstrap replications")
  cov(boot_coefs)
}

# --- Diagnostic tests ---

#' @keywords internal
.sargan_hansen_test <- function(result, y, X, Z) {
  n <- length(y); l <- ncol(Z); k <- ncol(X)
  e <- result$residuals
  S <- crossprod(Z * e) / n
  g <- as.numeric(t(Z) %*% e / n)
  S_inv <- tryCatch(solve(S), error = function(er) MASS::ginv(S))
  J <- n * as.numeric(t(g) %*% S_inv %*% g)
  df <- l - k
  list(statistic = J, df = df, p.value = 1 - pchisq(J, df))
}

#' @keywords internal
.weak_iv_test <- function(y, X, Z) {
  k <- ncol(X); l <- ncol(Z)
  # Cragg-Donald: first-stage F-statistics
  fs <- .first_stage_diagnostics(y, X, Z)
  list(first_stage_F = fs, critical_values = c("10%" = 16.38, "15%" = 8.96, "20%" = 6.66, "25%" = 5.53))
}

#' @keywords internal
.durbin_wu_hausman <- function(y, X, Z) {
  n <- length(y); k <- ncol(X)
  # OLS
  beta_ols <- solve(crossprod(X)) %*% crossprod(X, y)
  # 2SLS
  Pz <- Z %*% solve(crossprod(Z)) %*% t(Z)
  beta_iv <- solve(t(X) %*% Pz %*% X) %*% t(X) %*% Pz %*% y
  diff <- beta_iv - beta_ols
  e_iv <- as.numeric(y - X %*% beta_iv)
  sigma2 <- sum(e_iv^2) / (n - k)
  V_diff <- sigma2 * (solve(t(X) %*% Pz %*% X) - solve(crossprod(X)))
  V_diff <- tryCatch(V_diff, error = function(e) diag(k) * sigma2)
  chi2 <- tryCatch(as.numeric(t(diff) %*% solve(V_diff) %*% diff),
                   error = function(e) NA)
  list(statistic = chi2, df = k, p.value = if (!is.na(chi2)) 1 - pchisq(chi2, k) else NA)
}

#' @keywords internal
.first_stage_diagnostics <- function(y, X, Z) {
  k <- ncol(X)
  f_stats <- sapply(seq_len(k), function(j) {
    if (j == 1 && colnames(X)[1] == "(Intercept)") return(NA)
    x_j <- X[, j]
    X_exog <- X[, -j, drop = FALSE]
    resid_x <- residuals(lm.fit(Z, x_j))
    resid_z <- residuals(lm.fit(X_exog, x_j))
    n <- length(x_j); l <- ncol(Z)
    ssr_r <- sum(resid_z^2); ssr_u <- sum(resid_x^2)
    ((ssr_r - ssr_u) / (l - ncol(X_exog))) / (ssr_u / (n - l))
  })
  names(f_stats) <- colnames(X)
  f_stats
}

#' @export
predict.ivgmm_flex <- function(object, newdata = NULL, type = "response", ...) {
  if (is.null(newdata)) return(object$fitted.values)
  X_new <- model.matrix(object$formula, data = newdata, rhs = 1)
  as.numeric(X_new %*% object$coefficients)
}
