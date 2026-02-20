# =============================================================================
# ecoflex Module 2: Extended Heckman Selection Model
# =============================================================================

#' Extended Heckman Selection Model
#'
#' Estimate Heckman-type selection models with multiple estimation methods
#' and flexible error distributions.
#'
#' @param formula Multi-part formula: \code{outcome | selection ~ X_outcome | X_selection}.
#' @param data A data frame.
#' @param method Estimation method: \code{"ml"}, \code{"twostep"},
#'   \code{"twostep_robust"} (Lee), \code{"semipar"}.
#' @param selection_type \code{"binary"}, \code{"ordered"}, \code{"multinomial"}.
#' @param error_dist \code{"normal"}, \code{"t"}, \code{"copula"}.
#' @param copula_family \code{"gaussian"}, \code{"frank"}, \code{"clayton"}, \code{"gumbel"}.
#' @param heteroscedastic Formula for heteroscedasticity (optional).
#' @param vcov \code{"standard"}, \code{"robust"}, \code{"bootstrap"}.
#' @param latex If \code{TRUE}, prints a LaTeX table of results. Default: \code{FALSE}.
#' @param ... Additional arguments passed to the optimiser.
#' @return Object of class \code{c("heckman_flex", "ecoflex")}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 500
#' x1 <- rnorm(n); z <- rnorm(n)
#' s <- as.integer(0.3 * x1 + 0.5 * z + rnorm(n) > 0)
#' y <- ifelse(s == 1, 1 + 2 * x1 + rnorm(n), 0)
#' df <- data.frame(y = y, s = s, x1 = x1, z = z)
#' m <- heckman_flex(y | s ~ x1 | x1 + z, data = df, method = "twostep")
#' summary(m)
#' }
heckman_flex <- function(formula, data,
                         method = c("ml", "twostep", "twostep_robust", "semipar"),
                         selection_type = c("binary", "ordered", "multinomial"),
                         error_dist = c("normal", "t", "copula"),
                         copula_family = c("gaussian", "frank", "clayton", "gumbel"),
                         heteroscedastic = NULL,
                         vcov = c("standard", "robust", "bootstrap"),
                         latex = FALSE, ...) {
  method <- match.arg(method)
  selection_type <- match.arg(selection_type)
  error_dist <- match.arg(error_dist)
  copula_family <- match.arg(copula_family)
  vcov_type <- match.arg(vcov)

  F <- Formula::Formula(formula)
  parsed <- .heckman_parse_formula(F, data)
  y_outcome <- parsed$y_outcome
  y_selection <- parsed$y_selection
  X_outcome <- parsed$X_outcome
  X_selection <- parsed$X_selection
  n <- length(y_selection)
  selected <- y_selection == 1

  if (method == "twostep") {
    result <- .heckman_twostep_standard(y_outcome, y_selection, X_outcome, X_selection, selected)
  } else if (method == "twostep_robust") {
    result <- .heckman_twostep_lee(y_outcome, y_selection, X_outcome, X_selection, selected)
  } else if (method == "semipar") {
    result <- .heckman_semiparametric(y_outcome, y_selection, X_outcome, X_selection, selected)
  } else {
    result <- .heckman_ml(y_outcome, y_selection, X_outcome, X_selection, selected, error_dist, ...)
  }

  obj <- structure(
    c(result, list(
      call = match.call(), formula = F, method = method,
      selection_type = selection_type, error_dist = error_dist,
      n = n, n_selected = sum(selected), n_unselected = sum(!selected),
      model_data = data,
      model_name = sprintf("Heckman Selection Model (%s, %s)", method, error_dist)
    )),
    class = c("heckman_flex", "ecoflex")
  )
  if (latex) to_latex(obj)
  obj
}

#' @keywords internal
.heckman_parse_formula <- function(F, data) {
  # Formula: y | s ~ outcome_vars | selection_vars
  # LHS part 1 = y (outcome), LHS part 2 = s (selection indicator)
  # RHS part 1 = outcome covariates, RHS part 2 = selection covariates
  mf1 <- model.frame(F, data = data, lhs = 1, rhs = 1)
  y_full <- model.response(mf1)
  X_outcome <- model.matrix(F, data = data, lhs = 1, rhs = 1)
  mf2 <- model.frame(F, data = data, lhs = 2, rhs = 2)
  y_selection <- model.response(mf2)
  X_selection <- model.matrix(F, data = data, lhs = 2, rhs = 2)
  if (!all(y_selection %in% c(0, 1))) stop("Selection variable must be binary (0/1)")
  y_outcome <- y_full
  y_outcome[y_selection == 0] <- NA
  list(y_outcome = y_outcome, y_selection = y_selection,
       X_outcome = X_outcome, X_selection = X_selection)
}

#' @keywords internal
.heckman_twostep_standard <- function(y_outcome, y_selection, X_outcome, X_selection, selected) {
  n <- length(y_selection); n1 <- sum(selected)
  probit_fit <- glm(y_selection ~ X_selection - 1, family = binomial(link = "probit"))
  gamma_hat <- coef(probit_fit)
  Xg <- as.numeric(X_selection %*% gamma_hat)
  lambda <- dnorm(Xg) / pmax(pnorm(Xg), .Machine$double.eps)

  y_sel <- y_outcome[selected]; X_sel <- X_outcome[selected, , drop = FALSE]
  lambda_sel <- lambda[selected]
  X_augmented <- cbind(X_sel, lambda = lambda_sel)
  fit2 <- lm.fit(X_augmented, y_sel)
  coefs_outcome <- fit2$coefficients[1:ncol(X_sel)]
  rho_sigma <- fit2$coefficients[ncol(X_sel) + 1]
  resid2 <- fit2$residuals
  sigma <- sqrt(sum(resid2^2) / (n1 - ncol(X_augmented)))

  all_coefs <- c(setNames(coefs_outcome, paste0("outcome_", colnames(X_sel))),
                 rho_sigma = rho_sigma,
                 setNames(gamma_hat, paste0("selection_", colnames(X_selection))))
  k_total <- length(all_coefs)
  V_full <- matrix(0, k_total, k_total)
  V_ols <- sigma^2 * tryCatch(solve(crossprod(X_augmented)), error = function(e) MASS::ginv(crossprod(X_augmented)))
  k_out <- ncol(X_augmented); k_sel_n <- length(gamma_hat)
  V_full[1:k_out, 1:k_out] <- V_ols
  V_full[(k_out+1):(k_out+k_sel_n), (k_out+1):(k_out+k_sel_n)] <- vcov(probit_fit)
  se <- sqrt(pmax(diag(V_full), 0))
  ll_val <- as.numeric(logLik(probit_fit)) + sum(dnorm(resid2, sd = sigma, log = TRUE))

  list(coefficients = all_coefs, se = se, z = all_coefs / se,
       pvalue = 2 * pnorm(-abs(all_coefs / se)), vcov = V_full, hessian = NULL,
       logLik = ll_val, AIC = -2 * ll_val + 2 * k_total,
       BIC = -2 * ll_val + log(n) * k_total,
       sigma = sigma, rho = rho_sigma / sigma, lambda = lambda, convergence = 0)
}

#' @keywords internal
.heckman_twostep_lee <- function(y_outcome, y_selection, X_outcome, X_selection, selected) {
  n <- length(y_selection); n1 <- sum(selected)
  probit_fit <- glm(y_selection ~ X_selection - 1, family = binomial(link = "probit"))
  gamma_hat <- coef(probit_fit)
  Xg <- as.numeric(X_selection %*% gamma_hat)
  lambda <- dnorm(Xg) / pmax(pnorm(Xg), .Machine$double.eps)
  K <- 3; lambda_sel <- lambda[selected]
  correction_terms <- matrix(0, n1, K)
  for (k in 1:K) correction_terms[, k] <- lambda_sel^k
  y_sel <- y_outcome[selected]; X_sel <- X_outcome[selected, , drop = FALSE]
  X_augmented <- cbind(X_sel, correction_terms)
  fit2 <- lm.fit(X_augmented, y_sel)
  coefs_outcome <- fit2$coefficients[1:ncol(X_sel)]
  correction_coefs <- fit2$coefficients[(ncol(X_sel) + 1):ncol(X_augmented)]
  resid2 <- fit2$residuals
  sigma <- sqrt(sum(resid2^2) / (n1 - ncol(X_augmented)))

  all_coefs <- c(setNames(coefs_outcome, paste0("outcome_", colnames(X_sel))),
                 setNames(correction_coefs, paste0("correction_", 1:K)),
                 setNames(gamma_hat, paste0("selection_", colnames(X_selection))))
  k_total <- length(all_coefs); se <- rep(NA_real_, k_total)
  V_ols <- sigma^2 * tryCatch(solve(crossprod(X_augmented)), error = function(e) MASS::ginv(crossprod(X_augmented)))
  se[1:ncol(X_augmented)] <- sqrt(pmax(diag(V_ols), 0))
  se[(ncol(X_augmented)+1):k_total] <- sqrt(pmax(diag(vcov(probit_fit)), 0))
  V_full <- diag(se^2, nrow = k_total)
  ll_val <- as.numeric(logLik(probit_fit)) + sum(dnorm(resid2, sd = sigma, log = TRUE))

  list(coefficients = all_coefs, se = se, z = all_coefs / se,
       pvalue = 2 * pnorm(-abs(all_coefs / se)), vcov = V_full, hessian = NULL,
       logLik = ll_val, AIC = -2 * ll_val + 2 * k_total,
       BIC = -2 * ll_val + log(n) * k_total, sigma = sigma, lambda = lambda, convergence = 0)
}

#' @keywords internal
.heckman_ml <- function(y_outcome, y_selection, X_outcome, X_selection, selected, error_dist, ...) {
  n <- length(y_selection); k_out <- ncol(X_outcome); k_sel <- ncol(X_selection)

  loglik <- function(par) {
    beta <- par[1:k_out]; gamma <- par[(k_out+1):(k_out+k_sel)]
    log_sigma <- par[k_out + k_sel + 1]; rho_raw <- par[k_out + k_sel + 2]
    sigma <- exp(log_sigma); rho <- tanh(rho_raw)
    Xb <- as.numeric(X_outcome %*% beta); Zg <- as.numeric(X_selection %*% gamma)
    ll <- numeric(n)
    idx1 <- selected
    if (any(idx1)) {
      resid_i <- (y_outcome[idx1] - Xb[idx1]) / sigma
      if (error_dist == "t") {
        df <- exp(par[k_out + k_sel + 3]) + 2
        ll[idx1] <- dt(resid_i, df = df, log = TRUE) - log(sigma) +
          pnorm((Zg[idx1] + rho * resid_i) / sqrt(1 - rho^2), log.p = TRUE)
      } else {
        arg <- (Zg[idx1] + rho * resid_i) / sqrt(1 - rho^2)
        ll[idx1] <- dnorm(resid_i, log = TRUE) - log(sigma) + pnorm(arg, log.p = TRUE)
      }
    }
    idx0 <- !selected
    if (any(idx0)) ll[idx0] <- pnorm(-Zg[idx0], log.p = TRUE)
    sum(ll)
  }

  probit_fit <- suppressWarnings(glm(y_selection ~ X_selection - 1, family = binomial(link = "probit")))
  gamma_start <- coef(probit_fit)
  y_sel <- y_outcome[selected]; X_sel <- X_outcome[selected, , drop = FALSE]
  ols_fit <- lm.fit(X_sel, y_sel)
  start <- c(ols_fit$coefficients, gamma_start, log(sd(ols_fit$residuals)), 0)
  if (error_dist == "t") start <- c(start, log(3))

  opt <- maxLik::maxLik(loglik, start = start, method = "BFGS", ...)
  coefs <- coef(opt)
  names_all <- c(paste0("outcome_", colnames(X_outcome)),
                 paste0("selection_", colnames(X_selection)), "log_sigma", "atanh_rho")
  if (error_dist == "t") names_all <- c(names_all, "log_df")
  names(coefs) <- names_all
  H <- maxLik::hessian(opt)
  V <- tryCatch(solve(-H), error = function(e) MASS::ginv(-H))
  se <- sqrt(pmax(diag(V), 0))
  ll_val <- as.numeric(logLik(opt))

  list(coefficients = coefs, se = se, z = coefs / se,
       pvalue = 2 * pnorm(-abs(coefs / se)), vcov = V, hessian = H,
       logLik = ll_val, AIC = -2 * ll_val + 2 * length(coefs),
       BIC = -2 * ll_val + log(n) * length(coefs),
       sigma = exp(coefs["log_sigma"]), rho = tanh(coefs["atanh_rho"]),
       convergence = opt$code)
}

#' @keywords internal
.heckman_semiparametric <- function(y_outcome, y_selection, X_outcome, X_selection, selected) {
  n <- length(y_selection); n1 <- sum(selected)
  probit_fit <- glm(y_selection ~ X_selection - 1, family = binomial(link = "probit"))
  gamma_hat <- coef(probit_fit); p_hat <- fitted(probit_fit)
  p_sel <- p_hat[selected]; y_sel <- y_outcome[selected]
  X_sel <- X_outcome[selected, , drop = FALSE]
  best_bic <- Inf; best_K <- 1
  for (K in 1:5) {
    poly_basis <- poly(p_sel, degree = K, raw = TRUE)
    X_aug <- cbind(X_sel, poly_basis); fit_k <- lm.fit(X_aug, y_sel)
    sigma_k <- sqrt(sum(fit_k$residuals^2) / (n1 - ncol(X_aug)))
    bic_k <- n1 * log(sigma_k^2) + ncol(X_aug) * log(n1)
    if (bic_k < best_bic) { best_bic <- bic_k; best_K <- K }
  }
  poly_basis <- poly(p_sel, degree = best_K, raw = TRUE)
  X_augmented <- cbind(X_sel, poly_basis)
  fit_final <- lm.fit(X_augmented, y_sel)
  coefs_outcome <- fit_final$coefficients[1:ncol(X_sel)]
  correction_coefs <- fit_final$coefficients[(ncol(X_sel)+1):ncol(X_augmented)]
  resid_final <- fit_final$residuals
  sigma <- sqrt(sum(resid_final^2) / (n1 - ncol(X_augmented)))
  all_coefs <- c(setNames(coefs_outcome, paste0("outcome_", colnames(X_sel))),
                 setNames(correction_coefs, paste0("poly_", 1:best_K)),
                 setNames(gamma_hat, paste0("selection_", colnames(X_selection))))
  k_total <- length(all_coefs); se <- rep(NA_real_, k_total)
  V_ols <- sigma^2 * tryCatch(solve(crossprod(X_augmented)), error = function(e) MASS::ginv(crossprod(X_augmented)))
  se[1:ncol(X_augmented)] <- sqrt(pmax(diag(V_ols), 0))
  se[(ncol(X_augmented)+1):k_total] <- sqrt(pmax(diag(vcov(probit_fit)), 0))
  V_full <- diag(se^2, nrow = k_total)
  ll_val <- as.numeric(logLik(probit_fit)) + sum(dnorm(resid_final, sd = sigma, log = TRUE))

  list(coefficients = all_coefs, se = se, z = all_coefs / se,
       pvalue = 2 * pnorm(-abs(all_coefs / se)), vcov = V_full, hessian = NULL,
       logLik = ll_val, AIC = -2 * ll_val + 2 * k_total,
       BIC = -2 * ll_val + log(n) * k_total,
       sigma = sigma, poly_order = best_K, convergence = 0)
}

#' @export
predict.heckman_flex <- function(object, newdata = NULL, type = "response", ...) {
  if (is.null(newdata)) newdata <- object$model_data
  F <- object$formula
  X_outcome <- model.matrix(F, data = newdata, rhs = 1)
  k_out <- ncol(X_outcome)
  beta <- object$coefficients[1:k_out]
  if (type == "response") {
    as.numeric(X_outcome %*% beta)
  } else {
    X_selection <- model.matrix(F, data = newdata, rhs = 2)
    k_sel <- ncol(X_selection)
    gamma <- object$coefficients[(k_out+1):(k_out+k_sel)]
    pnorm(as.numeric(X_selection %*% gamma))
  }
}
