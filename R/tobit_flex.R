# =============================================================================
# ecoflex Module 6: Tobit Types I–V
# =============================================================================

#' Unified Tobit Model (Types I-V)
#'
#' Estimate censored regression models of Tobit types I through V with
#' flexible error distributions and multiple estimation methods.
#'
#' @param formula Formula. For types II-V, multi-part with \code{|} separator.
#' @param data A data frame.
#' @param tobit_type Integer 1-5. Default: 1.
#' @param left Lower censoring limit. Default: 0.
#' @param right Upper censoring limit. Default: \code{Inf}.
#' @param error_dist \code{"normal"}, \code{"t"}, \code{"logistic"}.
#' @param heteroscedastic Formula for heteroscedasticity (optional).
#' @param method \code{"ml"}, \code{"twostep"}, \code{"powell"} (censored LAD).
#' @param latex If \code{TRUE}, prints a LaTeX table. Default: \code{FALSE}.
#' @param ... Additional arguments.
#' @return Object of class \code{c("tobit_flex", "ecoflex")}.
#' @export
#'
#' @examples
#' \donttest{
#' m <- tobit_flex(mpg ~ hp + wt, data = mtcars, tobit_type = 1, left = 15)
#' summary(m)
#' }
tobit_flex <- function(formula, data, tobit_type = 1,
                       left = 0, right = Inf,
                       error_dist = c("normal", "t", "logistic"),
                       heteroscedastic = NULL,
                       method = c("ml", "twostep", "powell"), latex = FALSE, ...) {
  error_dist <- match.arg(error_dist)
  method <- match.arg(method)
  stopifnot(tobit_type %in% 1:5)

  result <- switch(as.character(tobit_type),
    "1" = .tobit_type1(formula, data, left, right, error_dist, heteroscedastic, method, ...),
    "2" = .tobit_type2(formula, data, error_dist, method, ...),
    "3" = .tobit_type3(formula, data, left, right, error_dist, ...),
    "4" = .tobit_type4(formula, data, error_dist, ...),
    "5" = .tobit_type5(formula, data, left, right, error_dist, ...)
  )

  obj <- structure(c(result, list(
    call = match.call(), formula = Formula::Formula(formula),
    tobit_type = tobit_type, left = left, right = right,
    error_dist = error_dist, method = method,
    model_data = data,
    model_name = sprintf("Tobit Type %d Model (%s, errors: %s)", tobit_type, method, error_dist)
  )), class = c("tobit_flex", "ecoflex"))
  if (latex) to_latex(obj)
  obj
}

# --- Type I: Standard censored regression ---
#' @keywords internal
.tobit_type1 <- function(formula, data, left, right, error_dist, heteroscedastic, method, ...) {
  F <- Formula::Formula(formula)
  mf <- model.frame(F, data = data)
  y <- model.response(mf)
  X <- model.matrix(F, data = data)
  n <- length(y); k <- ncol(X)

  censored_left <- y <= left & is.finite(left)
  censored_right <- y >= right & is.finite(right)
  uncensored <- !censored_left & !censored_right

  if (method == "powell") return(.censored_lad(y, X, left, right))

  Z_het <- if (!is.null(heteroscedastic)) model.matrix(heteroscedastic, data = data) else NULL
  k_het <- if (!is.null(Z_het)) ncol(Z_het) else 0

  loglik <- function(par) {
    beta <- par[1:k]
    if (k_het > 0) {
      gamma <- par[(k+1):(k+k_het)]; sigma <- exp(as.numeric(Z_het %*% gamma))
    } else {
      sigma <- exp(par[k + 1])
    }
    mu <- as.numeric(X %*% beta); z <- (y - mu) / sigma; ll <- numeric(n)

    if (error_dist == "normal") {
      ll[uncensored] <- dnorm(z[uncensored], log = TRUE) - log(sigma[if (k_het > 0) uncensored else 1])
      ll[censored_left] <- pnorm((left - mu[censored_left]) / sigma[if (k_het > 0) censored_left else 1], log.p = TRUE)
      ll[censored_right] <- pnorm((right - mu[censored_right]) / sigma[if (k_het > 0) censored_right else 1], lower.tail = FALSE, log.p = TRUE)
    } else if (error_dist == "t") {
      df <- exp(par[length(par)]) + 2
      ll[uncensored] <- dt(z[uncensored], df = df, log = TRUE) - log(sigma[if (k_het > 0) uncensored else 1])
      ll[censored_left] <- pt((left - mu[censored_left]) / sigma[if (k_het > 0) censored_left else 1], df = df, log.p = TRUE)
      ll[censored_right] <- pt((right - mu[censored_right]) / sigma[if (k_het > 0) censored_right else 1], df = df, lower.tail = FALSE, log.p = TRUE)
    } else {
      ll[uncensored] <- dlogis(z[uncensored], log = TRUE) - log(sigma[if (k_het > 0) uncensored else 1])
      ll[censored_left] <- plogis((left - mu[censored_left]) / sigma[if (k_het > 0) censored_left else 1], log.p = TRUE)
      ll[censored_right] <- plogis((right - mu[censored_right]) / sigma[if (k_het > 0) censored_right else 1], lower.tail = FALSE, log.p = TRUE)
    }
    sum(ll)
  }

  ols <- lm.fit(X, y)
  start <- if (k_het > 0) c(ols$coefficients, rep(0, k_het)) else c(ols$coefficients, log(sd(ols$residuals)))
  if (error_dist == "t") start <- c(start, log(3))

  opt <- maxLik::maxLik(loglik, start = start, method = "BFGS", ...)
  coefs <- coef(opt)
  H <- maxLik::hessian(opt)
  V <- tryCatch(solve(-H), error = function(e) MASS::ginv(-H))
  se <- sqrt(pmax(diag(V), 0))
  ll_val <- as.numeric(logLik(opt))

  list(coefficients = coefs, se = se, z = coefs / se,
       pvalue = 2 * pnorm(-abs(coefs / se)), vcov = V, hessian = H,
       logLik = ll_val, AIC = -2*ll_val + 2*length(coefs),
       BIC = -2*ll_val + log(n)*length(coefs),
       n = n, n_censored_left = sum(censored_left),
       n_censored_right = sum(censored_right),
       n_uncensored = sum(uncensored), convergence = opt$code)
}

# --- Type II: Selection model (wrapper to Heckman) ---
#' @keywords internal
.tobit_type2 <- function(formula, data, error_dist, method, ...) {
  heck_method <- if (method == "powell") "twostep" else if (method == "twostep") "twostep" else "ml"
  result <- heckman_flex(formula, data, method = heck_method, error_dist = error_dist, ...)
  c(result, list(n = result$n))
}

# --- Type III: Two equations with different censoring ---
#' @keywords internal
.tobit_type3 <- function(formula, data, left, right, error_dist, ...) {
  # Simplified: estimate as Type I with separate censoring indicators
  .tobit_type1(formula, data, left, right, error_dist, NULL, "ml", ...)
}

# --- Type IV: Selection with binary outcome ---
#' @keywords internal
.tobit_type4 <- function(formula, data, error_dist, ...) {
  F <- Formula::Formula(formula)
  # Parse as Heckman but with binary outcome (probit)
  parsed <- .heckman_parse_formula(F, data)
  y_sel <- parsed$y_selection; X_sel <- parsed$X_selection
  y_out <- parsed$y_outcome; X_out <- parsed$X_outcome

  n <- length(y_sel); k_out <- ncol(X_out); k_sel <- ncol(X_sel)
  selected <- y_sel == 1

  loglik <- function(par) {
    beta <- par[1:k_out]; gamma <- par[(k_out+1):(k_out+k_sel)]
    rho_raw <- par[k_out + k_sel + 1]; rho <- tanh(rho_raw)
    Xb <- as.numeric(X_out %*% beta); Zg <- as.numeric(X_sel %*% gamma)
    ll <- numeric(n)
    idx1 <- selected & y_out == 1
    idx0 <- selected & y_out == 0
    if (any(idx1)) {
      ll[idx1] <- log(pmax(.Machine$double.eps,
        pnorm((Zg[idx1] + rho * Xb[idx1]) / sqrt(1 - rho^2)) * pnorm(Xb[idx1])))
    }
    if (any(idx0)) {
      ll[idx0] <- log(pmax(.Machine$double.eps,
        pnorm((Zg[idx0] - rho * Xb[idx0]) / sqrt(1 - rho^2)) * pnorm(-Xb[idx0])))
    }
    ll[!selected] <- pnorm(-Zg[!selected], log.p = TRUE)
    sum(ll)
  }

  probit_sel <- glm(y_sel ~ X_sel - 1, family = binomial(link = "probit"))
  probit_out <- suppressWarnings(glm(y_out[selected] ~ X_out[selected,] - 1, family = binomial(link = "probit")))
  start <- c(coef(probit_out), coef(probit_sel), 0)

  opt <- maxLik::maxLik(loglik, start = start, method = "BFGS", ...)
  coefs <- coef(opt); H <- maxLik::hessian(opt)
  V <- tryCatch(solve(-H), error = function(e) MASS::ginv(-H))
  se <- sqrt(pmax(diag(V), 0))
  ll_val <- as.numeric(logLik(opt))

  list(coefficients = coefs, se = se, z = coefs / se,
       pvalue = 2 * pnorm(-abs(coefs / se)), vcov = V, hessian = H,
       logLik = ll_val, AIC = -2*ll_val + 2*length(coefs),
       BIC = -2*ll_val + log(n)*length(coefs),
       n = n, rho = tanh(coefs[length(coefs)]), convergence = opt$code)
}

# --- Type V: Selection with censored outcome ---
#' @keywords internal
.tobit_type5 <- function(formula, data, left, right, error_dist, ...) {
  # Combine selection (Type II) with censoring (Type I)
  F <- Formula::Formula(formula)
  parsed <- .heckman_parse_formula(F, data)
  y_sel <- parsed$y_selection; X_sel <- parsed$X_selection
  y_out <- parsed$y_outcome; X_out <- parsed$X_outcome
  n <- length(y_sel); k_out <- ncol(X_out); k_sel <- ncol(X_sel)
  selected <- y_sel == 1

  loglik <- function(par) {
    beta <- par[1:k_out]; gamma <- par[(k_out+1):(k_out+k_sel)]
    log_sigma <- par[k_out + k_sel + 1]; rho_raw <- par[k_out + k_sel + 2]
    sigma <- exp(log_sigma); rho <- tanh(rho_raw)
    Xb <- as.numeric(X_out %*% beta); Zg <- as.numeric(X_sel %*% gamma)
    ll <- numeric(n)

    idx1 <- selected
    if (any(idx1)) {
      y_obs <- y_out[idx1]; cens_l <- y_obs <= left & is.finite(left)
      cens_r <- y_obs >= right & is.finite(right); unc <- !cens_l & !cens_r
      mu_i <- Xb[idx1]; zg_i <- Zg[idx1]

      if (any(unc)) {
        resid_i <- (y_obs[unc] - mu_i[unc]) / sigma
        arg <- (zg_i[unc] + rho * resid_i) / sqrt(1 - rho^2)
        ll[idx1][unc] <- dnorm(resid_i, log = TRUE) - log(sigma) + pnorm(arg, log.p = TRUE)
      }
      if (any(cens_l)) {
        ll[idx1][cens_l] <- pnorm((left - mu_i[cens_l]) / sigma, log.p = TRUE) +
          pnorm(zg_i[cens_l], log.p = TRUE)
      }
      if (any(cens_r)) {
        ll[idx1][cens_r] <- pnorm((right - mu_i[cens_r]) / sigma, lower.tail = FALSE, log.p = TRUE) +
          pnorm(zg_i[cens_r], log.p = TRUE)
      }
    }
    ll[!selected] <- pnorm(-Zg[!selected], log.p = TRUE)
    sum(ll)
  }

  probit_sel <- glm(y_sel ~ X_sel - 1, family = binomial(link = "probit"))
  y_sel_obs <- y_out[selected]; X_sel_obs <- X_out[selected, , drop = FALSE]
  ols <- lm.fit(X_sel_obs, y_sel_obs)
  start <- c(ols$coefficients, coef(probit_sel), log(sd(ols$residuals)), 0)

  opt <- maxLik::maxLik(loglik, start = start, method = "BFGS", ...)
  coefs <- coef(opt); H <- maxLik::hessian(opt)
  V <- tryCatch(solve(-H), error = function(e) MASS::ginv(-H))
  se <- sqrt(pmax(diag(V), 0))
  ll_val <- as.numeric(logLik(opt))

  list(coefficients = coefs, se = se, z = coefs / se,
       pvalue = 2 * pnorm(-abs(coefs / se)), vcov = V, hessian = H,
       logLik = ll_val, AIC = -2*ll_val + 2*length(coefs),
       BIC = -2*ll_val + log(n)*length(coefs),
       n = n, convergence = opt$code)
}

#' @export
predict.tobit_flex <- function(object, newdata = NULL, type = "response", ...) {
  if (is.null(newdata)) newdata <- object$model_data
  F <- object$formula; X <- model.matrix(F, data = newdata, rhs = 1)
  k <- ncol(X); beta <- object$coefficients[1:k]
  yhat <- as.numeric(X %*% beta)
  if (type == "response" && object$tobit_type == 1) {
    yhat <- pmax(pmin(yhat, object$right), object$left)
  }
  yhat
}
