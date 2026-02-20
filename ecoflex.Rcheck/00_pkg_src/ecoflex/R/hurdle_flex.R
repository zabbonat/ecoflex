# =============================================================================
# ecoflex Module 1: Flexible Hurdle Models
# =============================================================================

#' Hurdle Model with Arbitrary Threshold
#'
#' Estimate a hurdle model where the threshold can be set at any non-negative
#' integer value (default: 0). Supports Poisson, Negative Binomial, and
#' Geometric count distributions with Logit or Probit binary components.
#'
#' @param formula Multi-part formula: \code{count_model | binary_model}.
#'   If only one part is specified, it is used for both components.
#' @param data A data frame.
#' @param threshold Non-negative integer threshold. Default: 0.
#' @param dist Count distribution: \code{"poisson"}, \code{"negbin"},
#'   \code{"geometric"}. Default: \code{"poisson"}.
#' @param zero_dist Binary distribution: \code{"binomial"} (logit),
#'   \code{"probit"}. Default: \code{"binomial"}.
#' @param method Estimation method: \code{"ml"} (joint maximum likelihood),
#'   \code{"two-step"}. Default: \code{"ml"}.
#' @param start Optional starting values.
#' @param ... Additional arguments passed to the optimiser.
#'
#' @return An object of class \code{c("hurdle_flex", "ecoflex")} with
#'   components: coefficients, se, z, pvalue, vcov, hessian, logLik, AIC, BIC,
#'   convergence, plus model data and metadata.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard hurdle (threshold = 0)
#' m1 <- hurdle_flex(y ~ x1 + x2, data = df)
#'
#' # Hurdle with threshold at 10
#' m2 <- hurdle_flex(citations ~ impact_factor + authors |
#'                   impact_factor + journal_rank,
#'                   data = df, threshold = 10, dist = "negbin")
#'
#' # Compare thresholds
#' sapply(c(0, 5, 10, 20), function(t) {
#'   m <- hurdle_flex(y ~ x1 + x2, data = df, threshold = t)
#'   c(threshold = t, AIC = AIC(m), BIC = BIC(m))
#' })
#' }
hurdle_flex <- function(formula, data, threshold = 0,
                        dist = c("poisson", "negbin", "geometric"),
                        zero_dist = c("binomial", "probit"),
                        method = c("ml", "two-step"),
                        start = NULL, latex = FALSE, ...) {

  # Validate inputs
  dist <- match.arg(dist)
  zero_dist <- match.arg(zero_dist)
  method <- match.arg(method)
  stopifnot(is.numeric(threshold), threshold >= 0, threshold == floor(threshold))

  # Parse multi-part formula
  F <- Formula::Formula(formula)
  n_parts <- length(F)[2]

  # Prepare model frames
  mf <- model.frame(F, data = data)
  y <- model.response(mf)

  if (n_parts == 1) {
    X_count <- model.matrix(F, data = data, rhs = 1)
    X_binary <- X_count
  } else {
    X_count <- model.matrix(F, data = data, rhs = 1)
    X_binary <- model.matrix(F, data = data, rhs = 2)
  }

  # Binary variable: exceeds the threshold?
  y_binary <- as.integer(y > threshold)
  # Counts above threshold (for truncated part)
  y_above <- y[y > threshold] - threshold  # shift

  n <- length(y)
  n_above <- sum(y_binary)
  n_below <- n - n_above

  if (method == "two-step") {
    result <- .hurdle_two_step(y_binary, y_above, X_binary, X_count,
                                y > threshold, dist, zero_dist, threshold)
  } else {
    result <- .hurdle_ml(y, y_binary, X_binary, X_count,
                          dist, zero_dist, threshold, start, ...)
  }

  # Build output
  obj <- structure(
    c(result, list(
      call = match.call(),
      formula = F,
      threshold = threshold,
      dist = dist,
      zero_dist = zero_dist,
      method = method,
      n = n,
      n_above = n_above,
      n_below = n_below,
      model_data = data,
      model_name = sprintf("Hurdle Model (threshold = %d, dist = %s)", threshold, dist),
      y = y,
      X_count = X_count,
      X_binary = X_binary
    )),
    class = c("hurdle_flex", "ecoflex")
  )
  if (latex) to_latex(obj)
  obj
}

# =============================================================================
# ML joint estimation
# =============================================================================

#' @keywords internal
.hurdle_ml <- function(y, y_binary, X_binary, X_count,
                        dist, zero_dist, threshold, start, ...) {

  k_binary <- ncol(X_binary)
  k_count <- ncol(X_count)
  n <- length(y)

  # Joint log-likelihood
  loglik <- function(par) {
    beta_b <- par[1:k_binary]
    beta_c <- par[(k_binary + 1):(k_binary + k_count)]
    alpha <- if (dist == "negbin") exp(par[k_binary + k_count + 1]) else NULL

    # Binary component
    eta_b <- as.numeric(X_binary %*% beta_b)
    if (zero_dist == "binomial") {
      p <- plogis(eta_b)
    } else {
      p <- pnorm(eta_b)
    }

    # Count component (left-truncated at 0, applied to y - threshold)
    mu <- exp(as.numeric(X_count %*% beta_c))

    ll <- numeric(n)

    # Observations <= threshold
    idx_below <- y <= threshold
    ll[idx_below] <- log(pmax(1 - p[idx_below], .Machine$double.eps))

    # Observations > threshold
    idx_above <- !idx_below
    y_shifted <- y[idx_above] - threshold
    mu_above <- mu[idx_above]

    # P(Y* > 0 | Y* ~ dist) for truncation
    if (dist == "poisson") {
      log_dens <- dpois(y_shifted, lambda = mu_above, log = TRUE)
      log_prob_positive <- log(pmax(1 - dpois(0, mu_above), .Machine$double.eps))
    } else if (dist == "negbin") {
      log_dens <- dnbinom(y_shifted, mu = mu_above, size = alpha, log = TRUE)
      log_prob_positive <- log(pmax(1 - dnbinom(0, mu = mu_above, size = alpha),
                                     .Machine$double.eps))
    } else {
      # geometric: NB with size = 1
      log_dens <- dnbinom(y_shifted, mu = mu_above, size = 1, log = TRUE)
      log_prob_positive <- log(pmax(1 - dnbinom(0, mu = mu_above, size = 1),
                                     .Machine$double.eps))
    }

    ll[idx_above] <- log(pmax(p[idx_above], .Machine$double.eps)) +
                     log_dens - log_prob_positive

    sum(ll)
  }

  # Starting values
  if (is.null(start)) {
    glm_binary <- suppressWarnings(
      glm(y_binary ~ X_binary - 1,
          family = if (zero_dist == "binomial") binomial() else
            binomial(link = "probit"))
    )
    y_shifted_all <- pmax(y - threshold, 0)
    above_idx <- y_binary == 1
    if (sum(above_idx) > 0) {
      glm_count <- suppressWarnings(
        glm(y_shifted_all[above_idx] ~ X_count[above_idx, ] - 1,
            family = poisson())
      )
      start <- c(coef(glm_binary), coef(glm_count))
    } else {
      start <- c(coef(glm_binary), rep(0, k_count))
    }
    if (dist == "negbin") start <- c(start, log(1))
  }

  # Optimisation
  opt <- maxLik::maxLik(loglik, start = start, method = "BFGS", ...)

  # Extract results
  coefs <- coef(opt)
  names_b <- paste0("binary_", colnames(X_binary))
  names_c <- paste0("count_", colnames(X_count))
  names_all <- c(names_b, names_c)
  if (dist == "negbin") names_all <- c(names_all, "log_alpha")
  names(coefs) <- names_all

  H <- maxLik::hessian(opt)
  V <- tryCatch(solve(-H), error = function(e) {
    warning("Hessian inversion failed; using generalised inverse.")
    MASS::ginv(-H)
  })
  se <- sqrt(pmax(diag(V), 0))
  names(se) <- names_all

  ll_val <- as.numeric(logLik(opt))

  list(
    coefficients = coefs,
    se = se,
    z = coefs / se,
    pvalue = 2 * pnorm(-abs(coefs / se)),
    vcov = V,
    hessian = H,
    logLik = ll_val,
    AIC = -2 * ll_val + 2 * length(coefs),
    BIC = -2 * ll_val + log(n) * length(coefs),
    convergence = opt$code
  )
}

# =============================================================================
# Two-step estimation
# =============================================================================

#' @keywords internal
.hurdle_two_step <- function(y_binary, y_above, X_binary, X_count,
                              above_idx, dist, zero_dist, threshold) {

  n <- length(y_binary)

  # Stage 1: binary model
  fam_binary <- if (zero_dist == "binomial") binomial() else binomial(link = "probit")
  fit_binary <- suppressWarnings(glm(y_binary ~ X_binary - 1, family = fam_binary))

  # Stage 2: truncated count (only observations above threshold)
  X_above <- X_count[above_idx, , drop = FALSE]

  if (dist == "poisson") {
    fit_count <- .truncated_poisson(y_above, X_above)
  } else if (dist == "negbin") {
    fit_count <- .truncated_negbin(y_above, X_above)
  } else {
    fit_count <- .truncated_negbin(y_above, X_above, size_fixed = 1)
  }

  coefs <- c(
    setNames(coef(fit_binary), paste0("binary_", colnames(X_binary))),
    setNames(fit_count$coefficients, paste0("count_", colnames(X_count)))
  )

  # Block-diagonal vcov (independence between stages)
  V_binary <- vcov(fit_binary)
  V_count <- fit_count$vcov
  V <- as.matrix(Matrix::bdiag(V_binary, V_count))

  se <- sqrt(pmax(diag(V), 0))

  ll_val <- as.numeric(logLik(fit_binary)) + as.numeric(fit_count$logLik)

  list(
    coefficients = coefs,
    se = se,
    z = coefs / se,
    pvalue = 2 * pnorm(-abs(coefs / se)),
    vcov = V,
    hessian = NULL,
    logLik = ll_val,
    AIC = -2 * ll_val + 2 * length(coefs),
    BIC = -2 * ll_val + log(n) * length(coefs),
    convergence = 0
  )
}

# =============================================================================
# Predict method for hurdle_flex
# =============================================================================

#' @export
predict.hurdle_flex <- function(object, newdata = NULL, type = "response", ...) {
  if (is.null(newdata)) newdata <- object$model_data

  F <- object$formula
  n_parts <- length(F)[2]

  if (n_parts == 1) {
    X_count <- model.matrix(F, data = newdata, rhs = 1)
    X_binary <- X_count
  } else {
    X_count <- model.matrix(F, data = newdata, rhs = 1)
    X_binary <- model.matrix(F, data = newdata, rhs = 2)
  }

  k_binary <- ncol(X_binary)
  k_count <- ncol(X_count)
  coefs <- object$coefficients

  beta_b <- coefs[1:k_binary]
  beta_c <- coefs[(k_binary + 1):(k_binary + k_count)]

  eta_b <- as.numeric(X_binary %*% beta_b)
  p <- if (object$zero_dist == "binomial") plogis(eta_b) else pnorm(eta_b)

  mu <- exp(as.numeric(X_count %*% beta_c))

  if (type == "response") {
    # E[Y] = P(Y > threshold) * E[Y | Y > threshold]
    # E[Y | Y > threshold] = threshold + E[Y* | Y* > 0]
    if (object$dist == "poisson") {
      e_trunc <- mu / (1 - exp(-mu))
    } else {
      # For negbin/geometric, use approximation
      e_trunc <- mu / (1 - dnbinom(0, mu = mu, size = 1))
    }
    fitted <- p * (object$threshold + e_trunc)
  } else if (type == "prob") {
    fitted <- p
  } else {
    fitted <- list(binary = eta_b, count = log(mu))
  }

  fitted
}
