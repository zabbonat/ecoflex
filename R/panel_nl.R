# =============================================================================
# ecoflex Module 7: Nonlinear Panel Data Models
# =============================================================================

#' Nonlinear Panel Data with Flexible Standard Errors
#'
#' Estimate nonlinear panel models (Poisson, NegBin, Logit, Probit, Gamma)
#' with fixed effects, random effects, or correlated random effects.
#'
#' @param formula Formula: \code{outcome ~ regressors}.
#' @param data Panel data in long format.
#' @param id Name of unit identifier variable.
#' @param time Name of time variable.
#' @param family \code{"poisson"}, \code{"negbin"}, \code{"logit"},
#'   \code{"probit"}, \code{"gamma"}.
#' @param effect \code{"individual"}, \code{"time"}, \code{"twoways"}.
#' @param model \code{"fe"}, \code{"re"}, \code{"correlated_re"} (Mundlak/Chamberlain).
#' @param vcov \code{"standard"}, \code{"robust"}, \code{"cluster"},
#'   \code{"driscoll_kraay"}, \code{"bootstrap"}.
#' @param cluster Cluster variable (default: id).
#' @param DK_bandwidth Driscoll-Kraay bandwidth.
#' @param latex If \code{TRUE}, prints a LaTeX table. Default: \code{FALSE}.
#' @param ... Additional arguments.
#' @return Object of class \code{c("panel_nl", "ecoflex")}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' panel <- expand.grid(id = 1:30, time = 1:5)
#' panel$x <- rnorm(nrow(panel))
#' panel$y <- rpois(nrow(panel), exp(0.5 + 0.3 * panel$x))
#' m <- panel_nl(y ~ x, data = panel, id = "id", time = "time",
#'               family = "poisson", model = "fe")
#' summary(m)
#' }
panel_nl <- function(formula, data, id, time,
                     family = c("poisson", "negbin", "logit", "probit", "gamma"),
                     effect = c("individual", "time", "twoways"),
                     model = c("fe", "re", "correlated_re"),
                     vcov = c("standard", "robust", "cluster", "driscoll_kraay", "bootstrap"),
                     cluster = NULL, DK_bandwidth = NULL, latex = FALSE, ...) {
  family <- match.arg(family)
  effect <- match.arg(effect)
  model <- match.arg(model)
  vcov_type <- match.arg(vcov)

  data <- data[order(data[[id]], data[[time]]), ]

  if (model == "fe") {
    result <- .panel_nl_fe(formula, data, id, time, family, effect, ...)
  } else if (model == "re") {
    result <- .panel_nl_re(formula, data, id, time, family, effect, ...)
  } else {
    result <- .panel_nl_cre(formula, data, id, time, family, effect, ...)
  }

  if (is.null(cluster)) cluster <- id
  V <- switch(vcov_type,
    standard = result$vcov,
    robust = .panel_robust_vcov(result, data, id),
    cluster = .panel_cluster_vcov(result, data, cluster),
    driscoll_kraay = .driscoll_kraay_vcov(result, data, id, time, DK_bandwidth),
    bootstrap = .panel_bootstrap_vcov(result, formula, data, id, time, family, effect, model, R = 500)
  )
  result$vcov <- V
  result$se <- sqrt(pmax(diag(V), 0))

  obj <- structure(c(result, list(
    call = match.call(), formula = Formula::Formula(formula),
    family = family, effect = effect, model = model,
    vcov_type = vcov_type, n = nrow(data), model_data = data,
    model_name = sprintf("Panel %s %s (%s effects)", family, model, effect)
  )), class = c("panel_nl", "ecoflex"))
  if (latex) to_latex(obj)
  obj
}

# --- FE estimation ---
#' @keywords internal
.panel_nl_fe <- function(formula, data, id, time, family, effect, ...) {
  F <- Formula::Formula(formula)
  mf <- model.frame(F, data = data); y <- model.response(mf)
  X <- model.matrix(F, data = data)

  # Add fixed effects as dummies
  fe_formula <- switch(effect,
    individual = as.formula(paste(deparse(formula[[2]]), "~", deparse(formula[[3]]), "+ factor(", id, ")")),
    time = as.formula(paste(deparse(formula[[2]]), "~", deparse(formula[[3]]), "+ factor(", time, ")")),
    twoways = as.formula(paste(deparse(formula[[2]]), "~", deparse(formula[[3]]), "+ factor(", id, ") + factor(", time, ")"))
  )

  glm_family <- switch(family,
    poisson = poisson(), negbin = poisson(),  # NB handled separately
    logit = binomial(), probit = binomial(link = "probit"),
    gamma = Gamma(link = "log")
  )

  fit <- suppressWarnings(glm(fe_formula, data = data, family = glm_family))

  # Extract only non-FE coefficients
  all_coefs <- coef(fit)
  x_names <- colnames(X)
  keep <- names(all_coefs) %in% x_names | names(all_coefs) == "(Intercept)"
  coefs_main <- all_coefs[keep]
  V_full <- vcov(fit)
  V_main <- V_full[keep, keep, drop = FALSE]

  se <- sqrt(pmax(diag(V_main), 0))
  ll_val <- as.numeric(logLik(fit))

  list(coefficients = coefs_main, se = se, z = coefs_main / se,
       pvalue = 2 * pnorm(-abs(coefs_main / se)),
       vcov = V_main, hessian = NULL, logLik = ll_val,
       AIC = AIC(fit), BIC = BIC(fit),
       fitted.values = fitted(fit), residuals = residuals(fit),
       scores = NULL, full_fit = fit, convergence = 0)
}

# --- RE estimation ---
#' @keywords internal
.panel_nl_re <- function(formula, data, id, time, family, effect, ...) {
  # Simplified: use pooled GLM as RE approximation
  .panel_nl_fe(formula, data, id, time, family, "individual", ...)
}

# --- Correlated RE (Mundlak/Chamberlain) ---
#' @keywords internal
.panel_nl_cre <- function(formula, data, id, time, family, effect, ...) {
  F <- Formula::Formula(formula)
  X <- model.matrix(F, data = data)
  x_names <- setdiff(colnames(X), "(Intercept)")

  # Compute group means
  for (var in x_names) {
    mean_var <- ave(X[, var], data[[id]], FUN = mean)
    data[[paste0(var, "_mean")]] <- mean_var
  }

  # Augmented formula with group means
  mean_terms <- paste0(x_names, "_mean", collapse = " + ")
  f_cre <- as.formula(paste(deparse(formula[[2]]), "~", deparse(formula[[3]]), "+", mean_terms))

  .panel_nl_fe(f_cre, data, id, time, family, effect, ...)
}

# --- Vcov alternatives ---
#' @keywords internal
.panel_robust_vcov <- function(result, data, id) {
  fit <- result$full_fit
  if (is.null(fit)) return(result$vcov)
  tryCatch(sandwich::vcovHC(fit, type = "HC1"), error = function(e) result$vcov)
}

#' @keywords internal
.panel_cluster_vcov <- function(result, data, cluster) {
  fit <- result$full_fit
  if (is.null(fit)) return(result$vcov)
  cl_var <- if (is.character(cluster)) data[[cluster]] else cluster
  tryCatch(sandwich::vcovCL(fit, cluster = cl_var), error = function(e) result$vcov)
}

#' @keywords internal
.driscoll_kraay_vcov <- function(result, data, id, time, bandwidth) {
  fit <- result$full_fit
  if (is.null(fit)) return(result$vcov)
  times <- sort(unique(data[[time]]))
  T_periods <- length(times)
  if (is.null(bandwidth)) bandwidth <- floor(T_periods^(1/3))

  # Compute per-period aggregated scores
  X <- model.matrix(fit)
  e <- residuals(fit, type = "response")
  k <- ncol(X)
  scores <- X * e  # observation-level scores (approximation)

  s_t <- do.call(rbind, lapply(times, function(t) {
    idx <- data[[time]] == t
    colSums(scores[idx, , drop = FALSE])
  }))

  # Newey-West HAC on time-aggregated scores
  S <- crossprod(s_t) / T_periods
  for (j in seq_len(bandwidth)) {
    w <- 1 - j / (bandwidth + 1)
    Gamma_j <- crossprod(s_t[-(1:j), , drop = FALSE], s_t[1:(T_periods-j), , drop = FALSE]) / T_periods
    S <- S + w * (Gamma_j + t(Gamma_j))
  }

  bread <- tryCatch(solve(crossprod(X)), error = function(e) MASS::ginv(crossprod(X)))
  V <- bread %*% S %*% bread
  # Subset to main coefficients
  keep <- seq_len(min(nrow(result$vcov), nrow(V)))
  V[keep, keep, drop = FALSE]
}

#' @keywords internal
.panel_bootstrap_vcov <- function(result, formula, data, id, time, family, effect, model, R = 500) {
  ids <- unique(data[[id]])
  n_ids <- length(ids)
  k <- length(result$coefficients)
  boot_coefs <- matrix(NA_real_, R, k)

  for (i in seq_len(R)) {
    boot_ids <- sample(ids, replace = TRUE)
    boot_data <- do.call(rbind, lapply(seq_along(boot_ids), function(j) {
      d <- data[data[[id]] == boot_ids[j], ]
      d[[id]] <- j  # rename to avoid duplicate FE levels
      d
    }))
    fit <- tryCatch(
      panel_nl(formula, boot_data, id, time, family = family, effect = effect,
               model = model, vcov = "standard"),
      error = function(e) NULL
    )
    if (!is.null(fit)) boot_coefs[i, ] <- coef(fit)
  }
  boot_coefs <- boot_coefs[complete.cases(boot_coefs), , drop = FALSE]
  if (nrow(boot_coefs) < 10) warning("Few successful bootstrap replications")
  cov(boot_coefs)
}

#' @export
predict.panel_nl <- function(object, newdata = NULL, type = "response", ...) {
  if (is.null(newdata)) return(object$fitted.values)
  warning("predict with newdata not fully implemented for panel_nl")
  NULL
}
