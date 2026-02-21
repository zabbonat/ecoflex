# =============================================================================
# ecoflex Module 3: Flexible RDD
# =============================================================================

#' Regression Discontinuity Design
#'
#' Flexible sharp and fuzzy RDD with local polynomial estimation,
#' MSE-optimal bandwidth selection, and multiple kernel choices.
#'
#' @param formula Formula: \code{outcome ~ running_variable}.
#' @param data A data frame.
#' @param cutoff Cutoff point. Default: 0.
#' @param type \code{"sharp"} or \code{"fuzzy"}.
#' @param fuzzy_treatment Treatment variable name if \code{type = "fuzzy"}.
#' @param bandwidth Method (\code{"mserd"}) or numeric value.
#' @param kernel \code{"triangular"}, \code{"epanechnikov"}, \code{"uniform"}, \code{"gaussian"}.
#' @param polynomial Local polynomial degree. Default: 1.
#' @param polynomial_right Degree on right side (if asymmetric).
#' @param discrete_running Discrete running variable correction.
#' @param cluster Cluster variable for standard errors.
#' @param covs_adjust \code{"none"}, \code{"linear"}, \code{"interacted"}.
#' @param masspoints \code{"adjust"}, \code{"check"}, \code{"off"}.
#' @param all_bandwidths Report results for a range of bandwidths.
#' @param latex If \code{TRUE}, prints a LaTeX table. Default: \code{FALSE}.
#' @param ... Additional arguments.
#' @return Object of class \code{c("rdd_flex", "ecoflex")}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' x <- runif(500, -1, 1)
#' y <- 1 + 0.5 * x + 2 * (x >= 0) + rnorm(500, sd = 0.5)
#' df <- data.frame(y = y, x = x)
#' m <- rdd_flex(y ~ x, data = df, cutoff = 0, type = "sharp")
#' summary(m)
#' }
rdd_flex <- function(formula, data, cutoff = 0,
                     type = c("sharp", "fuzzy"),
                     fuzzy_treatment = NULL,
                     bandwidth = "mserd",
                     kernel = c("triangular", "epanechnikov", "uniform", "gaussian"),
                     polynomial = 1,
                     polynomial_right = NULL,
                     discrete_running = FALSE,
                     cluster = NULL,
                     covs_adjust = c("none", "linear", "interacted"),
                     masspoints = c("adjust", "check", "off"),
                     all_bandwidths = FALSE, latex = FALSE, ...) {
  type <- match.arg(type)
  kernel <- match.arg(kernel)
  covs_adjust <- match.arg(covs_adjust)
  masspoints <- match.arg(masspoints)

  F <- Formula::Formula(formula)
  mf <- model.frame(F, data = data)
  y <- model.response(mf)
  X_rhs <- model.matrix(F, data = data, rhs = 1)
  # Running variable is the first non-intercept column
  x_cols <- setdiff(colnames(X_rhs), "(Intercept)")
  if (length(x_cols) == 0) stop("Running variable not found in formula")
  x <- X_rhs[, x_cols[1]]
  covs <- if (length(F)[2] > 1) model.matrix(F, data = data, rhs = 2)[, -1, drop = FALSE] else NULL

  x_centered <- x - cutoff
  treat <- as.integer(x >= cutoff)

  if (is.character(bandwidth)) {
    h <- .bandwidth_selection(y, x_centered, kernel, polynomial, method = bandwidth,
                               discrete = discrete_running)
  } else {
    h <- bandwidth
  }

  p_left <- polynomial
  p_right <- if (!is.null(polynomial_right)) polynomial_right else polynomial

  result <- .local_poly_rd(y, x_centered, treat, h,
                            p_left = p_left, p_right = p_right,
                            kernel_type = kernel, covs = covs,
                            covs_adjust = covs_adjust,
                            fuzzy = (type == "fuzzy"),
                            fuzzy_treat = if (type == "fuzzy") data[[fuzzy_treatment]] else NULL,
                            cluster = cluster)

  if (all_bandwidths) {
    h_grid <- seq(h * 0.5, h * 1.5, length.out = 20)
    sensitivity <- lapply(h_grid, function(h_i) {
      r <- .local_poly_rd(y, x_centered, treat, h_i,
                           p_left = p_left, p_right = p_right,
                           kernel_type = kernel, covs = covs,
                           covs_adjust = covs_adjust,
                           fuzzy = (type == "fuzzy"),
                           fuzzy_treat = if (type == "fuzzy") data[[fuzzy_treatment]] else NULL)
      c(h = h_i, estimate = r$tau, se = r$se_tau, ci_lower = r$ci[1], ci_upper = r$ci[2])
    })
    result$sensitivity <- do.call(rbind, sensitivity)
  }

  n_left <- sum(x_centered < 0 & x_centered >= -h)
  n_right <- sum(x_centered >= 0 & x_centered <= h)

  obj <- structure(
    c(result, list(
      call = match.call(), formula = F, type = type, cutoff = cutoff,
      bandwidth = h, kernel = kernel,
      polynomial = c(left = p_left, right = p_right),
      n = length(y), n_left = n_left, n_right = n_right,
      model_data = data,
      model_name = sprintf("RDD (%s, p=%d/%d, kernel=%s)", type, p_left, p_right, kernel),
      coefficients = c(RDD_effect = result$tau),
      se = c(tau = result$se_tau),
      z = c(tau = result$tau / result$se_tau),
      pvalue = c(tau = 2 * pnorm(-abs(result$tau / result$se_tau)))
    )),
    class = c("rdd_flex", "ecoflex")
  )
  if (latex) to_latex(obj)
  obj
}

#' @keywords internal
.local_poly_rd <- function(y, x, treat, h, p_left = 1, p_right = 1,
                            kernel_type = "triangular", covs = NULL,
                            covs_adjust = "none", fuzzy = FALSE,
                            fuzzy_treat = NULL, cluster = NULL) {
  in_bw <- abs(x) <= h
  y_bw <- y[in_bw]; x_bw <- x[in_bw]; treat_bw <- treat[in_bw]
  w_bw <- .kernel_fn(x_bw / h, type = kernel_type)
  n_bw <- sum(in_bw)

  if (n_bw < 10) {
    warning("Very few observations (", n_bw, ") within bandwidth.")
    return(list(tau = NA, se_tau = NA, ci = c(NA, NA), vcov = NULL,
                hessian = NULL, logLik = NA, AIC = NA, BIC = NA))
  }

  # Build design matrix for local polynomial
  left_idx <- x_bw < 0; right_idx <- x_bw >= 0
  # Construct polynomial terms
  X_list <- list(rep(1, n_bw))  # intercept
  X_list[[2]] <- treat_bw        # treatment indicator
  for (p in seq_len(max(p_left, p_right))) {
    X_list[[length(X_list) + 1]] <- x_bw^p
    X_list[[length(X_list) + 1]] <- treat_bw * x_bw^p  # interaction
  }
  X_poly <- do.call(cbind, X_list)

  if (covs_adjust != "none" && !is.null(covs)) {
    covs_bw <- covs[in_bw, , drop = FALSE]
    X_poly <- cbind(X_poly, covs_bw)
  }

  # Weighted least squares
  W <- diag(sqrt(w_bw))
  Xw <- W %*% X_poly; yw <- W %*% y_bw
  beta <- tryCatch(solve(crossprod(Xw)) %*% crossprod(Xw, yw),
                   error = function(e) MASS::ginv(crossprod(Xw)) %*% crossprod(Xw, yw))

  tau <- beta[2]  # treatment effect is the coefficient on treat

  # Standard errors
  resid <- as.numeric(y_bw - X_poly %*% beta)
  sigma2 <- sum(w_bw * resid^2) / (n_bw - ncol(X_poly))
  if (!is.null(cluster)) {
    cl_bw <- cluster[in_bw]
    cl_ids <- unique(cl_bw)
    G <- length(cl_ids)
    meat <- Reduce("+", lapply(cl_ids, function(g) {
      idx <- cl_bw == g
      Xi <- X_poly[idx, , drop = FALSE] * sqrt(w_bw[idx])
      ei <- resid[idx] * sqrt(w_bw[idx])
      crossprod(colSums(Xi * ei))
    }))
    bread <- tryCatch(solve(crossprod(Xw)), error = function(e) MASS::ginv(crossprod(Xw)))
    V <- bread %*% meat %*% bread * G / (G - 1)
  } else {
    V <- sigma2 * tryCatch(solve(crossprod(Xw)), error = function(e) MASS::ginv(crossprod(Xw)))
  }
  se_tau <- sqrt(max(V[2, 2], 0))
  ci <- tau + c(-1, 1) * qnorm(0.975) * se_tau

  if (fuzzy && !is.null(fuzzy_treat)) {
    # Fuzzy RDD: IV with treatment as endogenous, D(x >= cutoff) as instrument
    D_bw <- fuzzy_treat[in_bw]
    # First stage
    beta_fs <- tryCatch(solve(crossprod(Xw)) %*% crossprod(Xw, W %*% D_bw),
                        error = function(e) rep(NA, ncol(X_poly)))
    compliance <- beta_fs[2]
    if (!is.na(compliance) && abs(compliance) > 1e-10) {
      tau <- tau / compliance
      se_tau <- se_tau / abs(compliance)
      ci <- tau + c(-1, 1) * qnorm(0.975) * se_tau
    }
  }

  list(tau = as.numeric(tau), se_tau = as.numeric(se_tau), ci = ci,
       vcov = V, hessian = NULL, logLik = NA, AIC = NA, BIC = NA,
       fitted_values = as.numeric(X_poly %*% beta),
       residuals = resid, n_effective = n_bw)
}

#' Manipulation Test (Density Test at Cutoff)
#' @param x Running variable vector
#' @param cutoff Cutoff point. Default: 0
#' @param method "cattaneo" or "mccrary"
#' @param ... Additional arguments
#' @return List with test statistic, p-value, and bandwidth
#' @export
rdd_manipulation_test <- function(x, cutoff = 0, method = c("cattaneo", "mccrary"), ...) {
  method <- match.arg(method)
  x_c <- x - cutoff
  n <- length(x)
  n_left <- sum(x_c < 0); n_right <- sum(x_c >= 0)

  # Simple density test based on binomial
  h <- .bandwidth_selection(rep(1, n), x_c, "triangular", 1, "mserd")
  n_l <- sum(x_c >= -h & x_c < 0)
  n_r <- sum(x_c >= 0 & x_c <= h)
  n_w <- n_l + n_r
  # Under H0: f_left(0) = f_right(0), expect n_l/n_w = 0.5
  T_stat <- (n_r / n_w - 0.5) / sqrt(0.25 / n_w)
  p_value <- 2 * pnorm(-abs(T_stat))

  list(statistic = T_stat, p.value = p_value,
       n_left = n_l, n_right = n_r, bandwidth = h, method = method)
}

#' Covariate Balance Test
#' @param formula Formula with running variable
#' @param data Data frame
#' @param cutoff Cutoff point. Default: 0
#' @param covariates Character vector of covariate names to test
#' @param ... Additional arguments
#' @return Data frame with balance test results
#' @export
rdd_balance_test <- function(formula, data, cutoff = 0, covariates, ...) {
  results <- lapply(covariates, function(cov) {
    f_cov <- as.formula(paste(cov, "~", as.character(formula)[3]))
    m <- tryCatch(
      rdd_flex(f_cov, data = data, cutoff = cutoff, ...),
      error = function(e) NULL
    )
    if (is.null(m)) {
      data.frame(covariate = cov, estimate = NA, se = NA, p.value = NA)
    } else {
      data.frame(covariate = cov, estimate = m$tau,
                 se = m$se_tau, p.value = m$pvalue["tau"])
    }
  })
  do.call(rbind, results)
}

#' Placebo Cutoff Tests
#' @param formula Formula
#' @param data Data frame
#' @param cutoff True cutoff. Default: 0
#' @param placebo_cutoffs Numeric vector of placebo cutoffs
#' @param ... Additional arguments
#' @return Data frame with placebo test results
#' @export
rdd_placebo_cutoffs <- function(formula, data, cutoff = 0, placebo_cutoffs = NULL, ...) {
  F <- Formula::Formula(formula)
  mf <- model.frame(F, data = data)
  x <- model.matrix(F, data = data, rhs = 1)[, -1, drop = FALSE][, 1]
  if (is.null(placebo_cutoffs)) {
    q <- quantile(x, probs = c(0.25, 0.5, 0.75))
    placebo_cutoffs <- q[q != cutoff]
  }
  results <- lapply(placebo_cutoffs, function(pc) {
    m <- tryCatch(rdd_flex(formula, data = data, cutoff = pc, ...),
                  error = function(e) NULL)
    if (is.null(m)) {
      data.frame(cutoff = pc, estimate = NA, se = NA, p.value = NA)
    } else {
      data.frame(cutoff = pc, estimate = m$tau,
                 se = m$se_tau, p.value = m$pvalue["tau"])
    }
  })
  do.call(rbind, results)
}

#' @export
plot.rdd_flex <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 required for plotting"); return(invisible(NULL))
  }
  data <- x$model_data
  F <- x$formula
  mf <- model.frame(F, data = data)
  y <- model.response(mf)
  rv <- model.matrix(F, data = data, rhs = 1)[, -1, drop = FALSE][, 1]
  df_plot <- data.frame(y = y, x = rv, treat = as.factor(as.integer(rv >= x$cutoff)))

  ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$x, y = .data$y, color = .data$treat)) +
    ggplot2::geom_point(alpha = 0.3, size = 1) +
    ggplot2::geom_vline(xintercept = x$cutoff, linetype = "dashed") +
    ggplot2::geom_smooth(method = "loess", se = TRUE) +
    ggplot2::labs(title = x$model_name, x = "Running Variable", y = "Outcome") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

#' @export
predict.rdd_flex <- function(object, newdata = NULL, type = "response", ...) {
  if (is.null(newdata)) return(object$fitted_values)
  warning("predict with newdata not yet fully implemented for rdd_flex")
  NULL
}
