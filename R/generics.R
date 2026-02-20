# =============================================================================
# ecoflex: S3 Generic Methods
# =============================================================================

# --- summary ---

#' Summary method for ecoflex models
#'
#' @param object An ecoflex model object.
#' @param vcov Optional: a variance-covariance matrix or a function that
#'   computes one from the object. If supplied, standard errors are recalculated.
#' @param latex If \code{TRUE}, prints a LaTeX table instead of console output.
#'   Default: \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link{to_latex}} when
#'   \code{latex = TRUE}.
#' @return An object of class \code{"summary.ecoflex"} with a coefficient table
#'   (invisible when \code{latex = TRUE}).
#' @export
summary.ecoflex <- function(object, vcov = NULL, latex = FALSE, ...) {
  if (!is.null(vcov)) {
    V <- if (is.function(vcov)) vcov(object) else vcov
    object$se <- sqrt(diag(V))
    object$z <- object$coefficients / object$se
    object$pvalue <- 2 * pnorm(-abs(object$z))
  }

  if (latex) {
    to_latex(object, ...)
    return(invisible(object))
  }

  coef_table <- cbind(
    Estimate   = object$coefficients,
    `Std. Error` = object$se,
    `z value`    = object$z,
    `Pr(>|z|)`   = object$pvalue
  )
  object$coef_table <- coef_table
  class(object) <- c("summary.ecoflex", class(object))
  object
}

# --- print.summary ---

#' Print method for summary.ecoflex objects
#'
#' @param x A summary.ecoflex object
#' @param digits Number of significant digits
#' @param ... Additional arguments (ignored)
#' @export
print.summary.ecoflex <- function(x, digits = 4, ...) {
  cat("\n", x$model_name, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  printCoefmat(x$coef_table, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
  if (!is.null(x$logLik) && is.finite(x$logLik)) {
    cat("\nLog-Likelihood:", format(x$logLik, digits = digits), "\n")
  }
  if (!is.null(x$AIC) && is.finite(x$AIC)) {
    cat("AIC:", format(x$AIC, digits = digits),
        " BIC:", format(x$BIC, digits = digits), "\n")
  }
  if (!is.null(x$n)) {
    cat("N:", x$n, "\n")
  }
  cat(paste(rep("-", 60), collapse = ""), "\n")
  invisible(x)
}

# --- vcov ---

#' Variance-covariance matrix for ecoflex models
#'
#' @param object An ecoflex model object
#' @param type Type of vcov: "standard", "robust", "cluster", "bootstrap"
#' @param cluster Cluster variable (required if type = "cluster")
#' @param R Number of bootstrap replications (if type = "bootstrap")
#' @param ... Additional arguments
#' @export
vcov.ecoflex <- function(object, type = c("standard", "robust", "cluster", "bootstrap"),
                           cluster = NULL, R = 500, ...) {
  type <- match.arg(type)
  switch(type,
    standard = object$vcov,
    robust = {
      if (!is.null(object$hessian) && !is.null(object$scores)) {
        H_inv <- tryCatch(solve(-object$hessian),
                          error = function(e) MASS::ginv(-object$hessian))
        S <- crossprod(object$scores)
        H_inv %*% S %*% H_inv
      } else if (!is.null(object$bread) && !is.null(object$estfun)) {
        sandwich::vcovHC(object, type = "HC1")
      } else {
        warning("Robust vcov not available; returning standard vcov.")
        object$vcov
      }
    },
    cluster = {
      if (is.null(cluster)) stop("Specify 'cluster' variable")
      if (!is.null(object$hessian)) {
        scores <- .compute_scores(object)
        H_inv <- tryCatch(solve(-object$hessian),
                          error = function(e) MASS::ginv(-object$hessian))
        cluster_ids <- unique(cluster)
        G <- length(cluster_ids)
        n <- nrow(scores)
        k <- ncol(scores)
        meat <- Reduce("+", lapply(cluster_ids, function(g) {
          idx <- which(cluster == g)
          s_g <- colSums(scores[idx, , drop = FALSE])
          tcrossprod(s_g)
        }))
        correction <- G / (G - 1) * (n - 1) / (n - k)
        H_inv %*% (correction * meat) %*% H_inv
      } else {
        warning("Cluster vcov not available; returning standard vcov.")
        object$vcov
      }
    },
    bootstrap = .bootstrap_vcov(object, R = R, ...)
  )
}

# --- predict ---

#' Predict method for ecoflex models (generic dispatch)
#'
#' @param object An ecoflex model object
#' @param newdata Optional new data frame for predictions
#' @param type Type of prediction ("response", "link", "prob")
#' @param se.fit If TRUE, return standard errors of predictions
#' @param ... Additional arguments
#' @export
predict.ecoflex <- function(object, newdata = NULL, type = "response",
                              se.fit = FALSE, ...) {
  if (is.null(newdata)) {
    return(object$fitted.values)
  }
  # Subclass-specific methods should override this
  warning("predict not yet implemented for class ", class(object)[1])
  NULL
}

# --- coef ---

#' @export
coef.ecoflex <- function(object, ...) {
  object$coefficients
}

# --- logLik ---

#' @export
logLik.ecoflex <- function(object, ...) {
  val <- object$logLik
  attr(val, "df") <- length(object$coefficients)
  attr(val, "nobs") <- object$n
  class(val) <- "logLik"
  val
}

# --- nobs ---

#' @export
nobs.ecoflex <- function(object, ...) {
  object$n
}
