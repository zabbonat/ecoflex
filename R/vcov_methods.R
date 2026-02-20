# =============================================================================
# ecoflex: Flexible Variance-Covariance Methods
# =============================================================================

#' Flexible variance-covariance matrix for ecoflex models
#'
#' Compute variance-covariance matrices using different methods: inverse Hessian,
#' outer product of gradients (OPG), Huber-White sandwich, clustered, or bootstrap.
#'
#' @param object An ecoflex model object
#' @param type Type of matrix: \code{"hessian"}, \code{"opg"}, \code{"sandwich"},
#'   \code{"cluster"}, \code{"bootstrap"}
#' @param cluster Cluster variable (vector). Required if \code{type = "cluster"}.
#' @param R Number of bootstrap replications. Default: 500.
#' @param parallel Use parallel computation for bootstrap. Default: FALSE.
#' @param cores Number of cores for parallelisation. Default: 2.
#'
#' @return A k x k variance-covariance matrix
#' @export
#'
#' @examples
#' \dontrun{
#' m <- hurdle_flex(y ~ x1 + x2, data = df)
#' V_sand <- flex_vcov(m, type = "sandwich")
#' V_cl   <- flex_vcov(m, type = "cluster", cluster = df$group)
#' V_boot <- flex_vcov(m, type = "bootstrap", R = 1000)
#'
#' # Use with summary:
#' summary(m, vcov = V_cl)
#' }
flex_vcov <- function(object, type = "sandwich", cluster = NULL,
                      R = 500, parallel = FALSE, cores = 2) {

  type <- match.arg(type, c("hessian", "opg", "sandwich", "cluster", "bootstrap"))

  switch(type,
    hessian = {
      H <- object$hessian
      if (is.null(H)) stop("Hessian not available in model object.")
      tryCatch(solve(-H), error = function(e) MASS::ginv(-H))
    },

    opg = {
      scores <- .compute_scores(object)
      S <- crossprod(scores)
      tryCatch(solve(S), error = function(e) MASS::ginv(S))
    },

    sandwich = {
      H <- object$hessian
      if (is.null(H)) stop("Hessian not available in model object.")
      scores <- .compute_scores(object)
      H_inv <- tryCatch(solve(-H), error = function(e) MASS::ginv(-H))
      S <- crossprod(scores)
      H_inv %*% S %*% H_inv
    },

    cluster = {
      if (is.null(cluster)) stop("Specify cluster variable.")
      H <- object$hessian
      if (is.null(H)) stop("Hessian not available in model object.")
      scores <- .compute_scores(object)
      H_inv <- tryCatch(solve(-H), error = function(e) MASS::ginv(-H))

      cluster_ids <- unique(cluster)
      G <- length(cluster_ids)
      n <- nrow(scores)
      k <- ncol(scores)

      meat <- Reduce("+", lapply(cluster_ids, function(g) {
        idx <- which(cluster == g)
        s_g <- colSums(scores[idx, , drop = FALSE])
        tcrossprod(s_g)
      }))

      # Small-cluster correction
      correction <- G / (G - 1) * (n - 1) / (n - k)
      H_inv %*% (correction * meat) %*% H_inv
    },

    bootstrap = {
      .bootstrap_vcov(object, R = R, parallel = parallel, cores = cores)
    }
  )
}

#' @keywords internal
#' Bootstrap variance-covariance for any ecoflex model
.bootstrap_vcov <- function(object, R = 500, parallel = FALSE, cores = 2) {

  boot_fn <- function(data, indices) {
    d <- data[indices, ]
    fit <- tryCatch(
      update(object, data = d),
      error = function(e) NULL
    )
    if (is.null(fit)) return(rep(NA, length(coef(object))))
    coef(fit)
  }

  original_data <- object$model_data
  n <- nrow(original_data)
  k <- length(coef(object))

  if (parallel) {
    if (requireNamespace("parallel", quietly = TRUE)) {
      results <- parallel::mclapply(seq_len(R), function(i) {
        idx <- sample(n, replace = TRUE)
        boot_fn(original_data, idx)
      }, mc.cores = cores)
      boot_coefs <- do.call(rbind, results)
    } else {
      warning("Package 'parallel' not available; using sequential bootstrap.")
      parallel <- FALSE
    }
  }

  if (!parallel) {
    boot_coefs <- matrix(NA_real_, R, k)
    for (i in seq_len(R)) {
      idx <- sample(n, replace = TRUE)
      boot_coefs[i, ] <- boot_fn(original_data, idx)
    }
  }

  # Remove failed replications
  boot_coefs <- boot_coefs[complete.cases(boot_coefs), , drop = FALSE]

  if (nrow(boot_coefs) < 10) {
    warning("Too few successful bootstrap replications (", nrow(boot_coefs),
            "). Variance-covariance may be unreliable.")
  }

  cov(boot_coefs)
}
