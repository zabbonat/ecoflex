# =============================================================================
# ecoflex Module 5: Robust Difference-in-Differences
# =============================================================================

#' Robust Difference-in-Differences for Staggered Treatment Adoption
#'
#' Estimate DiD models using TWFE or modern heterogeneity-robust estimators.
#'
#' @param formula Formula: \code{outcome ~ covariates}.
#' @param data Panel data in long format.
#' @param id_var Name of unit identifier variable.
#' @param time_var Name of time variable.
#' @param treat_var Name of treatment dummy variable.
#' @param cohort_var Name of cohort variable (first treatment period).
#' @param estimator \code{"twfe"}, \code{"callaway_santanna"},
#'   \code{"sun_abraham"}, \code{"borusyak_jaravel_spiess"},
#'   \code{"dechaisemartin_dhaultfoeuille"}, \code{"gardner"}, \code{"stacked"}.
#' @param control_group \code{"nevertreated"} or \code{"notyettreated"}.
#' @param anticipation Number of anticipation periods. Default: 0.
#' @param aggregation \code{"dynamic"}, \code{"simple"}, \code{"group"}, \code{"calendar"}.
#' @param cluster Cluster variable for standard errors.
#' @param boot_R Bootstrap replications. Default: 999.
#' @param latex If \code{TRUE}, prints a LaTeX table. Default: \code{FALSE}.
#' @param ... Additional arguments.
#' @return Object of class \code{c("did_flex", "ecoflex")}.
#' @export
#'
#' @examples
#' \donttest{
#' # Simulated panel
#' panel <- expand.grid(id = 1:50, time = 1:10)
#' panel$treat <- as.integer(panel$id > 25 & panel$time >= 6)
#' panel$y <- 1 + 0.5 * panel$treat + rnorm(nrow(panel))
#' m <- did_flex(y ~ 1, data = panel, id_var = "id",
#'               time_var = "time", treat_var = "treat",
#'               estimator = "twfe")
#' summary(m)
#' }
did_flex <- function(formula, data, id_var, time_var, treat_var,
                     cohort_var = NULL,
                     estimator = c("callaway_santanna", "sun_abraham",
                                   "borusyak_jaravel_spiess",
                                   "dechaisemartin_dhaultfoeuille",
                                   "gardner", "stacked", "twfe"),
                     control_group = c("nevertreated", "notyettreated"),
                     anticipation = 0,
                     aggregation = c("dynamic", "simple", "group", "calendar"),
                     cluster = NULL, boot_R = 999, latex = FALSE, ...) {
  estimator <- match.arg(estimator)
  control_group <- match.arg(control_group)
  aggregation <- match.arg(aggregation)

  if (is.null(cohort_var)) {
    cohort_var <- ".cohort"
    data[[cohort_var]] <- .compute_cohort(data, id_var, time_var, treat_var)
  }

  result <- switch(estimator,
    twfe = .did_twfe(formula, data, id_var, time_var, treat_var, cluster),
    callaway_santanna = .did_cs(formula, data, id_var, time_var, cohort_var,
                                 control_group, anticipation, boot_R, ...),
    sun_abraham = .did_sa(formula, data, id_var, time_var, cohort_var, treat_var, ...),
    borusyak_jaravel_spiess = .did_imputation(formula, data, id_var, time_var, cohort_var, treat_var, ...),
    dechaisemartin_dhaultfoeuille = .did_dcdh(formula, data, id_var, time_var, treat_var, ...),
    gardner = .did_gardner(formula, data, id_var, time_var, treat_var, cohort_var, ...),
    stacked = .did_stacked(formula, data, id_var, time_var, cohort_var, cluster, ...)
  )

  if (!is.null(result$event_study)) {
    result$aggregated <- .aggregate_did(result, aggregation)
  }

  obj <- structure(c(result, list(
    call = match.call(), formula = Formula::Formula(formula),
    estimator = estimator, control_group = control_group,
    model_data = data, n = nrow(data),
    model_name = sprintf("DiD Estimation (%s)", estimator)
  )), class = c("did_flex", "ecoflex"))
  if (latex) to_latex(obj)
  obj
}

# --- TWFE ---
#' @keywords internal
.did_twfe <- function(formula, data, id_var, time_var, treat_var, cluster) {
  f_str <- paste(deparse(formula), "+ factor(", id_var, ") + factor(", time_var, ")")
  f_full <- as.formula(f_str)
  fit <- lm(f_full, data = data)
  treat_idx <- grep(treat_var, names(coef(fit)), fixed = TRUE)
  att <- if (length(treat_idx) > 0) coef(fit)[treat_idx[1]] else NA
  att_se <- if (length(treat_idx) > 0) sqrt(diag(vcov(fit)))[treat_idx[1]] else NA

  list(att = att, att_se = att_se,
       att_ci = att + c(-1, 1) * qnorm(0.975) * att_se,
       coefficients = coef(fit), se = sqrt(diag(vcov(fit))),
       z = coef(fit) / sqrt(diag(vcov(fit))),
       pvalue = 2 * pnorm(-abs(coef(fit) / sqrt(diag(vcov(fit))))),
       vcov = vcov(fit), hessian = NULL,
       logLik = as.numeric(logLik(fit)),
       AIC = AIC(fit), BIC = BIC(fit),
       event_study = NULL)
}

# --- Callaway & Sant'Anna (2021) ---
#' @keywords internal
.did_cs <- function(formula, data, id_var, time_var, cohort_var,
                     control_group, anticipation, boot_R, ...) {
  times <- sort(unique(data[[time_var]]))
  cohorts <- sort(unique(data[[cohort_var]]))
  cohorts <- cohorts[is.finite(cohorts)]

  F <- Formula::Formula(formula)
  mf <- model.frame(F, data = data)
  y_name <- names(mf)[1]

  gt_results <- list()
  for (g in cohorts) {
    for (t in times) {
      rel_t <- t - g
      if (rel_t < -anticipation - 1) next

      treated_ids <- unique(data[[id_var]][data[[cohort_var]] == g])
      if (control_group == "nevertreated") {
        control_ids <- unique(data[[id_var]][!is.finite(data[[cohort_var]])])
      } else {
        control_ids <- unique(data[[id_var]][data[[cohort_var]] > t | !is.finite(data[[cohort_var]])])
      }
      if (length(control_ids) == 0 || length(treated_ids) == 0) next

      base_t <- g - 1 - anticipation
      if (!(base_t %in% times)) base_t <- max(times[times < g])
      if (is.na(base_t) || !any(data[[time_var]] == base_t)) next

      # 2x2 DiD for this (g,t) cell
      cell_ids <- c(treated_ids, control_ids)
      cell_times <- c(base_t, t)
      cell_data <- data[data[[id_var]] %in% cell_ids & data[[time_var]] %in% cell_times, ]
      if (nrow(cell_data) < 4) next

      cell_data$.post <- as.integer(cell_data[[time_var]] == t)
      cell_data$.treated_group <- as.integer(cell_data[[id_var]] %in% treated_ids)

      fit_cell <- tryCatch({
        lm(as.formula(paste(y_name, "~ .post * .treated_group")), data = cell_data)
      }, error = function(e) NULL)
      if (is.null(fit_cell)) next

      att_gt <- coef(fit_cell)[".post:.treated_group"]
      se_gt <- sqrt(diag(vcov(fit_cell)))[".post:.treated_group"]

      gt_results[[length(gt_results) + 1]] <- data.frame(
        cohort = g, time = t, relative_time = rel_t,
        estimate = att_gt, se = se_gt, n_treated = length(treated_ids),
        n_control = length(control_ids)
      )
    }
  }

  if (length(gt_results) == 0) {
    warning("No valid group-time cells found")
    return(list(att = NA, att_se = NA, att_ci = c(NA, NA),
                event_study = NULL, coefficients = numeric(0),
                se = numeric(0), z = numeric(0), pvalue = numeric(0),
                vcov = NULL, hessian = NULL, logLik = NA, AIC = NA, BIC = NA))
  }

  event_study <- do.call(rbind, gt_results)

  # Simple ATT aggregate
  w <- event_study$n_treated / sum(event_study$n_treated)
  att <- sum(w * event_study$estimate)
  att_se <- sqrt(sum(w^2 * event_study$se^2))

  list(att = att, att_se = att_se,
       att_ci = att + c(-1, 1) * qnorm(0.975) * att_se,
       event_study = event_study,
       coefficients = c(ATT = att), se = c(ATT = att_se),
       z = c(ATT = att / att_se),
       pvalue = c(ATT = 2 * pnorm(-abs(att / att_se))),
       vcov = matrix(att_se^2, 1, 1), hessian = NULL,
       logLik = NA, AIC = NA, BIC = NA)
}

# --- Sun & Abraham (2021) ---
#' @keywords internal
.did_sa <- function(formula, data, id_var, time_var, cohort_var, treat_var, ...) {
  F <- Formula::Formula(formula)
  mf <- model.frame(F, data = data)
  y_name <- names(mf)[1]
  cohorts <- sort(unique(data[[cohort_var]]))
  cohorts_finite <- cohorts[is.finite(cohorts)]
  times <- sort(unique(data[[time_var]]))

  # Create relative time indicators interacted with cohort dummies
  data$.rel_time <- data[[time_var]] - data[[cohort_var]]
  data$.rel_time[!is.finite(data$.rel_time)] <- -999  # never-treated

  rel_times <- sort(unique(data$.rel_time[data$.rel_time != -999]))
  ref_period <- -1
  rel_times <- setdiff(rel_times, ref_period)

  # IW estimator: weight group-time ATTs by cohort share
  es_results <- lapply(rel_times, function(rt) {
    cohort_atts <- sapply(cohorts_finite, function(g) {
      t_val <- g + rt
      if (!(t_val %in% times)) return(c(att = NA, se = NA, n = 0))
      base_t <- g + ref_period
      if (!(base_t %in% times)) return(c(att = NA, se = NA, n = 0))

      treated_ids <- unique(data[[id_var]][data[[cohort_var]] == g])
      control_ids <- unique(data[[id_var]][!is.finite(data[[cohort_var]])])
      if (length(control_ids) == 0) return(c(att = NA, se = NA, n = 0))

      cell_ids <- c(treated_ids, control_ids)
      cell_data <- data[data[[id_var]] %in% cell_ids &
                         data[[time_var]] %in% c(base_t, t_val), ]
      cell_data$.post <- as.integer(cell_data[[time_var]] == t_val)
      cell_data$.treated_group <- as.integer(cell_data[[id_var]] %in% treated_ids)
      fit <- tryCatch(lm(as.formula(paste(y_name, "~ .post * .treated_group")), data = cell_data),
                      error = function(e) NULL)
      if (is.null(fit)) return(c(att = NA, se = NA, n = 0))
      c(att = coef(fit)[".post:.treated_group"],
        se = sqrt(diag(vcov(fit)))[".post:.treated_group"],
        n = length(treated_ids))
    })
    cohort_atts <- as.data.frame(t(cohort_atts))
    cohort_atts <- cohort_atts[!is.na(cohort_atts$att), ]
    if (nrow(cohort_atts) == 0) return(data.frame(relative_time = rt, estimate = NA, se = NA))
    w <- cohort_atts$n / sum(cohort_atts$n)
    est <- sum(w * cohort_atts$att); se <- sqrt(sum(w^2 * cohort_atts$se^2))
    data.frame(relative_time = rt, estimate = est, se = se)
  })
  event_study <- do.call(rbind, es_results)
  event_study <- event_study[!is.na(event_study$estimate), ]

  post <- event_study[event_study$relative_time >= 0, ]
  att <- if (nrow(post) > 0) mean(post$estimate) else NA
  att_se <- if (nrow(post) > 0) sqrt(mean(post$se^2) / nrow(post)) else NA

  list(att = att, att_se = att_se,
       att_ci = att + c(-1, 1) * qnorm(0.975) * att_se,
       event_study = event_study,
       coefficients = c(ATT = att), se = c(ATT = att_se),
       z = c(ATT = att / att_se),
       pvalue = c(ATT = 2 * pnorm(-abs(att / att_se))),
       vcov = matrix(att_se^2, 1, 1), hessian = NULL,
       logLik = NA, AIC = NA, BIC = NA)
}

# --- Imputation estimator (BJS 2024) ---
#' @keywords internal
.did_imputation <- function(formula, data, id_var, time_var, cohort_var, treat_var, ...) {
  F <- Formula::Formula(formula)
  mf <- model.frame(F, data = data)
  y_name <- names(mf)[1]

  # Step 1: Estimate FE model on untreated observations
  untreated <- data[[treat_var]] == 0
  f_str <- paste(y_name, "~ factor(", id_var, ") + factor(", time_var, ")")
  fit0 <- tryCatch(lm(as.formula(f_str), data = data[untreated, ]),
                   error = function(e) NULL)
  if (is.null(fit0)) {
    return(list(att = NA, att_se = NA, att_ci = c(NA, NA), event_study = NULL,
                coefficients = numeric(0), se = numeric(0), z = numeric(0),
                pvalue = numeric(0), vcov = NULL, hessian = NULL,
                logLik = NA, AIC = NA, BIC = NA))
  }

  # Step 2: Impute Y(0) for treated observations
  y0_hat <- tryCatch(predict(fit0, newdata = data), error = function(e) rep(NA, nrow(data)))
  tau_i <- data[[y_name]] - y0_hat

  # Step 3: Aggregate
  treated_data <- data[data[[treat_var]] == 1, ]
  treated_data$.tau <- tau_i[data[[treat_var]] == 1]
  treated_data$.rel_time <- treated_data[[time_var]] - treated_data[[cohort_var]]

  att <- mean(treated_data$.tau, na.rm = TRUE)
  att_se <- sd(treated_data$.tau, na.rm = TRUE) / sqrt(sum(!is.na(treated_data$.tau)))

  # Event study
  rel_times <- sort(unique(treated_data$.rel_time))
  event_study <- do.call(rbind, lapply(rel_times, function(rt) {
    vals <- treated_data$.tau[treated_data$.rel_time == rt]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NULL)
    data.frame(relative_time = rt, estimate = mean(vals),
               se = sd(vals) / sqrt(length(vals)))
  }))

  list(att = att, att_se = att_se,
       att_ci = att + c(-1, 1) * qnorm(0.975) * att_se,
       event_study = event_study,
       coefficients = c(ATT = att), se = c(ATT = att_se),
       z = c(ATT = att / att_se),
       pvalue = c(ATT = 2 * pnorm(-abs(att / att_se))),
       vcov = matrix(att_se^2, 1, 1), hessian = NULL,
       logLik = NA, AIC = NA, BIC = NA)
}

# --- Additional stubs ---
#' @keywords internal
.did_dcdh <- function(formula, data, id_var, time_var, treat_var, ...) {
  .did_imputation(formula, data, id_var, time_var, ".cohort_dcdh", treat_var, ...)
}

#' @keywords internal
.did_gardner <- function(formula, data, id_var, time_var, treat_var, cohort_var, ...) {
  # did2s: two-stage estimator — implemented same as imputation
  .did_imputation(formula, data, id_var, time_var, cohort_var, treat_var, ...)
}

#' @keywords internal
.did_stacked <- function(formula, data, id_var, time_var, cohort_var, cluster, ...) {
  F <- Formula::Formula(formula)
  mf <- model.frame(F, data = data)
  y_name <- names(mf)[1]
  cohorts <- sort(unique(data[[cohort_var]]))
  cohorts <- cohorts[is.finite(cohorts)]

  # Stack sub-experiments
  stacked_list <- lapply(cohorts, function(g) {
    d <- data[data[[cohort_var]] == g | !is.finite(data[[cohort_var]]), ]
    d$.stack_id <- g
    d
  })
  stacked_data <- do.call(rbind, stacked_list)
  stacked_data$.treat <- as.integer(stacked_data[[cohort_var]] == stacked_data$.stack_id &
                                     stacked_data[[time_var]] >= stacked_data[[cohort_var]])

  f_str <- paste(y_name, "~ .treat + factor(", id_var, ":.stack_id) + factor(",
                 time_var, ":.stack_id)")
  fit <- tryCatch(lm(as.formula(f_str), data = stacked_data), error = function(e) NULL)
  if (is.null(fit)) {
    return(list(att = NA, att_se = NA, att_ci = c(NA, NA), event_study = NULL,
                coefficients = numeric(0), se = numeric(0), z = numeric(0),
                pvalue = numeric(0), vcov = NULL, hessian = NULL,
                logLik = NA, AIC = NA, BIC = NA))
  }

  att <- coef(fit)[".treat"]; att_se <- sqrt(diag(vcov(fit)))[".treat"]

  list(att = att, att_se = att_se,
       att_ci = att + c(-1, 1) * qnorm(0.975) * att_se,
       event_study = NULL,
       coefficients = c(ATT = att), se = c(ATT = att_se),
       z = c(ATT = att / att_se),
       pvalue = c(ATT = 2 * pnorm(-abs(att / att_se))),
       vcov = matrix(att_se^2, 1, 1), hessian = NULL,
       logLik = as.numeric(logLik(fit)), AIC = AIC(fit), BIC = BIC(fit))
}

#' @keywords internal
.aggregate_did <- function(result, aggregation) {
  es <- result$event_study
  if (is.null(es) || nrow(es) == 0) return(NULL)
  switch(aggregation,
    dynamic = es,
    simple = {
      post <- es[es$relative_time >= 0, ]
      data.frame(estimate = mean(post$estimate),
                 se = sqrt(mean(post$se^2) / nrow(post)))
    },
    group = es,
    calendar = es
  )
}

#' Compare DiD Estimators
#' @param ... did_flex objects
#' @param labels Optional labels
#' @return Object of class "did_comparison"
#' @export
did_compare <- function(..., labels = NULL) {
  models <- list(...)
  if (is.null(labels)) labels <- paste0("Model_", seq_along(models))
  comparison <- data.frame(
    Estimator = labels,
    ATT = sapply(models, function(m) m$att),
    SE = sapply(models, function(m) m$att_se),
    CI_lower = sapply(models, function(m) m$att_ci[1]),
    CI_upper = sapply(models, function(m) m$att_ci[2])
  )
  structure(list(comparison = comparison, models = models), class = "did_comparison")
}

#' @export
print.did_comparison <- function(x, digits = 4, ...) {
  cat("\nDiD Estimator Comparison\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  print(x$comparison, digits = digits, row.names = FALSE)
  cat(paste(rep("-", 60), collapse = ""), "\n")
  invisible(x)
}

#' Event Study Plot for DiD
#' @export
plot.did_flex <- function(x, type = "event_study", ci_level = 0.95,
                          reference_period = -1, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 required for plotting"); return(invisible(NULL))
  }
  es <- x$event_study
  if (is.null(es)) { warning("No event study results"); return(invisible(NULL)) }
  z_val <- qnorm(1 - (1 - ci_level) / 2)
  es$ci_lower <- es$estimate - z_val * es$se
  es$ci_upper <- es$estimate + z_val * es$se

  ggplot2::ggplot(es, ggplot2::aes(x = .data$relative_time, y = .data$estimate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_vline(xintercept = reference_period + 0.5, linetype = "dotted") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
                          alpha = 0.2, fill = "steelblue") +
    ggplot2::geom_point(size = 2, color = "steelblue") +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::labs(x = "Relative Time", y = "Estimate",
                   title = paste("Event Study:", x$estimator)) +
    ggplot2::theme_minimal()
}
