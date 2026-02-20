# =============================================================================
# ecoflex: LaTeX Output
# =============================================================================

#' Export Model Results to LaTeX
#'
#' Generates publication-ready LaTeX tables. For two-stage models (hurdle,
#' heckman, tobit II+, switching), results are displayed side-by-side.
#'
#' @param object An ecoflex model object
#' @param file Optional file path to write the LaTeX code. If NULL, prints to console.
#' @param caption Table caption
#' @param label LaTeX label for referencing
#' @param digits Number of decimal places. Default: 4
#' @param stars Show significance stars. Default: TRUE
#' @param se_type "se" (standard errors), "t" (t-statistics), "p" (p-values)
#' @param notes Footnote text (character vector)
#' @param adjustbox Wrap table in adjustbox. Default: TRUE
#' @param adjustbox_width Width for adjustbox. Default: "\\textwidth"
#' @param font_size LaTeX font size command. Default: NULL (no change)
#' @param booktabs Use booktabs package rules. Default: TRUE
#' @param preamble Include \\usepackage directives. Default: FALSE
#' @param ... Additional arguments
#' @return Invisible character vector of LaTeX lines.
#' @export
#'
#' @examples
#' \donttest{
#' m <- tobit_flex(mpg ~ hp + wt, data = mtcars, tobit_type = 1, left = 15)
#' to_latex(m)
#' }
to_latex <- function(object, ...) UseMethod("to_latex")

#' @rdname to_latex
#' @export
to_latex.ecoflex <- function(object, file = NULL, caption = NULL, label = NULL,
                              digits = 4, stars = TRUE, se_type = "se",
                              notes = NULL, adjustbox = TRUE,
                              adjustbox_width = "\\textwidth",
                              font_size = NULL,
                              booktabs = TRUE, preamble = FALSE, ...) {

  model_class <- class(object)[1]

  # Dispatch to specialised formatter for two-stage models
  if (model_class %in% c("hurdle_flex", "heckman_flex", "switching_flex") ||
      (model_class == "tobit_flex" && object$tobit_type >= 2)) {
    lines <- .latex_two_stage(object, caption, label, digits, stars,
                               se_type, notes, adjustbox, adjustbox_width,
                               font_size, booktabs, preamble)
  } else {
    lines <- .latex_single(object, caption, label, digits, stars,
                            se_type, notes, adjustbox, adjustbox_width,
                            font_size, booktabs, preamble)
  }

  out <- paste(lines, collapse = "\n")
  if (!is.null(file)) {
    writeLines(out, file)
    message("LaTeX table written to: ", file)
  } else {
    cat(out, "\n")
  }
  invisible(lines)
}

# Aliases for each model class
#' @export
to_latex.hurdle_flex <- to_latex.ecoflex
#' @export
to_latex.heckman_flex <- to_latex.ecoflex
#' @export
to_latex.rdd_flex <- to_latex.ecoflex
#' @export
to_latex.ivgmm_flex <- to_latex.ecoflex
#' @export
to_latex.did_flex <- to_latex.ecoflex
#' @export
to_latex.tobit_flex <- to_latex.ecoflex
#' @export
to_latex.panel_nl <- to_latex.ecoflex
#' @export
to_latex.switching_flex <- to_latex.ecoflex


# ============================================================
# Internal: single-column table
# ============================================================
#' @keywords internal
.latex_single <- function(object, caption, label, digits, stars, se_type,
                           notes, adjustbox, adjustbox_width, font_size,
                           booktabs, preamble) {

  coefs <- object$coefficients
  se    <- object$se
  z     <- if (!is.null(object$z)) object$z else coefs / se
  pval  <- if (!is.null(object$pvalue)) object$pvalue else 2 * pnorm(-abs(z))

  nms   <- .latex_escape(names(coefs))
  if (is.null(nms)) nms <- paste0("$\\beta_{", seq_along(coefs), "}$")

  lines <- character()
  if (preamble) lines <- c(lines, .latex_preamble(adjustbox, booktabs))

  lines <- c(lines, .latex_table_open(caption, label, adjustbox,
                                       adjustbox_width, font_size, booktabs,
                                       ncol = 2, headers = c("", object$model_name %||% "Estimate")))

  for (i in seq_along(coefs)) {
    val <- .format_coef(coefs[i], pval[i], digits, stars)
    below <- .format_below(se[i], z[i], pval[i], se_type, digits)
    lines <- c(lines,
      paste0("  ", nms[i], " & ", val, " \\\\"),
      paste0("  & ", below, " \\\\")
    )
  }

  # Fit statistics
  lines <- c(lines, .latex_rule(booktabs))
  lines <- .add_fit_stats(lines, object, 1, booktabs)

  lines <- c(lines, .latex_table_close(notes, adjustbox, stars, booktabs))
  lines
}


# ============================================================
# Internal: side-by-side table for two-stage models
# ============================================================
#' @keywords internal
.latex_two_stage <- function(object, caption, label, digits, stars, se_type,
                              notes, adjustbox, adjustbox_width, font_size,
                              booktabs, preamble) {

  parts <- .split_stages(object)
  n_stages <- length(parts)

  all_names <- unique(unlist(lapply(parts, function(p) names(p$coefs))))
  all_names_esc <- .latex_escape(all_names)

  stage_labels <- names(parts)
  headers <- c("", stage_labels)

  lines <- character()
  if (preamble) lines <- c(lines, .latex_preamble(adjustbox, booktabs))

  lines <- c(lines, .latex_table_open(caption, label, adjustbox,
                                       adjustbox_width, font_size, booktabs,
                                       ncol = n_stages + 1, headers = headers))

  for (j in seq_along(all_names)) {
    nm <- all_names[j]
    nm_esc <- all_names_esc[j]
    vals <- sapply(parts, function(p) {
      idx <- match(nm, names(p$coefs))
      if (is.na(idx)) return("")
      .format_coef(p$coefs[idx], p$pval[idx], digits, stars)
    })
    belows <- sapply(parts, function(p) {
      idx <- match(nm, names(p$coefs))
      if (is.na(idx)) return("")
      .format_below(p$se[idx], p$z[idx], p$pval[idx], se_type, digits)
    })
    lines <- c(lines,
      paste0("  ", nm_esc, " & ", paste(vals, collapse = " & "), " \\\\"),
      paste0("  & ", paste(belows, collapse = " & "), " \\\\")
    )
  }

  lines <- c(lines, .latex_rule(booktabs))
  lines <- .add_fit_stats(lines, object, n_stages, booktabs)
  lines <- c(lines, .latex_table_close(notes, adjustbox, stars, booktabs))
  lines
}


# ============================================================
# Stage splitter — knows how to decompose each model type
# ============================================================
#' @keywords internal
.split_stages <- function(object) {
  cls <- class(object)[1]

  if (cls == "hurdle_flex") {
    coefs <- object$coefficients; se <- object$se
    z <- if (!is.null(object$z)) object$z else coefs / se
    pval <- if (!is.null(object$pvalue)) object$pvalue else 2 * pnorm(-abs(z))

    # Split by name pattern: binary_ vs count_
    is_binary <- grepl("^binary_", names(coefs))
    is_count  <- grepl("^count_", names(coefs))
    if (any(is_binary) || any(is_count)) {
      b_names <- sub("^binary_", "", names(coefs)[is_binary])
      c_names <- sub("^count_", "", names(coefs)[is_count])
      other <- !is_binary & !is_count
      return(list(
        `Zero/Binary` = list(coefs = setNames(coefs[is_binary], b_names),
                             se = se[is_binary], z = z[is_binary], pval = pval[is_binary]),
        `Count/Hurdle` = list(coefs = setNames(coefs[is_count], c_names),
                              se = se[is_count], z = z[is_count], pval = pval[is_count])
      ))
    }
    # Fallback: split in half
    n <- length(coefs); mid <- ceiling(n / 2)
    return(list(
      `Zero/Binary` = list(coefs = coefs[1:mid], se = se[1:mid],
                           z = z[1:mid], pval = pval[1:mid]),
      `Count/Hurdle` = list(coefs = coefs[(mid+1):n], se = se[(mid+1):n],
                            z = z[(mid+1):n], pval = pval[(mid+1):n])
    ))
  }

  if (cls == "heckman_flex") {
    coefs <- object$coefficients; se <- object$se
    z <- if (!is.null(object$z)) object$z else coefs / se
    pval <- if (!is.null(object$pvalue)) object$pvalue else 2 * pnorm(-abs(z))
    is_sel <- grepl("^selection_", names(coefs))
    is_out <- grepl("^outcome_", names(coefs))
    other  <- !is_sel & !is_out
    s_names <- sub("^selection_", "", names(coefs)[is_sel])
    o_names <- sub("^outcome_", "", names(coefs)[is_out])
    out <- list(
      Selection = list(coefs = setNames(coefs[is_sel], s_names),
                       se = se[is_sel], z = z[is_sel], pval = pval[is_sel]),
      Outcome   = list(coefs = setNames(coefs[is_out], o_names),
                        se = se[is_out], z = z[is_out], pval = pval[is_out])
    )
    if (any(other)) {
      out[["Ancillary"]] <- list(coefs = coefs[other], se = se[other],
                                  z = z[other], pval = pval[other])
    }
    return(out)
  }

  if (cls == "switching_flex") {
    coefs <- object$coefficients; se <- object$se
    z <- if (!is.null(object$z)) object$z else coefs / se
    pval <- if (!is.null(object$pvalue)) object$pvalue else 2 * pnorm(-abs(z))
    # Split by regime prefix
    regime_ids <- unique(sub("_.*", "", names(coefs)[grepl("^regime", names(coefs))]))
    if (length(regime_ids) >= 2) {
      out <- lapply(regime_ids, function(rid) {
        idx <- grepl(paste0("^", rid, "_"), names(coefs))
        clean_nm <- sub(paste0("^", rid, "_"), "", names(coefs)[idx])
        list(coefs = setNames(coefs[idx], clean_nm),
             se = se[idx], z = z[idx], pval = pval[idx])
      })
      names(out) <- gsub("regime", "Regime ", regime_ids)
      sel_idx <- grepl("^selection_|^atanh_|^log_sigma|^trans_", names(coefs))
      if (any(sel_idx)) {
        out[["Selection/Ancillary"]] <- list(coefs = coefs[sel_idx], se = se[sel_idx],
                                              z = z[sel_idx], pval = pval[sel_idx])
      }
      return(out)
    }
  }

  if (cls == "tobit_flex" && object$tobit_type >= 2) {
    coefs <- object$coefficients; se <- object$se
    z <- if (!is.null(object$z)) object$z else coefs / se
    pval <- if (!is.null(object$pvalue)) object$pvalue else 2 * pnorm(-abs(z))
    n <- length(coefs); mid <- ceiling(n / 2)
    return(list(
      Selection = list(coefs = coefs[1:mid], se = se[1:mid],
                       z = z[1:mid], pval = pval[1:mid]),
      Outcome   = list(coefs = coefs[(mid+1):n], se = se[(mid+1):n],
                        z = z[(mid+1):n], pval = pval[(mid+1):n])
    ))
  }

  # Generic fallback: single stage
  coefs <- object$coefficients; se <- object$se
  z <- if (!is.null(object$z)) object$z else coefs / se
  pval <- if (!is.null(object$pvalue)) object$pvalue else 2 * pnorm(-abs(z))
  list(Estimate = list(coefs = coefs, se = se, z = z, pval = pval))
}


# ============================================================
# Formatting helpers
# ============================================================

#' @keywords internal
.format_coef <- function(val, pval, digits, stars) {
  if (is.na(val)) return("")
  s <- formatC(val, digits = digits, format = "f")
  if (stars && !is.na(pval)) {
    star <- if (pval < 0.001) "^{***}" else if (pval < 0.01) "^{**}" else
            if (pval < 0.05) "^{*}" else ""
    s <- paste0("$", s, star, "$")
  }
  s
}

#' @keywords internal
.format_below <- function(se, z, pval, se_type, digits) {
  val <- switch(se_type,
    se = se, t = z, p = pval
  )
  if (is.na(val)) return("")
  paste0("(", formatC(val, digits = digits, format = "f"), ")")
}

#' @keywords internal
.latex_escape <- function(x) {
  if (is.null(x)) return(NULL)
  x <- gsub("_", "\\_", x, fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x <- gsub("&", "\\&", x, fixed = TRUE)
  x <- gsub("#", "\\#", x, fixed = TRUE)
  x
}

#' @keywords internal
.latex_preamble <- function(adjustbox, booktabs) {
  pkgs <- c("\\usepackage{booktabs}")
  if (adjustbox) pkgs <- c(pkgs, "\\usepackage{adjustbox}")
  pkgs <- c(pkgs, "\\usepackage{threeparttable}")
  pkgs
}

#' @keywords internal
.latex_table_open <- function(caption, label, adjustbox, adjustbox_width,
                               font_size, booktabs, ncol, headers) {
  lines <- character()
  lines <- c(lines, "\\begin{table}[htbp]", "\\centering")
  if (!is.null(caption)) lines <- c(lines, paste0("\\caption{", caption, "}"))
  if (!is.null(label)) lines <- c(lines, paste0("\\label{", label, "}"))
  if (!is.null(font_size)) lines <- c(lines, paste0("\\", font_size))
  lines <- c(lines, "\\begin{threeparttable}")
  if (adjustbox) {
    lines <- c(lines, paste0("\\begin{adjustbox}{max width=", adjustbox_width, "}"))
  }
  col_spec <- paste0("l", paste(rep("c", ncol - 1), collapse = ""))
  lines <- c(lines, paste0("\\begin{tabular}{", col_spec, "}"))
  if (booktabs) lines <- c(lines, "\\toprule") else lines <- c(lines, "\\hline")
  # Header row
  lines <- c(lines, paste0("  ", paste(headers, collapse = " & "), " \\\\"))
  if (booktabs) lines <- c(lines, "\\midrule") else lines <- c(lines, "\\hline")
  lines
}

#' @keywords internal
.latex_rule <- function(booktabs) {
  if (booktabs) "\\midrule" else "\\hline"
}

#' @keywords internal
.add_fit_stats <- function(lines, object, n_stages, booktabs) {
  mc <- function(val, ncol) {
    if (ncol > 1) paste0("\\multicolumn{", ncol, "}{c}{", val, "}")
    else val
  }
  nc <- n_stages

  if (!is.null(object$n) && !is.na(object$n))
    lines <- c(lines, paste0("  Observations & ", mc(object$n, nc), " \\\\"))
  if (!is.null(object$logLik) && !is.na(object$logLik))
    lines <- c(lines, paste0("  Log-Likelihood & ", mc(formatC(object$logLik, digits = 2, format = "f"), nc), " \\\\"))
  if (!is.null(object$AIC) && !is.na(object$AIC))
    lines <- c(lines, paste0("  AIC & ", mc(formatC(object$AIC, digits = 2, format = "f"), nc), " \\\\"))
  if (!is.null(object$BIC) && !is.na(object$BIC))
    lines <- c(lines, paste0("  BIC & ", mc(formatC(object$BIC, digits = 2, format = "f"), nc), " \\\\"))
  lines
}

#' @keywords internal
.latex_table_close <- function(notes, adjustbox, stars, booktabs) {
  lines <- character()
  if (booktabs) lines <- c(lines, "\\bottomrule") else lines <- c(lines, "\\hline")
  lines <- c(lines, "\\end{tabular}")
  if (adjustbox) lines <- c(lines, "\\end{adjustbox}")

  # Table notes
  tnotes <- character()
  if (stars) tnotes <- c(tnotes, "$^{***}$p$<$0.001; $^{**}$p$<$0.01; $^{*}$p$<$0.05")
  if (!is.null(notes)) tnotes <- c(tnotes, notes)
  if (length(tnotes) > 0) {
    lines <- c(lines, "\\begin{tablenotes}[flushleft]", "\\small")
    for (n in tnotes) lines <- c(lines, paste0("\\item ", n))
    lines <- c(lines, "\\end{tablenotes}")
  }

  lines <- c(lines, "\\end{threeparttable}", "\\end{table}")
  lines
}

# NULL-coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a
