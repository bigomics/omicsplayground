##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# MOFA AI Report — data extraction and template rendering
# =============================================================================
# All deterministic data functions consumed by the AI-report Shiny module
# (mofa_ai_text_server.R). No Shiny dependencies. No LLM calls.
#
# Layout (top to bottom):
#   - leaf renderers, ONE per {{placeholder}} in mofa_report_data.md:
#       mofa_data_overview, mofa_data_contrasts, mofa_data_variance,
#       mofa_data_factors_summary, mofa_data_factor_correlations, mofa_data_factor_detail
#   - per-factor block renderer:    .mofa_render_factor_block (via FACTOR_TEMPLATE)
#   - orchestrators:                 mofa_build_report_tables,
#                                    mofa_build_summary_params,
#                                    mofa_build_methods
#
# Verbal-label rendering uses omicsai::omicsai_verbalize_{r,q} — keep them
# in the prompt data block; do not re-introduce raw numbers in rendered
# text except where explicitly documented.

MOFA_PROMPTS_DIR <- file.path(OPG, "components/board.mofa/prompts/mofa")


# =============================================================================
# Verbalisers — keep together; promote to omicsai during the compaction pass
# =============================================================================
# Pattern matches omicsai::omicsai_verbalize_{r,q,logfc}: pure value → label
# functions with explicit threshold tables. Authoritative thresholds also
# appear in mofa_interpretation.md (any change there MUST be mirrored here
# and vice versa — the LLM is told what each label means).

#' Verbalise a MOFA feature weight (no direction; sign belongs with the trait).
#' Bins: dominant (|w| ≥ 1.0) / strong (≥ 0.5) / modest (≥ 0.25) / minimal.
mofa_verbalize_weight <- function(w, breaks = c(0.25, 0.5, 1.0),
                                  na_label = "not tested") {
  out <- rep(na_label, length(w))
  ok  <- !is.na(w)
  a   <- abs(w[ok])
  lab <- ifelse(a >= breaks[3], "dominant",
         ifelse(a >= breaks[2], "strong",
         ifelse(a >= breaks[1], "modest", "minimal")))
  out[ok] <- lab
  out
}

#' Verbalise a GSEA NES — direction inherits from sign; magnitude bucketed.
#' Bins: up/down-strong (|NES| ≥ 2.5) / -moderate (≥ 1.5) / -weak (≥ 0.5) /
#'       -nominal (< 0.5).
mofa_verbalize_nes <- function(nes, breaks = c(0.5, 1.5, 2.5),
                               na_label = "not tested") {
  out <- rep(na_label, length(nes))
  ok  <- !is.na(nes)
  a   <- abs(nes[ok])
  dir <- ifelse(nes[ok] >= 0, "up-", "down-")
  mag <- ifelse(a >= breaks[3], "strong",
         ifelse(a >= breaks[2], "moderate",
         ifelse(a >= breaks[1], "weak", "nominal")))
  out[ok] <- paste0(dir, mag)
  out
}

#' Verbalise a per-view variance-explained value. Accepts either a fraction
#' (0–1) or a percentage; auto-detects.
#' Bins: major (≥ 30%) / moderate (≥ 15%) / minor (≥ 5%) / negligible.
mofa_verbalize_variance <- function(v, breaks = c(5, 15, 30),
                                    na_label = "not tested") {
  out <- rep(na_label, length(v))
  ok  <- !is.na(v)
  pct <- v[ok]
  if (length(pct) > 0 && max(pct, na.rm = TRUE) <= 1) pct <- 100 * pct
  out[ok] <- ifelse(pct >= breaks[3], "major",
              ifelse(pct >= breaks[2], "moderate",
              ifelse(pct >= breaks[1], "minor", "negligible")))
  out
}

#' Verbalise a cross-view network centrality score in [0, 1].
#' Bins: hub (≥ 0.8) / central (≥ 0.6) / intermediate (≥ 0.3) / peripheral.
mofa_verbalize_centrality <- function(c_, breaks = c(0.3, 0.6, 0.8),
                                      na_label = "not tested") {
  out <- rep(na_label, length(c_))
  ok  <- !is.na(c_)
  out[ok] <- ifelse(c_[ok] >= breaks[3], "hub",
              ifelse(c_[ok] >= breaks[2], "central",
              ifelse(c_[ok] >= breaks[1], "intermediate", "peripheral")))
  out
}


# -----------------------------------------------------------------------------
# Factor data computation (gather everything needed per factor, once)
# -----------------------------------------------------------------------------

#' Compute per-factor data dicts for every requested factor.
#'
#' Returns a named list, one element per factor, each carrying the raw
#' values used by both Summary mode and Report mode.
.mofa_compute_factor_data <- function(mofa, pgx, factors, ntop = 10L) {
  out <- list()
  for (f in factors) {
    ctx <- mofa_ai_extract_factor_data(mofa, pgx, f, ntop = ntop)
    if (!is.null(ctx)) out[[f]] <- ctx
  }
  out
}

#' Tier per factor (strong/moderate/weak) from per-factor metrics.
.mofa_factor_tier <- function(metrics) {
  score <-
    0.45 * min(metrics$n_sig_pathways / 10, 1) +
    0.30 * min(metrics$max_abs_nes / 2.5, 1) +
    0.25 * min(metrics$max_abs_weight / 1.5, 1)
  if (score >= 0.70) "strong"
  else if (score >= 0.45) "moderate"
  else "weak"
}

#' Order factors by signal score, descending. Returns a character vector.
.mofa_factor_order <- function(factor_data) {
  if (length(factor_data) == 0) return(character(0))
  scores <- vapply(factor_data, function(fd) {
    0.45 * min(fd$metrics$n_sig_pathways / 10, 1) +
    0.30 * min(fd$metrics$max_abs_nes / 2.5, 1) +
    0.25 * min(fd$metrics$max_abs_weight / 1.5, 1)
  }, numeric(1))
  names(sort(scores, decreasing = TRUE))
}


# -----------------------------------------------------------------------------
# Leaf renderers — ONE per {{placeholder}} in mofa_report_data.md
# -----------------------------------------------------------------------------

#' {{experiment}}, {{organism}}, {{n_samples}}, {{n_views}}, {{view_names}},
#' {{n_factors_total}}, {{n_factors_used}} — the Overview metadata block.
mofa_data_overview <- function(mofa, pgx, n_factors_used) {
  experiment <- multiomics_ai_experiment_label(pgx, mofa$experiment)
  organism   <- pgx$organism %||% "unknown"
  n_samples  <- if (!is.null(mofa$F)) nrow(mofa$F) else NA_integer_

  view_names <- character(0)
  if (!is.null(mofa$ww) && is.list(mofa$ww)) {
    view_names <- names(mofa$ww)
  } else if (!is.null(mofa$views)) {
    view_names <- as.character(mofa$views)
  }
  n_views <- length(view_names)
  view_names_str <- if (n_views > 0) paste(view_names, collapse = ", ") else "unknown"

  n_factors_total <- length(mofa_ai_factor_choices(mofa))

  list(
    experiment       = experiment,
    organism         = organism,
    n_samples        = as.character(n_samples %||% "unknown"),
    n_views          = as.character(n_views),
    view_names       = view_names_str,
    n_factors_total  = as.character(n_factors_total),
    n_factors_used   = as.character(n_factors_used)
  )
}

#' {{contrasts_block}} — contrast names from PGX, one per line.
mofa_data_contrasts <- function(pgx) {
  contrasts <- tryCatch(colnames(playbase::pgx.getMetaMatrix(pgx)$fc),
                        error = function(e) character(0))
  if (length(contrasts) == 0) return("(no contrasts available)")
  paste(sprintf("- %s", contrasts), collapse = "\n")
}

#' {{variance_block}} — per-factor variance explained across views,
#' verbalised. Raw % is kept only as a parenthetical anchor on the lead
#' view (the view with the highest variance) per factor — mirrors the
#' WGCNA "raw r on the lead trait" convention.
mofa_data_variance <- function(mofa, factor_order) {
  if (is.null(mofa$variance) ||
      (!is.matrix(mofa$variance) && !is.data.frame(mofa$variance))) {
    return("(variance-explained matrix not available)")
  }
  v <- as.matrix(mofa$variance)
  keep <- intersect(rownames(v), factor_order)
  if (length(keep) == 0) {
    return("(variance-explained values not available for selected factors)")
  }
  v <- v[keep, , drop = FALSE]

  rows <- list()
  for (f in rownames(v)) {
    vec <- v[f, ]
    lead <- which.max(vec)
    cells <- mofa_verbalize_variance(vec)
    pct <- if (max(vec, na.rm = TRUE) <= 1) 100 * vec else vec
    cells[lead] <- sprintf("%s (%s%%)", cells[lead],
                           omicsai::omicsai_format_num(pct[lead], 0))
    rows[[f]] <- c(Factor = f, cells)
  }
  df <- as.data.frame(do.call(rbind, rows),
                      check.names = FALSE, stringsAsFactors = FALSE)
  paste(omicsai::omicsai_format_mdtable(df), collapse = "\n")
}

#' {{factors_summary_table}}, {{lead_factor}} — cross-factor summary
#' table + the lead factor identifier (template owns the prose framing).
mofa_data_factors_summary <- function(factor_data, factor_order, lead_factor) {
  if (length(factor_order) == 0) {
    return(list(table = "(no factors to summarise)", footnote = ""))
  }
  rows <- lapply(factor_order, function(f) {
    fd <- factor_data[[f]]
    tier <- .mofa_factor_tier(fd$metrics)
    top_pathway_theme <- if (!is.null(fd$top_pathways) && nrow(fd$top_pathways) > 0) {
      sub(".*:", "", fd$top_pathways$pathway[1])
    } else "—"
    data.frame(
      Factor             = f,
      Tier               = tier,
      `Top trait`        = fd$traits %||% "—",
      `Top pathway theme`= top_pathway_theme,
      `Sig. pathways`    = as.character(fd$metrics$n_sig_pathways),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  body <- paste(omicsai::omicsai_format_mdtable(df), collapse = "\n")

  ## Prose framing lives in mofa_report_data.md; we just return the
  ## lead-factor identifier (or "—" when none).
  list(table = body,
       lead_factor = if (!is.na(lead_factor)) lead_factor else "—")
}

#' {{factor_correlations}} — factor-factor score correlations |r| >= 0.30.
mofa_data_factor_correlations <- function(mofa, factor_order, thresh = 0.30) {
  if (is.null(mofa$F) || length(factor_order) < 2) {
    return("(insufficient factors for cross-factor correlation)")
  }
  keep <- intersect(factor_order, colnames(mofa$F))
  if (length(keep) < 2) return("(no factor scores available for correlation)")
  cors <- suppressWarnings(stats::cor(mofa$F[, keep, drop = FALSE],
                                      use = "pairwise.complete.obs"))
  pairs <- character(0)
  for (i in seq_len(nrow(cors) - 1L)) {
    for (j in (i + 1L):ncol(cors)) {
      r <- cors[i, j]
      if (!is.na(r) && abs(r) >= thresh) {
        pairs <- c(pairs, sprintf("- %s ↔ %s: %s (r = %s)",
          rownames(cors)[i], colnames(cors)[j],
          omicsai::omicsai_verbalize_r(r),
          omicsai::omicsai_format_num(r, 2)))
      }
    }
  }
  if (length(pairs) == 0) return("(no factor pairs with |r| ≥ 0.30)")
  paste(pairs, collapse = "\n")
}

#' {{factor_detail}} — concatenated per-factor blocks.
#'
#' `variance_mat` is the optional per-factor × per-view variance-explained
#' matrix (`mofa$variance`); when supplied, each factor block carries a
#' verbalised one-line summary under `{{variance_line}}`.
mofa_data_factor_detail <- function(factor_data, factor_order,
                               variance_mat = NULL,
                               ntop = 10L) {
  if (length(factor_order) == 0) return("")
  blocks <- vapply(factor_order, function(f) {
    .mofa_render_factor_block(factor_data[[f]], f,
                         variance_mat = variance_mat, ntop = ntop)
  }, character(1))
  paste(blocks, collapse = "\n\n")
}


# -----------------------------------------------------------------------------
# Per-factor block renderer (shared by Summary + Report)
# -----------------------------------------------------------------------------

## FACTOR_TEMPLATE — verbatim spec, loaded from disk.
## Edit prompts/mofa/mofa_factor_data.md directly to change wording;
## the R renderer only fills values.
FACTOR_TEMPLATE <- omicsai::omicsai_load_template(
  file.path(MOFA_PROMPTS_DIR, "mofa_factor_data.md")
)

#' Render one factor's data dict into a FACTOR_TEMPLATE block.
#'
#' All numeric quantities (weight, NES, padj, centrality, variance) are
#' verbalised via the verbalisers at the top of this file. The single
#' exception is the per-view variance percentage on the lead view, which
#' is kept as a parenthetical anchor (mirrors the WGCNA "raw r on lead
#' trait" convention).
.mofa_render_factor_block <- function(fd, factor_name,
                                 variance_mat = NULL, ntop = 10L) {
  tier <- .mofa_factor_tier(fd$metrics)

  view_mix <- if (length(fd$datatype_counts) > 0) {
    paste(sprintf("%s=%s", names(fd$datatype_counts),
                  as.integer(fd$datatype_counts)),
          collapse = ", ")
  } else "unknown"

  ## --- Per-factor variance line (verbalised, lead view gets raw %) ---
  variance_line <- "—"
  if (!is.null(variance_mat) && factor_name %in% rownames(variance_mat)) {
    vec <- variance_mat[factor_name, ]
    if (length(vec) > 0 && any(!is.na(vec))) {
      labels <- mofa_verbalize_variance(vec)
      pct <- if (max(vec, na.rm = TRUE) <= 1) 100 * vec else vec
      lead <- which.max(vec)
      cells <- sprintf("%s %s", names(vec), labels)
      cells[lead] <- sprintf("%s %s (%s%%)", names(vec)[lead], labels[lead],
                             omicsai::omicsai_format_num(pct[lead], 0))
      variance_line <- paste(cells, collapse = ", ")
    }
  }

  trait_summary <- fd$traits %||% "—"

  pathways_str <- if (!is.null(fd$top_pathways) && nrow(fd$top_pathways) > 0) {
    p <- head(fd$top_pathways, ntop)
    df <- data.frame(
      Rank                   = seq_len(nrow(p)),
      Pathway                = sub(".*:", "", p$pathway),
      `Direction & strength` = mofa_verbalize_nes(p$NES),
      Significance           = omicsai::omicsai_verbalize_q(p$padj),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    paste(omicsai::omicsai_format_mdtable(df), collapse = "\n")
  } else "—"

  features_str <- if (!is.null(fd$top_features) && nrow(fd$top_features) > 0) {
    f <- head(fd$top_features, ntop)
    df <- data.frame(
      Symbol         = paste0("*", f$symbol, "*"),
      Contribution   = mofa_verbalize_weight(f$weight),
      `Network role` = mofa_verbalize_centrality(f$centrality),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    paste(omicsai::omicsai_format_mdtable(df), collapse = "\n")
  } else "—"

  cross_view_str <- if (!is.null(fd$cross_view_features) &&
                       length(fd$cross_view_features) > 0) {
    paste(sprintf("*%s*", fd$cross_view_features), collapse = ", ")
  } else "—"

  omicsai::omicsai_substitute_template(FACTOR_TEMPLATE, list(
    factor               = factor_name,
    n_features           = as.character(fd$metrics$n_features),
    tier                 = tier,
    view_mix             = view_mix,
    variance_line        = variance_line,
    trait_summary        = trait_summary,
    n_sig_pathways       = as.character(fd$metrics$n_sig_pathways),
    n_total_pathways     = as.character(fd$metrics$n_total_pathways %||%
                                         fd$metrics$n_sig_pathways),
    pathway_themes_table = pathways_str,
    n_top                = as.character(min(ntop, fd$metrics$n_features)),
    top_features_table   = features_str,
    cross_view_features  = cross_view_str
  ))
}


# -----------------------------------------------------------------------------
# Orchestrator
# -----------------------------------------------------------------------------

#' Build structured report tables from MOFA results.
#'
#' Renders the canonical data block: a single markdown document substituted
#' into `prompts/mofa/mofa_report_data.md`, plus a structured `data` list
#' for any downstream callers that want the raw values.
#'
#' @return list(text = character, data = list, ranking = data.frame)
mofa_build_report_tables <- function(mofa, pgx,
                                     n_factors = 8L, ntop = 10L,
                                     include_variance = TRUE,
                                     include_contrasts = TRUE) {
  all_factors <- mofa_ai_factor_choices(mofa)
  if (length(all_factors) == 0) {
    return(list(text = "(no MOFA factors available)",
                data = list(), ranking = data.frame()))
  }

  factor_data <- .mofa_compute_factor_data(mofa, pgx, all_factors, ntop = ntop)
  if (length(factor_data) == 0) {
    return(list(text = "(no reportable MOFA factor contexts)",
                data = list(), ranking = data.frame()))
  }

  ordered <- .mofa_factor_order(factor_data)
  keep <- head(ordered, as.integer(n_factors))
  factor_data <- factor_data[keep]
  lead_factor <- if (length(keep) > 0) keep[1] else NA_character_

  overview <- mofa_data_overview(mofa, pgx, n_factors_used = length(keep))

  contrasts_block <- if (include_contrasts) mofa_data_contrasts(pgx)
                     else "(contrasts omitted)"

  variance_block <- if (include_variance) mofa_data_variance(mofa, keep)
                    else "(variance-explained omitted)"

  variance_mat <- if (!is.null(mofa$variance) &&
                      (is.matrix(mofa$variance) ||
                       is.data.frame(mofa$variance))) {
    as.matrix(mofa$variance)
  } else NULL

  factors_summary <- mofa_data_factors_summary(factor_data, keep, lead_factor)
  factor_corrs    <- mofa_data_factor_correlations(mofa, keep)
  per_factor      <- mofa_data_factor_detail(factor_data, keep,
                                        variance_mat = variance_mat,
                                        ntop = ntop)

  tmpl <- omicsai::omicsai_load_template(
    file.path(MOFA_PROMPTS_DIR, "mofa_report_data.md")
  )
  text <- omicsai::omicsai_substitute_template(tmpl, c(
    overview,
    list(
      contrasts_block          = contrasts_block,
      variance_block           = variance_block,
      factors_summary_table    = factors_summary$table,
      lead_factor              = factors_summary$lead_factor,
      factor_correlations      = factor_corrs,
      factor_detail            = per_factor
    )
  ))

  ranking_rows <- lapply(keep, function(f) {
    fd <- factor_data[[f]]
    data.frame(
      factor         = f,
      tier           = .mofa_factor_tier(fd$metrics),
      n_sig_pathways = fd$metrics$n_sig_pathways,
      max_abs_nes    = fd$metrics$max_abs_nes,
      max_abs_weight = fd$metrics$max_abs_weight,
      stringsAsFactors = FALSE
    )
  })
  ranking <- if (length(ranking_rows) > 0) do.call(rbind, ranking_rows)
             else data.frame()

  list(text = text, data = factor_data, ranking = ranking)
}


# -----------------------------------------------------------------------------
# Per-factor Summary mode (single-factor, used by the Summary tab)
# -----------------------------------------------------------------------------

#' Build prompt parameters for a single MOFA factor summary.
#'
#' Produces the same per-factor block shape as Report mode by routing
#' through `.mofa_render_factor_block()`. The downstream `{{module_detail}}`
#' placeholder in `wgcna_summary.md`-equivalent template receives the
#' rendered FACTOR_TEMPLATE block.
#'
#' @return Named list: experiment, factor, factor_detail (+ a few extra
#'   convenience fields kept for backward compatibility with the existing
#'   summary template wording).
mofa_build_summary_params <- function(mofa, factor_name, pgx, ntop = 12L) {
  fd <- mofa_ai_extract_factor_data(mofa, pgx, factor_name, ntop = ntop)
  if (is.null(fd)) return(NULL)

  variance_mat <- if (!is.null(mofa$variance) &&
                      (is.matrix(mofa$variance) ||
                       is.data.frame(mofa$variance))) {
    as.matrix(mofa$variance)
  } else NULL

  factor_detail <- .mofa_render_factor_block(fd, factor_name,
                                        variance_mat = variance_mat,
                                        ntop = ntop)

  ## Summary mode renders mofa_summary.md (outer) which embeds the
  ## inner FACTOR_TEMPLATE via {{factor_detail}}. All prose lives in
  ## the templates; this function returns raw values only.
  list(
    experiment    = fd$experiment,
    factor        = factor_name,
    factor_detail = factor_detail
  )
}


# -----------------------------------------------------------------------------
# Methods section (deterministic appendix)
# -----------------------------------------------------------------------------

#' Build deterministic methods section for MOFA report.
mofa_build_methods <- function(mofa, pgx) {
  template <- omicsai::omicsai_load_template(
    file.path(MOFA_PROMPTS_DIR, "mofa_methods.md")
  )

  report_date <- format(Sys.Date(), "%Y-%m-%d")
  params <- list(
    experiment = multiomics_ai_experiment_label(pgx, mofa$experiment),
    n_factors  = length(mofa_ai_factor_choices(mofa)),
    n_samples  = if (!is.null(mofa$F)) nrow(mofa$F) else NA_integer_,
    n_features = if (!is.null(mofa$W)) nrow(mofa$W) else NA_integer_,
    date       = report_date
  )

  omicsai::collapse_lines(
    omicsai::omicsai_substitute_template(template, params),
    sprintf("_This report was generated with OmicsPlayground (BigOmics, %s)._",
            report_date),
    "_Note: AI-generated interpretation may contain inaccuracies and must be independently verified._",
    sep = "\n\n"
  )
}
