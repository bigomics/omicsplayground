## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

# =============================================================================
# WGCNA AI Report — data extraction and template rendering
# =============================================================================
# Single home for all deterministic data functions consumed by the AI-report
# Shiny module (ai_text_server.R). No Shiny dependencies. No LLM calls.
#
# Layout (top to bottom):
#   - extraction primitives (resolve_*, display_*, extract_module_data,
#     .compute_module_data, .compute_families_text)
#   - section builders, one per template slot:
#       data_overview, data_contrast, data_modsummary,
#       data_eigen_cor, data_module_detail
#   - orchestrator:        wgcna_build_report_tables
#   - per-module summary:  wgcna_build_summary_params
#   - methods section:     wgcna_build_methods
#
# Verbal-label rendering uses omicsai::omicsai_verbalize_{r,q,logfc} — keep
# them in the prompt data block; do not re-introduce raw numbers in rendered
# text except where explicitly documented (lead module's top correlated trait).

BOARD_PROMPTS_DIR <- file.path(OPG, "components/board.wgcna/prompts")


# -----------------------------------------------------------------------------
# Extraction primitives
# -----------------------------------------------------------------------------

#' Look up function descriptions for features (NULL-safe)
resolve_functions <- function(features, annot, max_chars = 60L) {
  funcs <- rep("", length(features))
  if (is.null(annot)) return(funcs)
  func_col <- intersect(c("gene_title", "gene_name", "description"), colnames(annot))
  if (length(func_col) == 0) return(funcs)
  idx <- match(features, rownames(annot))
  valid <- !is.na(idx)
  funcs[valid] <- as.character(annot[idx[valid], func_col[1]])
  substr(funcs, 1, max_chars)
}

#' Render module identifiers in canonical MEcolor form.
#'
#' Datasets keyed by integer suffix (`ME0`/`ME1`/...) get remapped via
#' `WGCNA::labels2colors`; datasets already in MEcolor form pass through.
#' Identifiers not present in `wgcna$me.genes` are returned unchanged.
display_module <- function(x, wgcna) {
  me_names <- names(wgcna$me.genes)
  suf <- suppressWarnings(as.integer(sub("^ME", "", me_names)))
  if (length(suf) == 0 || any(is.na(suf))) return(x)
  map <- setNames(paste0("ME", WGCNA::labels2colors(suf)), me_names)
  out <- map[x]
  out[is.na(out)] <- x[is.na(out)]
  unname(out)
}


#' Extract all data for a single WGCNA module.
#'
#' Returns a list with raw values used by both the summary mode and the
#' multi-module report mode. The caller is expected to cache `top` and `M`
#' across modules.
#'
#' Handles two `wgcna$gse` shapes: per-module list (mox-brca, Bruker_Generoso)
#' and flat data.frame with a `module` column (example-data and similar).
extract_module_data <- function(wgcna, module, pgx,
                                top = NULL, M = NULL,
                                max_traits = 3L) {
  annot <- wgcna$annot

  if (is.null(top)) {
    top <- tryCatch(
      playbase::wgcna.getTopGenesAndSets(wgcna, annot = annot, ntop = 40,
                                         level = "gene", rename = "gene_title"),
      error = function(e) list(pheno = list(), genes = list(), sets = list())
    )
  }

  phenotypes <- "None"
  if (module %in% names(top$pheno)) {
    phenotypes <- paste(top$pheno[[module]], collapse = ", ")
  }
  fallback_genes <- top$genes[[module]]
  if (is.null(fallback_genes)) fallback_genes <- character(0)
  fallback_sets <- top$sets[[module]]
  if (is.null(fallback_sets)) fallback_sets <- character(0)

  size <- length(wgcna$me.genes[[module]])

  if (is.null(M)) {
    M <- tryCatch(playbase:::wgcna.get_modTraits(wgcna), error = function(e) NULL)
  }

  top_trait <- ""
  top_r <- NA_real_
  top_pos_trait <- NA_character_
  top_pos_r <- NA_real_
  top_neg_trait <- NA_character_
  top_neg_r <- NA_real_
  if (!is.null(M) && module %in% rownames(M)) {
    cors <- M[module, ]
    top_idx <- which.max(abs(cors))
    if (length(top_idx) > 0) {
      top_trait <- names(cors)[top_idx]
      top_r <- cors[top_idx]
    }
    ## Dual-trait split: best positive + best negative trait, each gated by
    ## |r| >= 0.5 so direction cannot be inferred from a single absolute-value
    ## pick (the silent contradiction we hit in v01).
    pos_idx <- which.max(cors)
    if (length(pos_idx) > 0 && !is.na(cors[pos_idx]) && cors[pos_idx] >= 0.5) {
      top_pos_trait <- names(cors)[pos_idx]
      top_pos_r     <- cors[pos_idx]
    }
    neg_idx <- which.min(cors)
    if (length(neg_idx) > 0 && !is.na(cors[neg_idx]) && cors[neg_idx] <= -0.5) {
      top_neg_trait <- names(cors)[neg_idx]
      top_neg_r     <- cors[neg_idx]
    }
  }

  traits <- character(0)
  if (phenotypes != "None" && nzchar(phenotypes)) {
    traits <- trimws(strsplit(phenotypes, ",\\s*")[[1]])
    traits <- traits[nzchar(traits) & traits != "None"]
    traits <- head(traits, max_traits)
  }

  eigengene_profile <- NULL
  datME <- if (!is.null(wgcna$datME)) wgcna$datME else wgcna$net$MEs
  groups <- if (!is.null(pgx$samples$group)) pgx$samples$group else NULL
  if (!is.null(datME) && module %in% colnames(datME) && !is.null(groups)) {
    eigengene_profile <- tapply(datME[, module], groups, mean)
  }

  ## Enrichment: resolve both gse shapes (per-module list vs flat df).
  n_sig <- 0L
  n_total <- 0L
  gse_raw <- wgcna$gse
  gse <- NULL
  if (is.data.frame(gse_raw) && "module" %in% colnames(gse_raw)) {
    gse <- gse_raw[!is.na(gse_raw$module) & gse_raw$module == module, , drop = FALSE]
  } else if (is.list(gse_raw)) {
    gse <- gse_raw[[module]]
  }
  if (is.data.frame(gse) && nrow(gse) > 0 && "q.value" %in% colnames(gse)) {
    n_sig <- sum(gse$q.value < 0.05, na.rm = TRUE)
    n_total <- nrow(gse)
  }

  trait_stats <- list()
  stat_traits <- if (length(traits) > 0) traits else if (nzchar(top_trait)) top_trait
  for (tr in stat_traits) {
    gs <- tryCatch(
      playbase::wgcna.getGeneStats(wgcna, module = module, trait = tr, plot = FALSE),
      error = function(e) NULL
    )
    if (is.data.frame(gs) && nrow(gs) > 0) {
      trait_stats[[tr]] <- gs
    }
  }

  list(
    size = size,
    phenotypes = phenotypes,
    traits = traits,
    fallback_genes = fallback_genes,
    fallback_sets = fallback_sets,
    eigengene_profile = eigengene_profile,
    top_trait = top_trait,
    top_r = top_r,
    top_pos_trait = top_pos_trait,
    top_pos_r = top_pos_r,
    top_neg_trait = top_neg_trait,
    top_neg_r = top_neg_r,
    n_sig = n_sig,
    n_total = n_total,
    gse = gse,
    trait_stats = trait_stats
  )
}


#' Compute per-module structured data for a set of modules.
#'
#' Wraps extract_module_data() over a list of modules, attaching hub-gene
#' tables (resolved symbols + descriptions) and the trait selected for the
#' single-trait report view.
.compute_module_data <- function(wgcna, pgx, modules,
                                 ntop_enrichment = 20L, ntop_genes = 50L) {
  annot <- if (!is.null(wgcna$annot)) wgcna$annot else pgx$genes

  M <- tryCatch(playbase:::wgcna.get_modTraits(wgcna), error = function(e) NULL)
  top <- tryCatch(
    playbase::wgcna.getTopGenesAndSets(wgcna, annot = annot, ntop = 40,
                                       level = "gene", rename = "gene_title"),
    error = function(e) list(pheno = list(), genes = list(), sets = list())
  )

  module_data <- list()
  for (mod in modules) {
    md <- extract_module_data(wgcna, mod, pgx, top = top, M = M)

    ## Hub genes for the single-trait report view (first phenotype-derived
    ## trait if available, else the absolute-top trait).
    hub_genes_df <- NULL
    trait_for_genes <- if (length(md$traits) > 0) md$traits[1] else md$top_trait

    gs <- if (nzchar(trait_for_genes) && trait_for_genes %in% names(md$trait_stats)) {
      md$trait_stats[[trait_for_genes]]
    } else {
      NULL
    }

    if (is.data.frame(gs) && nrow(gs) > 0) {
      if ("moduleMembership" %in% colnames(gs)) {
        gs <- gs[order(-abs(gs$moduleMembership)), , drop = FALSE]
      }
      top_gs <- head(gs, ntop_genes)

      symbols <- playbase::probe2symbol(top_gs$feature, annot, "symbol", fill_na = TRUE)
      funcs <- resolve_functions(top_gs$feature, annot)

      mm_vals <- if ("moduleMembership" %in% colnames(top_gs)) {
        top_gs$moduleMembership
      } else {
        rep(NA_real_, nrow(top_gs))
      }
      logfc_vals <- if ("foldChange" %in% colnames(top_gs)) {
        top_gs$foldChange
      } else {
        rep(NA_real_, nrow(top_gs))
      }

      hub_genes_df <- data.frame(
        symbol = symbols,
        MM = mm_vals,
        logFC = logfc_vals,
        func = funcs,
        stringsAsFactors = FALSE
      )
    }

    module_data[[mod]] <- list(
      size = md$size,
      eigengene_profile = md$eigengene_profile,
      top_trait = md$top_trait,
      top_r = md$top_r,
      top_pos_trait = md$top_pos_trait,
      top_pos_r = md$top_pos_r,
      top_neg_trait = md$top_neg_trait,
      top_neg_r = md$top_neg_r,
      n_sig = md$n_sig,
      n_total = md$n_total,
      full_gse = md$gse,
      hub_genes = hub_genes_df,
      trait_for_genes = trait_for_genes
    )
  }

  module_data
}


#' Per-module gene-family enrichment (≥3 members in a family of 5–500).
#' Returns a named list of character vectors, keyed by module.
.compute_families_text <- function(wgcna, pgx, modules) {
  if (is.null(pgx$families)) return(NULL)
  fam_names <- names(pgx$families)
  fam_names <- fam_names[fam_names != "<all>"]
  fam_sizes <- lengths(pgx$families[fam_names])
  fam_names <- fam_names[fam_sizes >= 5 & fam_sizes <= 500]
  if (length(fam_names) == 0) return(NULL)

  out <- list()
  for (mod in modules) {
    mod_genes <- wgcna$me.genes[[mod]]
    if (is.null(mod_genes) || length(mod_genes) == 0) next

    overlaps <- vapply(fam_names, function(fn) {
      length(intersect(mod_genes, pgx$families[[fn]]))
    }, integer(1))

    sig_fam <- overlaps[overlaps >= 3]
    if (length(sig_fam) == 0) next
    sig_fam <- head(sort(sig_fam, decreasing = TRUE), 5)

    out[[mod]] <- vapply(names(sig_fam), function(fn) {
      sprintf("%s (%d of %d)", fn, sig_fam[fn], length(pgx$families[[fn]]))
    }, character(1))
  }
  out
}


# -----------------------------------------------------------------------------
# Section builders — one per template slot
# -----------------------------------------------------------------------------

#' Build the `## Overview` placeholder values.
#'
#' Returns a list of named substitutions for the Overview block of the
#' data template. Caller merges these into the substitute() arg list.
data_overview <- function(wgcna, pgx) {
  experiment <- if (!is.null(wgcna$experiment)) wgcna$experiment
                else if (!is.null(pgx$description)) pgx$description
                else "omics experiment"
  organism <- if (!is.null(pgx$organism)) pgx$organism else "unknown"
  n_samples <- tryCatch(nrow(pgx$samples), error = function(e) NA_integer_)
  n_features_total <- tryCatch(nrow(pgx$X), error = function(e) NA_integer_)
  n_features_used <- length(unlist(wgcna$me.genes))
  if (is.null(n_features_used) || n_features_used == 0) {
    n_features_used <- tryCatch(ncol(wgcna$datExpr), error = function(e) NA_integer_)
  }

  power <- if (!is.null(wgcna$power)) wgcna$power
           else if (!is.null(wgcna$net$power)) wgcna$net$power
           else "NA"
  min_mod_size <- if (!is.null(wgcna$minModSize)) wgcna$minModSize else "20"
  merge_cut_height <- if (!is.null(wgcna$mergeCutHeight)) wgcna$mergeCutHeight else "0.15"

  list(
    experiment       = as.character(experiment),
    organism         = as.character(organism),
    n_samples        = as.character(n_samples),
    n_features_total = as.character(n_features_total),
    n_features_used  = as.character(n_features_used),
    power            = as.character(power),
    min_mod_size     = as.character(min_mod_size),
    merge_cut_height = as.character(merge_cut_height)
  )
}


#' Build the `## Experimental contrasts` block.
#'
#' Returns a markdown bullet list (one line per contrast) or a placeholder
#' string if no contrast matrix is available.
data_contrast <- function(pgx) {
  cm <- pgx$model.parameters$contr.matrix
  if (is.null(cm)) return("(no contrasts available)")

  group_names <- rownames(cm)
  lines <- vapply(colnames(cm), function(cn) {
    pos <- group_names[cm[, cn] > 0]
    neg <- group_names[cm[, cn] < 0]
    sprintf("- %s: %s vs %s", cn,
            paste(pos, collapse = "+"), paste(neg, collapse = "+"))
  }, character(1))
  paste(lines, collapse = "\n")
}


#' Build the `## Modules summary` table block.
#'
#' Dual-trait columns: best positive (>= 0.5) and best negative (<= -0.5).
#' Lead module (first row) also carries `(r = ±0.NN)` next to the verbal
#' label on its top correlated trait. Grey is appended as a final row.
data_modsummary <- function(wgcna, module_data, lead_module) {
  module_order <- names(module_data)
  if (length(module_order) == 0) return(list(table = "", footnote = ""))

  fmt_trait_col <- function(trait, r, is_lead_pos = FALSE) {
    if (is.na(trait) || !nzchar(trait)) return("—")
    verbal <- omicsai::omicsai_verbalize_r(r)
    if (is_lead_pos && !is.na(r)) {
      sprintf("%s (%s, r = %+.2f)", trait, verbal, r)
    } else {
      sprintf("%s (%s)", trait, verbal)
    }
  }

  ov_rows <- lapply(module_order, function(m) {
    md <- module_data[[m]]
    is_lead <- identical(m, lead_module)
    data.frame(
      Module = display_module(m, wgcna),
      Genes  = as.character(md$size),
      `Top correlated trait`      = fmt_trait_col(md$top_pos_trait, md$top_pos_r, is_lead),
      `Top anti-correlated trait` = fmt_trait_col(md$top_neg_trait, md$top_neg_r),
      `Enrichment hits`           = as.character(md$n_sig),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  grey_size <- length(wgcna$me.genes[["MEgrey"]])
  if (is.null(grey_size) || grey_size == 0) {
    grey_size <- length(wgcna$me.genes[["grey"]])
  }
  if (!is.null(grey_size) && grey_size > 0) {
    ov_rows[[length(ov_rows) + 1]] <- data.frame(
      Module = "grey", Genes = as.character(grey_size),
      `Top correlated trait` = "—",
      `Top anti-correlated trait` = "—",
      `Enrichment hits` = "—",
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }

  overview_df <- do.call(rbind, ov_rows)
  table_md <- paste(omicsai::omicsai_format_mdtable(overview_df), collapse = "\n")

  ## Footnote: "Showing X of N modules."
  raw <- names(wgcna$me.genes)
  suf <- suppressWarnings(as.integer(sub("^ME", "", raw)))
  display <- raw
  if (length(suf) > 0 && all(!is.na(suf))) {
    display <- paste0("ME", WGCNA::labels2colors(suf))
  }
  is_grey <- raw %in% c("MEgrey", "grey") | display %in% c("MEgrey", "grey")
  n_total_non_grey <- sum(!is_grey)
  footnote <- sprintf("Showing %d of %d modules.",
                      length(module_order), n_total_non_grey)

  list(table = table_md, footnote = footnote)
}


#' Build the `## Module-module eigengene correlations` block.
#'
#' Verbal labels only; raw r stripped at the rendering layer.
data_eigen_cor <- function(wgcna, modules) {
  me_data <- if (!is.null(wgcna$datME)) wgcna$datME else wgcna$net$MEs
  if (is.null(me_data)) return("(no strong correlations detected)")

  me_cols <- intersect(modules, colnames(me_data))
  if (length(me_cols) < 2) return("(no strong correlations detected)")

  me_cor <- cor(me_data[, me_cols])
  pairs <- character(0)
  for (i in seq_len(length(me_cols) - 1)) {
    for (j in (i + 1):length(me_cols)) {
      r <- me_cor[i, j]
      if (!is.na(r) && abs(r) >= 0.7) {
        pairs <- c(pairs, sprintf("%s ↔ %s: %s",
                                  display_module(me_cols[i], wgcna),
                                  display_module(me_cols[j], wgcna),
                                  omicsai::omicsai_verbalize_r(r)))
      }
    }
  }
  if (length(pairs) == 0) return("(no strong correlations detected)")
  paste(pairs, collapse = "\n")
}


#' Build per-module detail blocks (v03 layout).
#'
#' For each module: heading, eigengene profile, dual-trait coordination,
#' hub genes (BEFORE enrichment), enrichment top-N with Hub overlap column,
#' optional gene-family enrichment.
data_module_detail <- function(wgcna, pgx, module_data,
                               families_text = NULL,
                               ntop_enrichment = 20L,
                               ntop_genes = 10L) {
  module_order <- names(module_data)
  if (length(module_order) == 0) return("")

  blocks <- vapply(module_order, function(mod) {
    .render_module_block(module_data[[mod]], mod, wgcna,
                         families_text   = families_text,
                         ntop_enrichment = ntop_enrichment,
                         ntop_genes      = ntop_genes)
  }, character(1))

  paste(blocks, collapse = "\n\n")
}


# -----------------------------------------------------------------------------
# Orchestrator
# -----------------------------------------------------------------------------

#' Build structured report tables from WGCNA results.
#'
#' Renders the v03 data block: a single markdown document substituted into
#' `prompts/wgcna_report_data.md`, plus a structured `data` list for any
#' downstream callers that want the raw values.
#'
#' @return list(text = character, data = list)
wgcna_build_report_tables <- function(wgcna, pgx,
                                      n_modules = 8L,
                                      ntop_enrichment = 20L,
                                      ntop_genes = 50L,
                                      include_module_cors = TRUE,
                                      include_families = TRUE,
                                      include_contrasts = TRUE) {
  ## Lazy-fill labels + stats so older PGX files extract cleanly.
  wgcna <- playbase::wgcna.ensureStats(wgcna)

  modules <- playbase::wgcna.getTopModules(wgcna, min_modules = 5L)
  modules <- head(modules, n_modules)

  module_data <- .compute_module_data(wgcna, pgx, modules,
                                      ntop_enrichment = ntop_enrichment,
                                      ntop_genes = ntop_genes)

  ## Order modules by significant-enrichment count, descending.
  module_order <- names(sort(
    vapply(module_data, function(x) x$n_sig, integer(1)),
    decreasing = TRUE
  ))
  module_data <- module_data[module_order]
  lead_module <- if (length(module_order) > 0) module_order[1] else NA_character_

  families_text <- if (include_families) {
    .compute_families_text(wgcna, pgx, module_order)
  } else NULL

  contrasts_block <- if (include_contrasts) data_contrast(pgx)
                     else "(contrasts omitted)"
  module_cors <- if (include_module_cors) data_eigen_cor(wgcna, module_order)
                 else "(module-module correlations omitted)"

  overview_params <- data_overview(wgcna, pgx)
  modsum <- data_modsummary(wgcna, module_data, lead_module)
  per_module <- data_module_detail(wgcna, pgx, module_data,
                                   families_text = families_text,
                                   ntop_enrichment = ntop_enrichment,
                                   ntop_genes = ntop_genes)

  tmpl <- omicsai::omicsai_load_template(
    file.path(BOARD_PROMPTS_DIR, "wgcna_report_data.md")
  )

  text <- omicsai::omicsai_substitute_template(tmpl, c(
    overview_params,
    list(
      contrasts_block            = contrasts_block,
      modules_summary_table      = modsum$table,
      modules_summary_footnote   = modsum$footnote,
      module_module_correlations = module_cors,
      module_detail              = per_module
    )
  ))

  list(
    text = text,
    data = list(
      experiment       = overview_params$experiment,
      organism         = overview_params$organism,
      n_samples        = overview_params$n_samples,
      n_features_total = overview_params$n_features_total,
      n_features_used  = overview_params$n_features_used,
      modules          = module_data
    )
  )
}


## MODULE_SUMMARY — verbatim v03 prompt-optim spec, loaded from disk.
## Source of truth lives at `prompts/wgcna_module_data.md`, copied
## word-for-word from `tmp/.../candidate/wgcna/data/v03_prompttemplates/data.md`.
## Edit the file directly to change wording; the R renderer only fills values.
MODULE_SUMMARY <- omicsai::omicsai_load_template(
  file.path(BOARD_PROMPTS_DIR, "wgcna_module_data.md")
)


# -- Leaf renderers: one per template placeholder that needs derivation. -----
# Each returns a single string. No prose, only data → string formatting.

.eigengene_profile_str <- function(md) {
  if (is.null(md$eigengene_profile)) return("—")
  paste(sprintf("%s=%s", names(md$eigengene_profile),
                omicsai::omicsai_format_num(md$eigengene_profile, 2)),
        collapse = ", ")
}

.hub_genes_sorted <- function(md, ntop) {
  if (is.null(md$hub_genes) || nrow(md$hub_genes) == 0) return(NULL)
  hg <- md$hub_genes
  mm  <- if ("MM"    %in% colnames(hg)) hg$MM    else rep(NA_real_, nrow(hg))
  lfc <- if ("logFC" %in% colnames(hg)) hg$logFC else rep(NA_real_, nrow(hg))
  head(hg[order(-abs(mm), -abs(lfc), na.last = TRUE), , drop = FALSE], ntop)
}

.hub_genes_table_str <- function(hg) {
  if (is.null(hg) || nrow(hg) == 0) return("—")
  paste(omicsai::omicsai_format_mdtable(data.frame(
    Gene             = paste0("*", hg$symbol, "*"),
    `Known function` = hg$func,
    stringsAsFactors = FALSE, check.names = FALSE
  )), collapse = "\n")
}

.enrichment_table_str <- function(md, hub_symbols, ntop) {
  full_gse <- md$full_gse
  if (!is.data.frame(full_gse) || nrow(full_gse) == 0) return("—")

  ## Sort by overlap (preferred), then q-value, then score.
  full_gse <- if ("overlap" %in% colnames(full_gse)) {
    full_gse[order(full_gse$overlap, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  } else if ("overlap_count" %in% colnames(full_gse)) {
    full_gse[order(full_gse$overlap_count, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  } else if ("q.value" %in% colnames(full_gse)) {
    full_gse[order(full_gse$q.value, na.last = TRUE), , drop = FALSE]
  } else if ("score" %in% colnames(full_gse)) {
    full_gse[order(full_gse$score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  } else full_gse

  top_n <- head(full_gse, ntop)
  hub_overlap_col <- if (length(hub_symbols) == 0) {
    rep("—", nrow(top_n))
  } else if ("genes" %in% colnames(top_n)) {
    vapply(top_n$genes, function(g) {
      gset <- strsplit(as.character(g), "|", fixed = TRUE)[[1]]
      sprintf("%d/%d", length(intersect(gset, hub_symbols)), length(hub_symbols))
    }, character(1), USE.NAMES = FALSE)
  } else {
    rep(sprintf("0/%d", length(hub_symbols)), nrow(top_n))
  }

  paste(omicsai::omicsai_format_mdtable(data.frame(
    Rank         = seq_len(nrow(top_n)),
    Geneset      = as.character(top_n$geneset),
    `Hub overlap` = hub_overlap_col,
    stringsAsFactors = FALSE, check.names = FALSE
  )), collapse = "\n")
}


#' Render one module's data dict into a MODULE_SUMMARY block.
#'
#' Shared by `data_module_detail()` (Report mode, looped) and
#' `wgcna_build_summary_params()` (Summary mode, single module). The `md`
#' argument is one element from `.compute_module_data()`'s output. This
#' function builds only the value dict; MODULE_SUMMARY owns all prose.
.render_module_block <- function(md, mod_name, wgcna,
                                 families_text = NULL,
                                 ntop_enrichment = 20L,
                                 ntop_genes = 10L) {
  hg <- .hub_genes_sorted(md, ntop_genes)
  hub_symbols <- if (!is.null(hg)) as.character(hg$symbol) else character(0)

  fam <- if (!is.null(families_text) && !is.null(families_text[[mod_name]])) {
    paste(families_text[[mod_name]], collapse = "; ")
  } else "—"

  overlaps <- if (!is.null(md$enrichment_overlap) && length(md$enrichment_overlap) > 0) {
    paste(md$enrichment_overlap, collapse = "; ")
  } else "—"

  na_dash <- function(x) if (is.na(x) || !nzchar(x)) "—" else x

  omicsai::omicsai_substitute_template(MODULE_SUMMARY, list(
    ME_color                      = display_module(mod_name, wgcna),
    n_genes                       = as.character(md$size),
    tier                          = "—",                ## tier classification dropped (epic playbase-fad)
    eigengene_profile_qualitative = .eigengene_profile_str(md),
    top_pos_trait                 = na_dash(md$top_pos_trait),
    top_pos_verbal                = omicsai::omicsai_verbalize_r(md$top_pos_r),
    top_neg_trait                 = na_dash(md$top_neg_trait),
    top_neg_verbal                = omicsai::omicsai_verbalize_r(md$top_neg_r),
    n_sig_terms                   = as.character(md$n_sig),
    n_total_terms                 = as.character(md$n_total),
    enrichment_themes_table       = .enrichment_table_str(md, hub_symbols, ntop_enrichment),
    n_hub                         = as.character(if (is.null(hg)) 0L else nrow(hg)),
    hub_genes_table               = .hub_genes_table_str(hg),
    gene_families_summary         = fam,
    enrichment_overlaps           = overlaps
  ))
}


# -----------------------------------------------------------------------------
# Per-module summary mode (single-module, used by the Summary tab)
# -----------------------------------------------------------------------------

#' Build prompt parameters for a single WGCNA module summary.
#'
#' Produces the same per-module block shape as Report mode by routing through
#' `.compute_module_data()` and `.render_module_block()`. The downstream
#' `{{module_detail}}` placeholder in `wgcna_summary.md` receives the
#' rendered MODULE_SUMMARY block.
#'
#' @return Named list: experiment, module, module_detail.
wgcna_build_summary_params <- function(wgcna, module, pgx) {
  wgcna <- playbase::wgcna.ensureStats(wgcna)

  module_data <- .compute_module_data(wgcna, pgx, modules = module)
  md <- module_data[[module]]

  module_detail <- .render_module_block(md, module, wgcna)

  list(
    experiment    = if (!is.null(wgcna$experiment)) wgcna$experiment else "",
    module        = module,
    module_detail = module_detail
  )
}


# -----------------------------------------------------------------------------
# Methods section (deterministic appendix)
# -----------------------------------------------------------------------------

#' Build deterministic methods section for WGCNA report.
wgcna_build_methods <- function(wgcna, pgx) {
  template_path <- file.path(BOARD_PROMPTS_DIR, "wgcna_methods.md")
  template <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

  n_wgcna <- length(unlist(wgcna$me.genes))
  if (is.null(n_wgcna) || n_wgcna == 0) {
    n_wgcna <- tryCatch(ncol(wgcna$datExpr), error = function(e) NA)
  }

  feature_type <- "features"
  if (!is.null(pgx$datatype)) {
    if (grepl("prot", pgx$datatype, ignore.case = TRUE)) {
      feature_type <- "proteins"
    } else if (grepl("rna|transcript|gene", pgx$datatype, ignore.case = TRUE)) {
      feature_type <- "genes"
    }
  }

  n_samples <- tryCatch(nrow(wgcna$datExpr), error = function(e) {
    tryCatch(nrow(pgx$samples), error = function(e2) NA)
  })

  network_type <- if (!is.null(wgcna$networktype)) wgcna$networktype
                  else if (!is.null(wgcna$net$networkType)) wgcna$net$networkType
                  else "signed"
  power <- if (!is.null(wgcna$power)) wgcna$power
           else if (!is.null(wgcna$net$power)) wgcna$net$power
           else "NA"
  min_mod_size <- if (!is.null(wgcna$minModSize)) wgcna$minModSize else "20"
  merge_cut_height <- if (!is.null(wgcna$mergeCutHeight)) wgcna$mergeCutHeight else "0.15"
  min_kme <- if (!is.null(wgcna$minKME)) wgcna$minKME else "0.3"

  modules_no_grey <- setdiff(names(wgcna$me.genes), c("grey", "MEgrey"))
  n_modules <- length(modules_no_grey)
  grey_size <- length(wgcna$me.genes[["MEgrey"]]) + length(wgcna$me.genes[["grey"]])

  n_genesets <- tryCatch({
    first_mod <- modules_no_grey[1]
    nrow(wgcna$gse[[first_mod]])
  }, error = function(e) "NA")

  params <- list(
    n_features_wgcna = as.character(n_wgcna),
    feature_type = feature_type,
    n_samples = as.character(n_samples),
    network_type = as.character(network_type),
    power = as.character(power),
    min_mod_size = as.character(min_mod_size),
    merge_cut_height = as.character(merge_cut_height),
    min_kme = as.character(min_kme),
    n_modules = as.character(n_modules),
    grey_size = as.character(grey_size),
    n_genesets_tested = as.character(n_genesets),
    date = format(Sys.Date(), "%Y-%m-%d")
  )

  omicsai::omicsai_substitute_template(template, params)
}
