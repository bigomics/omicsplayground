## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

# =============================================================================
# AI Text Server for WGCNA
# =============================================================================
# This module handles text generation for both summary and report modes in
# WGCNA AI reports. It provides functions to extract parameters, generate
# summaries, and create integrated reports.

# =============================================================================
# SECTION 1: Pure Functions (No Shiny Dependencies)
# =============================================================================

#' Build AI prompt parameters for WGCNA module summary
#'
#' Extracts parameters from WGCNA results to populate the prompt template.
#' Replicates exact extraction logic from playbase::wgcna.describeModules()
#' (see playbase/R/pgx-wgcna.R lines 5340-5361)
#'
#' @param wgcna WGCNA results object
#' @param module Character; module name (e.g., "blue", "turquoise")
#' @param annot Data frame; gene annotation (optional)
#' @param ntop Integer; number of top genes/sets to include (default 40)
#' @param multi Logical; TRUE for multi-omics WGCNA
#'
#' @return Named list with template parameters (board_params):
#'   module, phenotypes, experiment, genesets, keygenes_section
wgcna_build_ai_params <- function(wgcna,
                                  module,
                                  annot = NULL,
                                  ntop = 40,
                                  multi = FALSE) {

  # Get top genes and sets using playbase function

# (mirrors playbase/R/pgx-wgcna.R lines 5280-5288)
  if (multi) {
    top <- playbase::wgcna.getMultiTopGenesAndSets(
      wgcna,
      annot = annot,
      ntop = ntop,
      level = "gene",
      rename = "gene_title"
    )
  } else {
    top <- playbase::wgcna.getTopGenesAndSets(
      wgcna,
      annot = annot,
      ntop = ntop,
      level = "gene",
      rename = "gene_title"
    )
  }

  # Extract phenotypes - format as comma-separated list
  pp <- "None"
  if (module %in% names(top$pheno)) {
    pp <- paste(top$pheno[[module]], collapse = ", ")
  }

  # Extract genesets with quantitative metrics (fallback to names only)
  ss <- ""
  gse <- NULL
  if (!is.null(wgcna$gse) && !is.null(wgcna$gse[[module]])) {
    gse <- wgcna$gse[[module]]
  }
  if (is.data.frame(gse) && nrow(gse) > 0 &&
    all(c("geneset", "score", "q.value", "overlap") %in% colnames(gse))) {
    gse <- gse[order(gse$q.value, gse$score, na.last = TRUE), , drop = FALSE]
    sig <- gse[!is.na(gse$q.value) & gse$q.value < 0.05, , drop = FALSE]
    top_gse <- head(sig, 8)
    if (nrow(top_gse) > 0) {
      gset_df <- data.frame(
        Pathway = top_gse$geneset,
        Score = top_gse$score,
        `Q-value` = top_gse$q.value,
        Overlap = top_gse$overlap,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      gset_table <- omicsai::omicsai_format_mdtable(gset_df, formatters = list(
        Score = function(x) omicsai::omicsai_format_num(x, 2),
        `Q-value` = function(x) omicsai::omicsai_format_num(x, 3)
      ))

      score_vals <- gse$score[!is.na(gse$score)]
      score_range <- if (length(score_vals) > 0) {
        paste0(
          omicsai::omicsai_format_num(min(score_vals), 2), "-",
          omicsai::omicsai_format_num(max(score_vals), 2),
          " (median: ", omicsai::omicsai_format_num(median(score_vals), 2), ")"
        )
      } else {
        "NA"
      }
      n_sig <- sum(gse$q.value < 0.05, na.rm = TRUE)

      ss <- omicsai::collapse_lines(
        "**Top Enriched Pathways (q < 0.05):**",
        gset_table,
        paste0(
          "**Score range:** ", score_range, "  \n",
          "**Total significant pathways:** ", n_sig, " of ", nrow(gse), " tested"
        ),
        sep = "\n\n"
      )
    }
  }
  if (ss == "" && !is.null(top$sets[[module]])) {
    pathways <- top$sets[[module]]
    ss <- paste0("- ", pathways, collapse = "\n")
  }

  # Extract key genes with metrics — multi-trait approach
  # WGCNA modules are trait-independent (MM, centrality are network properties).
  # TS and logFC are trait-dependent, so we compute them for ALL correlated traits
  # to give the AI a complete picture across phenotypes.
  keygenes_section <- ""
  traits <- NULL
  if (is.character(pp) && pp != "None" && nzchar(pp)) {
    traits <- trimws(strsplit(pp, ",\\s*")[[1]])
    traits <- traits[nzchar(traits) & traits != "None"]
    traits <- head(traits, 3) # cap to avoid very wide tables
  }

  # Collect gene stats for each correlated trait
  trait_stats_list <- list()
  if (length(traits) > 0) {
    for (tr in traits) {
      gs <- tryCatch(
        playbase::wgcna.getGeneStats(
          wgcna,
          module = module,
          trait = tr,
          plot = FALSE
        ),
        error = function(e) NULL
      )
      if (is.data.frame(gs) && nrow(gs) > 0) {
        trait_stats_list[[tr]] <- gs
      }
    }
  }

  if (length(trait_stats_list) > 0) {
    # Use first trait's stats as base for shared columns (MM, centrality)
    base_stats <- trait_stats_list[[1]]
    if (!"moduleMembership" %in% colnames(base_stats)) {
      base_stats$moduleMembership <- NA_real_
    }
    if (!"centrality" %in% colnames(base_stats)) {
      base_stats$centrality <- NA_real_
    }

    # Sort by module membership (trait-independent)
    if ("moduleMembership" %in% colnames(base_stats)) {
      base_stats <- base_stats[order(-abs(base_stats$moduleMembership)), , drop = FALSE]
    } else if ("score" %in% colnames(base_stats)) {
      base_stats <- base_stats[order(-base_stats$score), , drop = FALSE]
    }

    top_genes <- head(base_stats, 8)
    features <- top_genes$feature

    # Resolve gene symbols
    symbols <- features
    if (!is.null(annot)) {
      symbols <- tryCatch(
        playbase::probe2symbol(features, annot, "symbol"),
        error = function(e) features
      )
      symbols <- ifelse(is.na(symbols) | symbols == "", features, symbols)
    }

    # --- Table 1: trait-independent network metrics ---
    network_df <- data.frame(
      Gene = symbols,
      MM = top_genes$moduleMembership,
      Centrality = top_genes$centrality,
      stringsAsFactors = FALSE
    )
    network_table <- omicsai::collapse_lines(
      "### Trait-independent metrics (network structure)",
      omicsai::omicsai_format_mdtable(network_df, formatters = list(
        MM = function(x) omicsai::omicsai_format_num(x, 2),
        Centrality = function(x) omicsai::omicsai_format_num(x, 2)
      )),
      sep = "\n\n"
    )

    # --- Tables 2..N: per-trait gene-trait associations ---
    trait_names <- names(trait_stats_list)
    M <- tryCatch(playbase::wgcna.get_modTraits(wgcna), error = function(e) NULL)

    trait_tables <- vapply(trait_names, function(tr) {
      gs <- trait_stats_list[[tr]]

      # Module-trait correlation for header
      mod_cor <- NA_real_
      if (is.matrix(M) || is.data.frame(M)) {
        mi <- which(rownames(M) %in% c(module, paste0("ME", module)))
        ti <- which(colnames(M) == tr)
        if (length(mi) > 0 && length(ti) > 0) mod_cor <- M[mi[1], ti[1]]
      }
      cor_label <- if (!is.na(mod_cor)) {
        paste0(" (module r = ", omicsai::omicsai_format_num(mod_cor, 2), ")")
      } else {
        ""
      }

      # Build data for each hub gene
      ts_vals <- vapply(features, function(f) {
        idx <- match(f, gs$feature)
        if (!is.na(idx) && "traitSignificance" %in% colnames(gs)) gs$traitSignificance[idx] else NA_real_
      }, numeric(1))
      fc_vals <- vapply(features, function(f) {
        idx <- match(f, gs$feature)
        if (!is.na(idx) && "foldChange" %in% colnames(gs)) gs$foldChange[idx] else NA_real_
      }, numeric(1))

      trait_df <- data.frame(
        Gene = symbols, TS = ts_vals, logFC = fc_vals,
        stringsAsFactors = FALSE
      )
      omicsai::collapse_lines(
        paste0("### Gene-trait associations: ", tr, cor_label),
        omicsai::omicsai_format_mdtable(trait_df, formatters = list(
          TS = function(x) omicsai::omicsai_format_num(x, 2),
          logFC = function(x) omicsai::omicsai_format_num(x, 1)
        )),
        sep = "\n\n"
      )
    }, character(1))

    # Module size
    module_size <- NA_integer_
    if (!is.null(wgcna$me.genes) && !is.null(wgcna$me.genes[[module]])) {
      module_size <- length(wgcna$me.genes[[module]])
    } else {
      module_size <- nrow(base_stats)
    }

    keygenes_section <- omicsai::collapse_lines(
      paste0(
        "**Hub Genes (ranked by module membership, top ",
        nrow(top_genes), " of ",
        if (!is.na(module_size)) module_size else nrow(base_stats),
        " module genes):**"
      ),
      network_table,
      paste(trait_tables, collapse = "\n\n"),
      paste(
        "**MM:** Module Membership (correlation with eigengene) — trait-independent  ",
        "**Centrality:** Intramodular connectivity — trait-independent  ",
        "**TS:** Trait Significance (gene-trait correlation) — specific to each trait  ",
        "**logFC:** Log2 fold change — specific to each trait",
        sep = "\n"
      ),
      sep = "\n\n"
    )
  } else if (!is.null(top$genes[[module]]) && length(top$genes[[module]]) > 0) {
    genes <- head(top$genes[[module]], 15)
    keygenes_section <- paste0(
      "The following hub genes show high intramodular connectivity:\n\n",
      paste(genes, collapse = ", ")
    )
  }

  # Module-level summary statistics (multi-trait)
  module_stats <- ""
  module_size <- NA_integer_
  if (!is.null(wgcna$me.genes) && !is.null(wgcna$me.genes[[module]])) {
    module_size <- length(wgcna$me.genes[[module]])
  } else if (length(trait_stats_list) > 0) {
    module_size <- nrow(trait_stats_list[[1]])
  }

  stat_lines <- c("**Module Statistics:**")
  if (!is.na(module_size)) {
    stat_lines <- c(stat_lines, paste0("- **Size:** ", module_size, " genes"))
  }

  # Per-trait: module-trait correlation and mean fold change
  if (length(traits) > 0) {
    M <- tryCatch(playbase::wgcna.get_modTraits(wgcna), error = function(e) NULL)
    for (tr in traits) {
      trait_cor <- NA_real_
      if (is.matrix(M) || is.data.frame(M)) {
        module_idx <- which(rownames(M) %in% c(module, paste0("ME", module)))
        trait_idx <- which(colnames(M) == tr)
        if (length(module_idx) > 0 && length(trait_idx) > 0) {
          trait_cor <- M[module_idx[1], trait_idx[1]]
        }
      }
      if (!is.na(trait_cor)) {
        stat_lines <- c(stat_lines, paste0(
          "- **Trait correlation with ", tr, ":** ", omicsai::omicsai_format_num(trait_cor, 2)
        ))
      }
      # Mean fold change for this trait
      if (tr %in% names(trait_stats_list)) {
        gs <- trait_stats_list[[tr]]
        if ("foldChange" %in% colnames(gs)) {
          mean_fc <- mean(gs$foldChange, na.rm = TRUE)
          if (!is.na(mean_fc)) {
            fc_dir <- ifelse(mean_fc > 0, "upregulated", "downregulated")
            stat_lines <- c(stat_lines, paste0(
              "- **Mean expression change (", tr, "):** ",
              omicsai::omicsai_format_num(abs(mean_fc), 1), "-fold ", fc_dir
            ))
          }
        }
      }
    }
  }

  module_stats <- paste(stat_lines, collapse = "\n")

  # Get experiment description
  experiment <- wgcna$experiment %||% ""

  list(
    module = module,
    phenotypes = pp,
    experiment = experiment,
    genesets = ss,
    keygenes_section = keygenes_section,
    module_stats = module_stats
  )
}

#' Generate summaries for all top WGCNA modules
#'
#' Loops over the top N modules, extracts params via wgcna_build_ai_params(),
#' renders each through omicsai_gen_text(), and returns a named list of
#' summary texts.
#'
# =============================================================================
# SECTION 2: Shiny Module
# =============================================================================

#' WGCNA AI Text Server
#'
#' Handles text generation for both summary and report modes.
#' In summary mode, generates a single-module summary via omicsai_gen_text().
#' In report mode, generates multi-module summaries then integrates via omicsai_create_report().
#'
#' @param id Module namespace ID
#' @param wgcna Reactive returning WGCNA results object
#' @param pgx PGX object (non-reactive)
#' @param controls List of control reactives (from ai_report_controls_server)
#' @param parent_session Parent Shiny session (for getUserOption)
#'
#' @return List with reactives: text (current mode), report_text (report-only)
wgcna_ai_text_server <- function(id, wgcna, pgx, controls, parent_session) {
  moduleServer(id, function(input, output, session) {

    # ---- Shared: templates, context, model ----

    ai_summary_template <- omicsai::omicsai_load_template(
      file.path(OPG, "components/board.wgcna/prompts/wgcna_summary_template.md")
    )

    ai_context_template <- omicsai::omicsai_load_template(
      file.path(OPG, "components/board.wgcna/prompts/WGCNA_methods.md")
    )

    ai_context <- shiny::reactive({
      w <- wgcna()
      shiny::req(w)
      omicsai::omicsai_substitute_template(
        ai_context_template,
        list(experiment = w$experiment %||% "omics experiment")
      )
    })

    ai_model <- shiny::reactive({
      m <- getUserOption(parent_session, "llm_model")
      shiny::req(m, m != "")
      m
    })

    # ---- Prompt caches (for instant toggle without regeneration) ----
    summary_prompt_cache <- shiny::reactiveVal(NULL)
    report_prompt_cache <- shiny::reactiveVal(NULL)

    # ---- SUMMARY MODE: per-module summary ----

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module()),
    {
      mode <- controls$mode()
      if (mode != "summary" || controls$trigger() < 1) return(NULL)

      w <- wgcna()
      module <- controls$selected_module()
      shiny::req(w, module)
      model <- ai_model()
      shiny::req(model)

      style <- controls$summary_style() %||% "short"
      params <- wgcna_build_ai_params(w, module, annot = w$annot, multi = FALSE)
      params$methods_context <- ai_context()
      params$style_instructions <- omicsai::omicsai_instructions(paste0("format_", style))

      ## Cache the prompt for instant toggle
      prompt <- omicsai::omicsai_substitute_template(ai_summary_template, params)
      summary_prompt_cache(paste0("# Summary Prompt\n\n", prompt))

      result <- omicsai::omicsai_gen_text(
        template = ai_summary_template,
        params = params,
        config = omicsai::omicsai_config(model = model %||% Sys.getenv("OMICS_AI_MODEL", "ollama:llama3.2"))
      )
      result$text
    }, ignoreNULL = FALSE)

    # ---- REPORT MODE: single-shot integrated report ----

    ai_report <- shiny::eventReactive(controls$trigger(), {
      mode <- controls$mode()
      if (mode != "report" || controls$trigger() < 1) return(NULL)

      w <- wgcna()
      shiny::req(w)
      model <- ai_model()
      shiny::req(model)

      progress <- shiny::Progress$new()
      on.exit(progress$close())

      ## Step 1: Build structured data tables
      progress$set(message = "Extracting module data...", value = 0.1)
      tables <- wgcna_build_report_tables(w, pgx)

      ## Step 2: Classify module signal strength
      progress$set(message = "Classifying modules...", value = 0.2)
      ranking <- wgcna_rank_modules(w)

      ## Step 3: Load board rules
      board_rules <- paste(readLines(
        file.path(BOARD_PROMPTS_DIR, "wgcna_report_rules.md"),
        warn = FALSE
      ), collapse = "\n")

      ## Step 4: Species context (graceful: empty string if unknown)
      organism <- pgx$organism %||% NULL
      species_text <- omicsai::omicsai_species_prompt(organism)

      ## Step 5: Assemble user message (board owns all domain logic)
      user_message <- omicsai::collapse_lines(
        board_rules,
        species_text,
        "---",
        "## Module Signal Classification",
        omicsai::omicsai_format_ranking(ranking),
        "---",
        "## Input Data",
        tables$text,
        sep = "\n\n"
      )

      ## Cache the full prompt for instant toggle
      sys_prompt <- tryCatch({
        fp <- omicsai::omicsai_prompt_path("report_format.md")
        txt <- paste(readLines(fp, warn = FALSE), collapse = "\n")
        omicsai::omicsai_substitute_template(txt, list(max_words = "1500"))
      }, error = function(e) "(report_format.md not found)")

      report_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", sys_prompt,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", user_message
      ))

      ## Step 6: Generate report via single LLM call
      progress$set(message = "Generating report...", value = 0.3)

      result <- omicsai::omicsai_gen_report(
        template = user_message,
        model = model
      )

      ## Step 7: Append deterministic methods section
      methods <- wgcna_build_methods(w, pgx)
      full_report <- paste(result$text, methods, sep = "\n\n")

      progress$set(message = "Done!", value = 1)
      full_report
    }, ignoreNULL = FALSE)

    # ---- Unified text content (instant toggle between result and prompt) ----

    ai_text <- shiny::reactive({
      mode <- controls$mode()
      show <- isTRUE(controls$show_prompt())

      if (mode == "summary") {
        if (show) summary_prompt_cache() else ai_summary()
      } else {
        if (show) report_prompt_cache() else ai_report()
      }
    })

    list(
      text = ai_text,
      report_text = ai_report
    )
  })
}
