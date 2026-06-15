##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' WGCNA module AI summary UI
#'
#' Renders an on-demand AI summary card for the selected WGCNA module.
#'
#' @param id Shiny module namespace ID.
#' @param label Optional PlotModule label.
#' @param title Card title.
#' @param info.text Help text.
#' @param caption Caption text.
#' @param height Card height vector c(default, fullscreen).
#' @param width Card width vector c(default, fullscreen).
#'
#' @return Shiny UI tags.
wgcna_module_ai_summary_ui <- function(id,
                                       label = "",
                                       title = "",
                                       info.text = "",
                                       caption = "",
                                       height,
                                       width) {
  AiTextCardUI(
    id = id,
    title = title,
    caption = caption,
    info.text = info.text,
    height = height,
    width = width
  )
}

#' WGCNA module AI summary server
#'
#' Builds current-module evidence from the live WGCNA object and generates a
#' fresh AI summary on demand. Results are not written back to PGX.
#'
#' @param id Shiny module namespace ID.
#' @param wgcna Reactive returning the current WGCNA object.
#' @param pgx Current PGX object.
#' @param r_module Reactive returning the selected module name.
#' @param parent_session Parent Shiny session with AI user options.
#' @param watermark Logical; add watermark in PlotModule.
#'
#' @return Reactive returning the latest omicsai result, or NULL.
wgcna_module_ai_summary_server <- function(id,
                                           wgcna,
                                           pgx,
                                           r_module,
                                           parent_session,
                                           watermark = FALSE) {
  # Build prompt data from the currently selected module and live WGCNA object.
  summary_params <- shiny::reactive({
    res <- wgcna()
    module <- r_module()
    shiny::req(res, module)
    wgcna_module_summary_params(res, module, pgx)
  })

  # Assemble the structured omicsai prompt into system and user-message parts.
  prompt_parts <- shiny::reactive({
    params <- summary_params()
    organism <- tryCatch(pgx$organism, error = function(e) NULL)
    prompt <- omicsai::summary_prompt(
      role = omicsai::frag("system_base"),
      task = NULL,
      species = omicsai::omicsai_species_prompt(organism),
      context = omicsai::frag("wgcna/wgcna_interpretation"),
      data = omicsai::frag("wgcna/wgcna_module_data", params)
    )
    omicsai::build_prompt(prompt)
  })

  # Pass the already-built user message through the generic AI text card.
  template <- shiny::reactive(prompt_parts()$board)

  # Keep model selection tied to app Settings while using the WGCNA system prompt.
  config <- shiny::reactive({
    model <- get_ai_model(parent_session)
    system_prompt <- prompt_parts()$system
    omicsai::omicsai_config(
      model = model,
      system_prompt = system_prompt
    )
  })

  # Render the serial on-demand AI summary card; no PGX write-back happens here.
  AiTextCardServer(
    id = id,
    params_reactive = shiny::reactive(list()),
    template_reactive = template,
    config_reactive = config,
    cache = omicsai::omicsai_cache_init("mem"),
    watermark = watermark
  )
}

#' Build prompt parameters for one WGCNA module
#'
#' @param wgcna WGCNA result object.
#' @param module Character module name such as `"MEblue"`.
#' @param pgx Current PGX object.
#'
#' @return Named list matching `omicsai` prompt `wgcna/wgcna_module_data`.
wgcna_module_summary_params <- function(wgcna, module, pgx) {
  # Normalize source tables before mapping them into the omicsai prompt schema.
  input <- .wgcna_module_input_data(wgcna, module)

  # Keep this list aligned with inst/prompts/wgcna/wgcna_module_data.md.
  list(
    ME_color = module,
    n_genes = length(input$module_genes),
    tier = .wgcna_module_tier(input$module_traits, module, input$gse),
    eigengene_profile_qualitative = .wgcna_module_eigengene_profile(wgcna, module, pgx),
    top_pos_trait = .wgcna_module_trait(input$module_traits, module, "positive"),
    top_pos_verbal = .wgcna_module_trait(input$module_traits, module, "positive", verbal = TRUE),
    top_neg_trait = .wgcna_module_trait(input$module_traits, module, "negative"),
    top_neg_verbal = .wgcna_module_trait(input$module_traits, module, "negative", verbal = TRUE),
    n_sig_terms = sum(input$gse$q.value < 0.05, na.rm = TRUE),
    n_total_terms = nrow(input$gse),
    enrichment_themes_table = .wgcna_module_enrichment_table(input$gse),
    n_hub = 50L,
    hub_genes_table = .wgcna_module_hub_table(
      wgcna,
      module,
      pgx,
      input$module_traits,
      input$top
    ),
    gene_families_summary = "",
    enrichment_overlaps = .wgcna_module_enrichment_overlaps(input$gse)
  )
}

# Normalize WGCNA module aliases and table shapes for prompt assembly.
# This is necessary because precomputed and live WGCNA objects expose different shapes.
.wgcna_module_input_data <- function(wgcna, module) {
  module_key <- module
  module_alt <- sub("^ME", "", module)
  module_genes <- wgcna$me.genes[[module_key]]
  if (is.null(module_genes) && module_alt != module_key) {
    module_genes <- wgcna$me.genes[[module_alt]]
  }
  if (is.null(module_genes)) {
    module_genes <- character(0)
  }
  module_traits <- tryCatch(playbase:::wgcna.get_modTraits(wgcna), error = function(e) NULL)
  top <- tryCatch(
    playbase::wgcna.getTopGenesAndSets(
      wgcna,
      annot = wgcna$annot,
      ntop = 40,
      level = "gene",
      rename = "gene_title"
    ),
    error = function(e) list(pheno = list(), genes = list(), sets = list())
  )
  top_genes <- top$genes[[module_key]]
  top_sets <- top$sets[[module_key]]
  if (is.null(top_genes) && module_alt != module_key) {
    top_genes <- top$genes[[module_alt]]
  }
  if (is.null(top_sets) && module_alt != module_key) {
    top_sets <- top$sets[[module_alt]]
  }
  top$genes[[module_key]] <- if (is.null(top_genes)) character(0) else top_genes
  top$sets[[module_key]] <- if (is.null(top_sets)) character(0) else top_sets

  if (is.data.frame(wgcna$gse)) {
    gse <- wgcna$gse[wgcna$gse$module %in% c(module_key, module_alt), , drop = FALSE]
  } else {
    gse <- wgcna$gse[[module_key]]
    if (is.null(gse) && module_alt != module_key) {
      gse <- wgcna$gse[[module_alt]]
    }
  }
  if (!is.data.frame(gse)) {
    gse <- data.frame()
  }
  if (nrow(gse) > 0 && !"score" %in% colnames(gse) &&
      all(c("odd.ratio", "p.value") %in% colnames(gse))) {
    gse$odd.ratio[is.infinite(gse$odd.ratio)] <- 99
    gse$score <- gse$odd.ratio * -log10(gse$p.value)
  }

  list(
    module_genes = module_genes,
    module_traits = module_traits,
    top = top,
    gse = gse
  )
}

# Classify module signal from trait correlation and significant enrichment count.
# This helps us give the prompt a compact strength label instead of raw thresholds.
.wgcna_module_tier <- function(module_traits, module, gse) {
  max_cor <- 0
  if (!is.null(module_traits) && module %in% rownames(module_traits)) {
    max_cor <- max(abs(module_traits[module, ]), na.rm = TRUE)
  }
  n_sig <- if ("q.value" %in% colnames(gse)) sum(gse$q.value < 0.05, na.rm = TRUE) else 0L
  if (max_cor >= 0.7 || n_sig >= 10L) return("strong signal")
  if (max_cor >= 0.5 || n_sig >= 3L) return("moderate signal")
  "weak signal"
}

# Summarize eigengene group means for the selected module.
# This helps us expose module directionality without passing full sample matrices.
.wgcna_module_eigengene_profile <- function(wgcna, module, pgx) {
  datME <- if (!is.null(wgcna$datME)) wgcna$datME else wgcna$net$MEs
  groups <- pgx$samples$group
  if (is.null(datME) || !module %in% colnames(datME) || is.null(groups)) {
    return("not available")
  }
  group_means <- sort(tapply(datME[, module], groups, mean), decreasing = TRUE)
  paste(
    sprintf("%s: %+.2f", names(group_means), group_means),
    collapse = "; "
  )
}

# Return the strongest positive or negative trait above the reporting threshold.
# Separating the code here keeps trait picking consistent across summary fields.
.wgcna_module_trait <- function(module_traits,
                                module,
                                direction = c("positive", "negative"),
                                verbal = FALSE) {
  direction <- match.arg(direction)
  if (is.null(module_traits) || !module %in% rownames(module_traits)) {
    return("")
  }
  cors <- module_traits[module, ]
  cors <- cors[!is.na(cors)]
  cors <- if (identical(direction, "positive")) cors[cors >= 0.5] else cors[cors <= -0.5]
  if (length(cors) == 0) {
    return("")
  }
  value <- if (identical(direction, "positive")) max(cors) else min(cors)
  trait <- names(cors)[which(cors == value)[1]]
  if (isTRUE(verbal)) sprintf("r = %+.2f", value) else trait
}

# Format the top significant enrichment terms for the omicsai data fragment.
# This helps us keep WGCNA-specific prompt wording out of generic formatters.
.wgcna_module_enrichment_table <- function(gse) {
  needed <- c("geneset", "score", "q.value")
  if (!is.data.frame(gse) || nrow(gse) == 0 || !all(needed %in% colnames(gse))) {
    return("No enrichment terms available.")
  }
  sig <- gse[!is.na(gse$q.value) & gse$q.value < 0.05, , drop = FALSE]
  if (nrow(sig) == 0) {
    return("No enrichment terms passed q < 0.05.")
  }
  sig <- sig[order(sig$q.value, -sig$score), , drop = FALSE]
  df <- data.frame(
    theme = head(sig$geneset, 20),
    score = head(sig$score, 20),
    q = head(sig$q.value, 20),
    stringsAsFactors = FALSE
  )
  paste(omicsai::omicsai_format_mdtable(df, formatters = list(
    score = function(x) omicsai::omicsai_format_num(x, 2),
    q = omicsai::omicsai_format_pvalue
  )), collapse = "\n")
}

# Format hub genes using gene stats when a correlated trait is available.
# This is necessary because hub ranking depends on the selected module and trait.
.wgcna_module_hub_table <- function(wgcna, module, pgx, module_traits, top) {
  trait <- .wgcna_module_trait(module_traits, module, "positive")
  if (!nzchar(trait)) {
    trait <- .wgcna_module_trait(module_traits, module, "negative")
  }
  genes <- top$genes[[module]]
  if (length(genes) == 0) {
    module_alt <- sub("^ME", "", module)
    genes <- wgcna$me.genes[[module]]
    if (is.null(genes) && module_alt != module) {
      genes <- wgcna$me.genes[[module_alt]]
    }
  }
  if (is.null(genes)) {
    genes <- character(0)
  }
  stats <- if (nzchar(trait)) {
    tryCatch(
      playbase::wgcna.getGeneStats(wgcna, module = module, trait = trait, plot = FALSE),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  if (is.data.frame(stats) && nrow(stats) > 0 && "moduleMembership" %in% colnames(stats)) {
    stats <- stats[order(-abs(stats$moduleMembership)), , drop = FALSE]
    genes <- head(stats$feature, 50)
  } else {
    genes <- head(genes, 50)
  }
  symbols <- .wgcna_module_symbols(genes, pgx$genes)
  funcs <- .wgcna_module_functions(genes, pgx$genes)
  paste(omicsai::omicsai_format_mdtable(data.frame(
    gene = symbols,
    known_function = funcs,
    stringsAsFactors = FALSE
  )), collapse = "\n")
}

# Resolve feature IDs to gene symbols.
# This helps us let the AI summary mention readable annotated gene names.
.wgcna_module_symbols <- function(features, annot) {
  symbols <- tryCatch(
    playbase::probe2symbol(features, annot, "symbol"),
    error = function(e) features
  )
  symbols[is.na(symbols) | symbols == ""] <- features[is.na(symbols) | symbols == ""]
  symbols
}

# Resolve feature IDs to short gene descriptions.
# This helps us include compact biological context for hub genes.
.wgcna_module_functions <- function(features, annot) {
  out <- rep("function not annotated", length(features))
  if (is.null(annot) || is.null(rownames(annot))) {
    return(out)
  }
  col <- intersect(c("gene_title", "gene_name", "description"), colnames(annot))
  if (length(col) == 0) {
    return(out)
  }
  idx <- match(features, rownames(annot))
  ok <- !is.na(idx)
  out[ok] <- as.character(annot[idx[ok], col[1]])
  out[is.na(out) | out == ""] <- "function not annotated"
  substr(out, 1, 90)
}

# Format overlap gene lists for top enrichment terms.
# This helps us keep overlap wording tailored to the WGCNA prompt.
.wgcna_module_enrichment_overlaps <- function(gse) {
  if (!is.data.frame(gse) || nrow(gse) == 0 || !"genes" %in% colnames(gse)) {
    return("")
  }
  sig <- if ("q.value" %in% colnames(gse)) {
    gse[!is.na(gse$q.value) & gse$q.value < 0.05, , drop = FALSE]
  } else {
    gse
  }
  if (nrow(sig) == 0) {
    return("")
  }
  sig <- head(sig, 5)
  paste(
    sprintf("- %s: %s", sig$geneset, sig$genes),
    collapse = "\n"
  )
}
