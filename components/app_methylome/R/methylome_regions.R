##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Regions and pathways: a real DMR caller, and gene-set testing corrected for
## the probes-per-gene bias.
##
## dmrff rather than DMRcate: dmrff runs on EWAS summary statistics, which is
## exactly what the model step already produces, and its entire dependency
## list is parallel + ggplot2. DMRcate would drag in bsseq, Gviz, biomaRt,
## missMethyl and an ExperimentHub package that downloads at runtime. A 2024
## BMC Bioinformatics comparison also put dmrff at the highest power with a
## controlled false-positive rate.

## ------------------------------------------------------------------- DMRs --

mp_can_dmrff <- function() requireNamespace("dmrff", quietly = TRUE)

#' Call DMRs from the fitted model.
#'
#' dmrff needs the effect estimate, its standard error and the p-value per
#' probe, plus chromosome and position. limma does not return the standard
#' error directly, but se = logFC / t and t = logFC / se, so it comes back
#' out of the moderated t statistic.
mp_call_dmrs <- function(pgx, res, maxgap = 500) {
  shiny::validate(shiny::need(mp_can_dmrff(),
    "The dmrff package is not installed, so region calling is unavailable."))
  d <- res$data
  fit <- res$meta
  shiny::validate(shiny::need(!is.null(fit$t),
    "This model did not return moderated t statistics."))

  keep <- !is.na(d$chr) & !is.na(d$pos) & is.finite(d$pos) &
    is.finite(d$p) & is.finite(fit$t) & fit$t != 0
  d <- d[keep, , drop = FALSE]
  tt <- fit$t[keep]
  shiny::validate(shiny::need(nrow(d) > 100, "Too few positioned probes to call regions."))

  est <- d$logFC_M
  se <- abs(est / tt)
  se[!is.finite(se) | se <= 0] <- NA

  beta <- mp_beta(pgx)
  grp <- mp_contrast_groups(pgx, res$contrast)
  ss <- intersect(names(grp), colnames(beta))
  M <- playbase::betaToM(beta[d$probe, ss, drop = FALSE])

  ok <- !is.na(se)
  out <- tryCatch(
    dmrff::dmrff(estimate = est[ok], se = se[ok], p.value = d$p[ok],
                 methylation = M[ok, , drop = FALSE],
                 chr = paste0("chr", d$chr[ok]), pos = d$pos[ok],
                 maxgap = maxgap, verbose = FALSE),
    error = function(e) {
      warning("[methylome] dmrff failed: ", conditionMessage(e)); NULL
    })
  shiny::validate(shiny::need(!is.null(out), "Region calling failed on this dataset."))
  ## Single-CpG "regions" are just DMPs and are already in the hits table.
  out <- out[!is.na(out$n) & out$n > 1, , drop = FALSE]
  out
}

## Genes overlapping each region, from the annotation already loaded.
mp_dmr_genes <- function(pgx, dmrs, d) {
  if (!nrow(dmrs)) return(character(0))
  vapply(seq_len(nrow(dmrs)), function(i) {
    hit <- d$chr == sub("^chr", "", dmrs$chr[i]) &
      d$pos >= dmrs$start[i] & d$pos <= dmrs$end[i]
    g <- unique(d$gene[hit])
    g <- g[!is.na(g) & g != ""]
    if (!length(g)) "-" else paste(utils::head(g, 4), collapse = ", ")
  }, character(1))
}

methylome_table_dmr_ui <- function(id, title, info.text, caption, height, width) {
  ns <- shiny::NS(id)
  TableModuleUI(ns("tblmod"), title = title, info.text = info.text,
                caption = caption, height = height, width = width, label = "a")
}

methylome_table_dmr_server <- function(id, pgx, r.regions, scrollY = "26vh") {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      dd <- r.regions()
      shiny::validate(shiny::need(!is.null(dd) && nrow(dd$dmrs) > 0,
        "No differentially methylated regions found. Press Call regions in the settings panel, or relax the model."))
      x <- dd$dmrs
      data.frame(
        Region = sprintf("%s:%s-%s", x$chr, format(x$start, big.mark = ","),
                         format(x$end, big.mark = ",")),
        Genes = dd$genes,
        `CpGs` = x$n,
        `Width (bp)` = x$end - x$start,
        `Effect (M)` = signif(x$estimate, 3),
        `P value` = signif(x$p.value, 3),
        `Adjusted p` = signif(x$p.adjust, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })
    render <- function(sy) {
      dt <- table_data(); shiny::req(dt)
      DT::datatable(dt, class = "compact hover", rownames = FALSE, selection = "none",
        extensions = c("Buttons", "Scroller"), plugins = "scrollResize",
        options = list(dom = "lfrtip", scroller = TRUE, scrollX = TRUE,
                       scrollY = sy, scrollResize = TRUE, deferRender = TRUE,
                       order = list(list(5, "asc")))) |>
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }
    TableModuleServer("tblmod", func = function() render(scrollY),
                      func2 = function() render("60vh"),
                      csvFunc = table_data, selector = "none")
  })
}

## ------------------------------------------------------------- enrichment --

mp_can_gometh <- function() requireNamespace("missMethyl", quietly = TRUE)

#' Gene-set testing on the hit list, corrected for probes per gene.
#'
#' Genes carry very different numbers of probes, so a naive hypergeometric
#' test returns the same probe-dense developmental gene sets whatever the
#' biology (Geeleher 2013). gometh reweights with Wallenius' noncentral
#' hypergeometric and also handles CpGs mapping to several genes.
mp_run_gometh <- function(sig_cpg, all_cpg, collection = "GO", array_type = "450K") {
  shiny::validate(shiny::need(mp_can_gometh(),
    "The missMethyl package is not installed, so bias-corrected enrichment is unavailable."))
  shiny::validate(shiny::need(length(sig_cpg) >= 10,
    "Too few hits at this threshold for a meaningful gene-set test."))
  out <- tryCatch(
    missMethyl::gometh(sig.cpg = sig_cpg, all.cpg = all_cpg,
                       collection = collection, array.type = array_type,
                       prior.prob = TRUE),
    error = function(e) {
      warning("[methylome] gometh failed: ", conditionMessage(e)); NULL
    })
  shiny::validate(shiny::need(!is.null(out), "Gene-set testing failed on this dataset."))
  out$TERM_ID <- rownames(out)
  out[order(out$P.DE), , drop = FALSE]
}

methylome_table_enrich_ui <- function(id, title, info.text, caption, height, width) {
  ns <- shiny::NS(id)
  TableModuleUI(ns("tblmod"), title = title, info.text = info.text,
                caption = caption, height = height, width = width, label = "b")
}

methylome_table_enrich_server <- function(id, r.enrich, scrollY = "26vh") {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      e <- r.enrich()
      shiny::validate(shiny::need(!is.null(e) && nrow(e) > 0,
        "No gene-set result. Press Test gene sets in the settings panel."))
      keep <- utils::head(e, 300)
      data.frame(
        Term = if ("TERM" %in% colnames(keep)) keep$TERM else keep$TERM_ID,
        Ontology = if ("ONTOLOGY" %in% colnames(keep)) keep$ONTOLOGY else "-",
        `Genes in set` = keep$N,
        `Hits` = keep$DE,
        `P value` = signif(keep$P.DE, 3),
        FDR = signif(keep$FDR, 3),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })
    render <- function(sy) {
      dt <- table_data(); shiny::req(dt)
      DT::datatable(dt, class = "compact hover", rownames = FALSE, selection = "none",
        extensions = c("Buttons", "Scroller"), plugins = "scrollResize",
        options = list(dom = "lfrtip", scroller = TRUE, scrollX = TRUE,
                       scrollY = sy, scrollResize = TRUE, deferRender = TRUE)) |>
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }
    TableModuleServer("tblmod", func = function() render(scrollY),
                      func2 = function() render("60vh"),
                      csvFunc = table_data, selector = "none")
  })
}
