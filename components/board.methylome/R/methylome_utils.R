##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Methylome Profiler — shared helpers for the screen modules.
##
## Everything here reads the loaded pgx and nothing else: pgx$X (beta),
## pgx$genes (450K/EPIC annotation) and pgx$gx.meta (the EWAS already computed
## at upload). No screen needs a contrast except EWAS.

MP_PAL <- list(
  hypo = "#4575b4", hyper = "#d73027", grey = "#8697a6",
  ok = "#18753a", bad = "#a8202c", warn = "#8a5a06"
)

## ---------------------------------------------------------------- helpers ---

## Every screen calls this first. Without it the app renders confident
## nonsense on an RNA-seq or proteomics dataset - beta "densities" of log
## counts, a DNAm age for samples that have no CpGs.
mp_require_methylomics <- function(pgx) {
  shiny::validate(shiny::need(
    !is.null(pgx) && !is.null(pgx$datatype) &&
      tolower(pgx$datatype) == "methylomics",
    paste0(
      "Methylome Profiler needs a methylation dataset. ",
      "The loaded dataset is '",
      if (is.null(pgx$datatype)) "unknown" else pgx$datatype,
      "'. Load a methylomics dataset from the Library."
    )
  ))
  invisible(TRUE)
}

mp_beta <- function(pgx) {
  mp_require_methylomics(pgx)
  X <- pgx$X
  if (max(X, na.rm = TRUE) > 1 || min(X, na.rm = TRUE) < 0) X <- playbase::mToBeta(X)
  X
}

## Fraction of probes outside the intermediate band. A healthy methylome is
## bimodal at ~0.2/0.8, so this sits near 0.85-0.95; a flat sample drops.
mp_bimodality <- function(X) {
  round(1 - colMeans(X > 0.3 & X < 0.7, na.rm = TRUE), 3)
}

## Imprinted DMRs sit near 0.5. Deviation means normalisation went wrong.
mp_imprint_drift <- function(X) {
  imp <- tryCatch({
    e <- new.env(); utils::data("iDMR", package = "wateRmelon", envir = e)
    get("iDMR", envir = e)
  }, error = function(e) character(0))
  ii <- intersect(rownames(X), imp)
  if (length(ii) < 10) return(rep(NA_real_, ncol(X)))
  round(colMeans(abs(X[ii, , drop = FALSE] - 0.5), na.rm = TRUE), 3)
}

## method: "all" fits the five wateRmelon clocks, which is what the Epigenetic
## age fallback needs. The ledger shows one number and used to pay for five, so
## it asks for one - the dominant cost on that screen.
mp_clocks <- function(X, method = "all") {
  ## agep() calls data("age_coefficients") with no package= argument, so the
  ## coefficients only resolve when wateRmelon is attached. Calling it
  ## namespace-qualified silently yields NA ages.
  if (!"package:wateRmelon" %in% search()) {
    suppressPackageStartupMessages(library(wateRmelon))
  }
  A <- tryCatch(wateRmelon::agep(X, method = method), error = function(e) NULL)
  if (is.null(A)) return(NULL)
  A <- as.data.frame(A)
  keep <- grep("\\.age$", colnames(A), value = TRUE)
  ## A single-method call returns one unsuffixed column, so fall back to
  ## whatever came back rather than returning an empty frame.
  if (!length(keep)) {
    colnames(A) <- method
    return(round(A, 1))
  }
  out <- A[, keep, drop = FALSE]
  colnames(out) <- sub("\\..*", "", keep)
  round(out, 1)
}

mp_sex <- function(X) {
  chr <- NULL
  s <- tryCatch(wateRmelon::estimateSex(X, do_plot = FALSE), error = function(e) NULL)
  if (is.null(s) || !"predicted_sex" %in% names(s)) return(NULL)
  s
}

## Does the dataset even carry sex chromosomes? GEO matrices are often shipped
## with filterXY already applied, in which case we refuse to predict.
mp_has_xy <- function(pgx) {
  cc <- grep("^chr$|chrom", tolower(colnames(pgx$genes)))
  if (!length(cc)) return(FALSE)
  v <- as.character(pgx$genes[[cc[1]]])
  any(grepl("^(chr)?[XY]$", v, ignore.case = TRUE))
}

mp_context_col <- function(pgx, what = c("island", "gene")) {
  what <- match.arg(what)
  cn <- colnames(pgx$genes)
  hit <- if (what == "island") grep("Relation_to_Island", cn, ignore.case = TRUE)
         else grep("genomic_location|RefGene_Group", cn, ignore.case = TRUE)
  if (!length(hit)) return(NULL)
  as.character(pgx$genes[[cn[hit[1]]]])
}

mp_card <- function(title, ..., right = NULL) {
  shiny::div(
    class = "mp-card",
    shiny::div(class = "mp-card-head", title, shiny::span(class = "mp-sp"), right),
    shiny::div(class = "mp-card-body", ...)
  )
}

mp_tile <- function(k, v, d) {
  shiny::div(class = "mp-tile",
    shiny::div(class = "mp-k", k), shiny::div(class = "mp-v", v), shiny::div(class = "mp-d", d))
}

