##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Methylome Profiler — shared helpers for the screen modules.
##
## Everything here reads the loaded pgx and nothing else: pgx$X (beta),
## pgx$genes (450K/EPIC annotation) and pgx$gx.meta (the EWAS already computed
## at upload). No screen needs a contrast except EWAS.

## Shown by every panel on the EWAS screen until a model has been fitted. One
## string in one place, so the whole screen says the same thing and says what
## to press rather than rendering an empty card.
MP_NO_EWAS <- paste(
  "No model has been fitted yet.",
  "Open the settings panel on the right, choose an Outcome,",
  "and press Run EWAS. Everything on this screen follows that fit."
)

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

## pgx$X is beta for every methylomics dataset the pipeline produces. The
## conversion below is a safety net for a matrix that is genuinely on the M
## scale, and the threshold matters: "anything outside [0,1]" is not that test.
## Batch correction and SVD imputation both push a real beta matrix a little
## past its bounds - measured at -0.24 and 1.28 on this demo - and converting
## then destroys it silently, collapsing sd 0.35 to 0.06 inside [0,1] where
## nothing can detect the damage afterwards. M-values live on roughly
## [-20, 20], so anything beyond [-1, 2] is unambiguous and a smaller
## excursion is a bounded beta that spilled; those are clamped, not converted.
## validate() inside an observeEvent is swallowed: Shiny treats it as a
## silent error, the observer just stops, and the button reads as having
## worked - the withProgress flash even looks like success. Every refusal a
## button can hit has to be turned back into something the user sees.
mp_on_click <- function(expr) {
  tryCatch(expr, validation = function(e) {
    msg <- conditionMessage(e)
    ## try(): showNotification needs a reactive domain and errors without one,
    ## which would turn a refusal into a crash when this is called from a test.
    try(shiny::showNotification(
      if (nzchar(msg)) msg else "There is nothing to run yet.",
      type = "warning", duration = 12), silent = TRUE)
    invisible(msg)
  })
}

mp_beta <- function(pgx) {
  mp_require_methylomics(pgx)
  X <- pgx$X
  ## Not range(X, na.rm = TRUE): range.default drops the NAs by SUBSETTING,
  ## which copies the whole matrix - measured at +3.9 GB on an 850k x 200 EPIC
  ## cohort, and every panel in this app enters through here. min and max take
  ## na.rm in C and allocate nothing.
  rng <- c(min(X, na.rm = TRUE), max(X, na.rm = TRUE))
  if (!all(is.finite(rng))) return(X)
  if (rng[2] > 2 || rng[1] < -1) {
    return(playbase.epigenetics::mToBeta(X))
  }
  ## The pipeline no longer produces a spilled beta - batch correction runs on
  ## M and returns through 2^m/(1+2^m), which cannot leave (0,1) - but this
  ## clamp is not therefore dead. Every methylomics pgx written before that
  ## change can still carry values outside [0,1], and SVD2 imputation reaches
  ## them by another route. It is a compatibility guard now, not a live one.
  if (rng[2] > 1 || rng[1] < 0) {
    X[] <- pmin(pmax(X, 0), 1)
  }
  X
}

## The per-sample QC computation moved to playbase.epigenetics::sample_qc(),
## which takes a matrix, a probe annotation and a sample sheet rather than a
## pgx - the upload step must show the same table before any pgx exists. What
## stays here are the pgx-unpacking wrappers the screens and tests call.

## Does the dataset even carry sex chromosomes? GEO matrices are often shipped
## with filterXY already applied, in which case we refuse to predict.
mp_has_xy <- function(pgx, min_probes = 100) {
  playbase.epigenetics::has_sexchr(pgx$X, pgx$genes, min_probes = min_probes)
}

## Recorded sex from the sample sheet, normalised to "F"/"M".
mp_recorded_sex <- function(pgx) {
  playbase.epigenetics::recorded_sex(pgx$samples)
}

## Per-sample verdict on predicted vs recorded sex, aligned to `samples`.
mp_sex_check <- function(pgx, predicted, samples) {
  playbase.epigenetics::sex_check(pgx$samples, predicted, samples)
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

