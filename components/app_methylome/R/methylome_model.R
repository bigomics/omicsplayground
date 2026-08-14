##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## The differential-methylation model, fitted in-app so covariates can be
## chosen at analysis time.
##
## Why this exists: the result stored on the pgx at upload is an unadjusted
## two-group comparison. In bulk tissue that is usually a test of cell
## composition rather than of the phenotype - most published blood CpG-age
## associations lose significance once composition is adjusted for (Jaffe &
## Irizarry 2014, Genome Biol 15:R31). Anything that shifts the neutrophil to
## lymphocyte ratio (infection, obesity, smoking, age) produces a large,
## highly significant and meaningless hit list. So the board refits limma
## itself with whatever covariates the user selects.

## --------------------------------------------------------------- masking ---

## Probes whose signal is genotype rather than methylation: a common SNP in
## the probe body, at the target CpG, or at the single-base extension site.
## Roughly 10% of the 450K array. Read straight from the manifest we already
## load, so this costs nothing.
mp_snp_probes <- function(pgx, maf = 0.05) {
  ann <- mp_manifest(pgx)
  if (is.null(ann)) return(character(0))
  cols <- intersect(c("Probe_maf", "CpG_maf", "SBE_maf"), colnames(ann))
  if (!length(cols)) return(character(0))
  hit <- Reduce(`|`, lapply(cols, function(k) {
    v <- suppressWarnings(as.numeric(ann[[k]]))
    !is.na(v) & v > maf
  }))
  out <- rownames(ann)[hit]
  ## ~10% of the 450K carries a common SNP; zero means the manifest did not
  ## load and the user would otherwise get an unmasked run believing it masked.
  if (!length(out)) warning("[methylome] SNP mask matched no probes - check the manifest")
  out
}

## Cross-reactive / multi-mapping probes, from the bundled published lists.
## The board is sourced, not installed, so search the plausible working
## directories: the app runs from components/app/R, tests from
## app_methylome/test.
MP_XREACTIVE_CACHE <- new.env(parent = emptyenv())

mp_xreactive_probes <- function() {
  if (!is.null(MP_XREACTIVE_CACHE$p)) return(MP_XREACTIVE_CACHE$p)
  rel <- file.path("masking", "cross_reactive_probes.csv")
  cand <- c(
    system.file(rel, package = "app_methylome"),
    file.path("../../app_methylome/inst", rel),
    file.path("app_methylome/inst", rel),
    file.path("../inst", rel),
    file.path("components/app_methylome/inst", rel)
  )
  f <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(f)) {
    warning("[methylome] cross-reactive probe list not found; that mask is a no-op")
    MP_XREACTIVE_CACHE$p <- character(0)
    return(character(0))
  }
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  MP_XREACTIVE_CACHE$p <- unique(as.character(d$probe))
  MP_XREACTIVE_CACHE$p
}

## getAnnotation() dispatches through the annotation object, which must be
## ATTACHED, not merely loaded - with only requireNamespace() it fails with
## "no item called package:... on the search list" and the mask silently
## becomes a no-op. Same trap as agep() and data("age_coefficients").
mp_manifest <- function(pgx) {
  pkg <- if (!is.null(pgx$meth_type) && grepl("EPIC", pgx$meth_type, ignore.case = TRUE)) {
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
  } else {
    "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    warning("[methylome] ", pkg, " is not installed; SNP masking disabled")
    return(NULL)
  }
  if (!paste0("package:", pkg) %in% search()) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  tryCatch(as.data.frame(minfi::getAnnotation(getExportedValue(pkg, pkg))),
           error = function(e) {
             warning("[methylome] could not read the ", pkg, " manifest: ",
                     conditionMessage(e))
             NULL
           })
}

## Which probes survive the selected masks.
mp_masked_probes <- function(pgx, mask = c("snp", "xreactive")) {
  out <- character(0)
  if ("snp" %in% mask) out <- c(out, mp_snp_probes(pgx))
  if ("xreactive" %in% mask) out <- c(out, mp_xreactive_probes())
  unique(out)
}

## ----------------------------------------------------------------- model ---

## Sample columns usable as covariates: not the contrast itself, not free text.
mp_model_vars <- function(pgx) {
  s <- pgx$samples
  if (is.null(s)) return(character(0))
  keep <- vapply(colnames(s), function(k) {
    v <- s[[k]]
    if (is.numeric(v)) return(stats::sd(v, na.rm = TRUE) > 0)
    u <- unique(as.character(v[!is.na(v) & v != ""]))
    length(u) >= 2 && length(u) <= 12
  }, logical(1))
  ## dot-prefixed columns are app-derived metrics, still valid covariates
  sort(colnames(s)[keep])
}

#' Fit the differential-methylation model.
#'
#' @return list(table, formula, n, groups, dropped_covars, masked, ok)
mp_fit_ewas <- function(pgx, contrast, covars = character(0),
                        cellfracs = NULL, mask = character(0)) {
  mp_require_methylomics(pgx)
  shiny::validate(shiny::need(!is.null(contrast), "No contrast selected."))

  beta <- mp_beta(pgx)
  grp <- mp_contrast_groups(pgx, contrast)
  shiny::validate(shiny::need(!is.null(grp), "This contrast has no group labels."))

  ss <- intersect(names(grp), colnames(beta))
  grp <- droplevels(factor(grp[ss]))
  shiny::validate(shiny::need(nlevels(grp) == 2,
    "This model handles two-group contrasts."))

  ## Covariate frame, aligned to the contrast's samples.
  dat <- data.frame(.group = grp, row.names = ss, stringsAsFactors = FALSE)
  used <- character(0); dropped <- character(0)
  for (k in covars) {
    v <- pgx$samples[ss, k]
    if (is.numeric(v)) {
      if (stats::sd(v, na.rm = TRUE) > 0) { dat[[k]] <- v; used <- c(used, k) }
      else dropped <- c(dropped, paste0(k, " (constant)"))
    } else {
      f <- droplevels(factor(as.character(v)))
      if (nlevels(f) >= 2) { dat[[k]] <- f; used <- c(used, k) }
      else dropped <- c(dropped, paste0(k, " (single level)"))
    }
  }
  if (!is.null(cellfracs)) {
    cf <- cellfracs[ss, , drop = FALSE]
    ## One fraction is redundant with the rest (they sum to 1); drop the
    ## largest so the design stays full rank.
    cf <- cf[, order(-colMeans(cf, na.rm = TRUE)), drop = FALSE][, -1, drop = FALSE]
    for (k in colnames(cf)) { dat[[k]] <- as.numeric(cf[, k]); used <- c(used, k) }
  }

  ## Drop samples missing any model variable, then re-check the contrast.
  keep <- stats::complete.cases(dat)
  dat <- dat[keep, , drop = FALSE]
  shiny::validate(shiny::need(nrow(dat) >= 6, "Too few complete samples for this model."))
  dat$.group <- droplevels(dat$.group)
  shiny::validate(shiny::need(nlevels(dat$.group) == 2,
    "One group has no samples left once incomplete covariates are dropped."))

  design <- stats::model.matrix(stats::as.formula(
    paste("~ 0 + .group", if (length(used)) paste("+", paste(used, collapse = " + ")) else "")
  ), data = dat)
  colnames(design) <- make.names(colnames(design))

  ## Refuse a rank-deficient design rather than returning silent nonsense.
  ne <- limma::nonEstimable(design)
  shiny::validate(shiny::need(is.null(ne),
    paste0("These terms are collinear with the contrast and cannot be fitted: ",
           paste(ne, collapse = ", "), ". Remove one of them.")))
  shiny::validate(shiny::need(nrow(dat) - ncol(design) >= 2,
    "Not enough residual degrees of freedom - drop a covariate."))

  probes <- rownames(beta)
  masked <- intersect(mp_masked_probes(pgx, mask), probes)
  probes <- setdiff(probes, masked)
  shiny::validate(shiny::need(length(probes) > 100, "Masking removed nearly every probe."))

  ## Fit on M-values, report on beta. M is the right scale for testing
  ## (Du 2010) but is not interpretable as a change in methylation.
  B <- beta[probes, rownames(dat), drop = FALSE]
  M <- playbase::betaToM(B)

  lv <- levels(dat$.group)
  g1 <- make.names(paste0(".group", lv[1]))
  g2 <- make.names(paste0(".group", lv[2]))
  fit <- limma::lmFit(M, design)
  cm <- limma::makeContrasts(contrasts = paste0(g2, "-", g1), levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm), trend = TRUE)
  tt <- limma::topTable(fit2, number = Inf, sort.by = "none")

  dbeta <- rowMeans(B[, dat$.group == lv[2], drop = FALSE], na.rm = TRUE) -
           rowMeans(B[, dat$.group == lv[1], drop = FALSE], na.rm = TRUE)

  list(
    table = data.frame(
      probe = rownames(tt),
      logFC_M = tt$logFC,
      dbeta = round(dbeta[rownames(tt)], 4),
      p = tt$P.Value, q = tt$adj.P.Val,
      stringsAsFactors = FALSE
    ),
    formula = paste0("~ ", lv[2], " vs ", lv[1],
                     if (length(used)) paste0(" + ", paste(used, collapse = " + ")) else
                       "   (unadjusted)"),
    n = nrow(dat), groups = table(dat$.group),
    dropped_covars = dropped, masked = length(masked),
    covars = used
  )
}
