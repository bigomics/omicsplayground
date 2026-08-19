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
##
## Anchored on OPG rather than on the working directory: the app runs with its
## cwd in components/app/R, and a candidate list of relative paths alone once
## resolved to nothing there, leaving this mask a silent no-op - 53,498 probes
## the checkbox claimed to exclude were being tested, while the summary line
## still reported the SNP mask's count and so looked like both had applied.
## pgx.system.file() is not usable here: it only ever resolves
## components/board.<name>, and this is an app, not a board.
##
## The relative fallbacks are for the headless test harness, which runs from
## the test/ directory without the app's globals.
MP_XREACTIVE_CACHE <- new.env(parent = emptyenv())

mp_xreactive_probes <- function() {
  if (!is.null(MP_XREACTIVE_CACHE$p)) return(MP_XREACTIVE_CACHE$p)
  rel <- file.path("masking", "cross_reactive_probes.csv")
  cand <- character(0)
  if (exists("OPG")) cand <- file.path(OPG, "components", "app_methylome", "inst", rel)
  cand <- c(cand, file.path("../inst", rel),
            file.path("components/app_methylome/inst", rel))
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
## Array type, inferred from the probe IDs rather than from metadata.
##
## playbase takes meth_type as an argument to pgx.createPGX and uses it when
## annotating, but does not store it on the pgx - so pgx$meth_type is NULL and
## anything keyed on it silently falls through to 450K. Inferring from the
## probes themselves is both more robust and works on subsets, where a probe
## count would tell you nothing: EPIC keeps roughly 453k of the 450K probes
## and adds ~413k of its own, so any real EPIC dataset carries a large
## fraction of probes that do not exist on the 450K at all.
MP_ARRAY_CACHE <- new.env(parent = emptyenv())

MP_ANNO_PKG <- c(`450K` = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                 EPIC = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

## Attach-and-read a manifest. getAnnotation() dispatches through the
## annotation object, which must be ATTACHED, not merely loaded - with only
## requireNamespace() it fails with "no item called package:... on the search
## list" and callers silently get nothing.
mp_read_manifest <- function(pkg) {
  key <- paste0("anno_", pkg)
  if (!is.null(MP_ARRAY_CACHE[[key]])) return(MP_ARRAY_CACHE[[key]])
  if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)
  if (!paste0("package:", pkg) %in% search()) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  ann <- tryCatch(as.data.frame(minfi::getAnnotation(getExportedValue(pkg, pkg))),
                  error = function(e) {
                    warning("[methylome] could not read ", pkg, ": ",
                            conditionMessage(e)); NULL
                  })
  MP_ARRAY_CACHE[[key]] <- ann
  ann
}

## EPICv2 probe ids carry a replicate suffix - cg00000029_TC21 - so they match
## no entry in either manifest we ship. Left undetected the EPIC-only fraction
## comes out at zero and the array is reported as 450K, which is the worst
## outcome: annotation, SNP masking and the gometh background all degrade with
## nothing said. Detect it by the suffix and refuse.
MP_EPICV2_RE <- "^cg[0-9]+_[A-Z]{2}[0-9]{2}$"

#' Infer the array from the probes present. Returns "450K", "EPIC" or "EPICv2".
mp_array_type <- function(pgx, min_frac = 0.01) {
  ids <- rownames(pgx$X)
  key <- paste0("type_", length(ids), "_", digest_ids(ids))
  if (!is.null(MP_ARRAY_CACHE[[key]])) return(MP_ARRAY_CACHE[[key]])

  if (mean(grepl(MP_EPICV2_RE, ids)) > min_frac) {
    MP_ARRAY_CACHE[[key]] <- "EPICv2"
    return("EPICv2")
  }

  epic <- mp_read_manifest(MP_ANNO_PKG[["EPIC"]])
  k450 <- mp_read_manifest(MP_ANNO_PKG[["450K"]])
  type <- "450K"
  if (is.null(epic)) {
    warning("[methylome] the EPIC manifest is not installed; assuming 450K")
  } else if (is.null(k450)) {
    type <- "EPIC"
  } else {
    epic_only <- setdiff(rownames(epic), rownames(k450))
    frac <- mean(ids %in% epic_only)
    if (frac > min_frac) type <- "EPIC"
  }
  MP_ARRAY_CACHE[[key]] <- type
  type
}

## Cheap stable key for a probe set - avoids re-inferring on every call
## without pulling in a hashing dependency.
digest_ids <- function(ids) {
  n <- length(ids)
  paste0(ids[1], "_", ids[max(1, n %/% 2)], "_", ids[n])
}

MP_EPICV2_MSG <- paste(
  "This looks like an EPIC v2 array. The EPICv2 annotation package",
  "(IlluminaHumanMethylationEPICv2anno.20a1.hg38) is not installed, so probe",
  "annotation, SNP masking and gene-set testing are not available for it -",
  "and its probe ids carry replicate suffixes that match neither the 450K nor",
  "the EPIC manifest. Everything that needs annotation is refused rather than",
  "silently computed against the wrong array."
)

## Refuse an array we cannot annotate. Called from the manifest and from the
## model, so an EWAS with masking switched off still stops here rather than
## running to completion with no annotation behind it.
mp_require_array <- function(pgx) {
  shiny::validate(shiny::need(mp_array_type(pgx) != "EPICv2", MP_EPICV2_MSG))
}

mp_manifest <- function(pgx) {
  type <- mp_array_type(pgx)
  mp_require_array(pgx)
  ann <- mp_read_manifest(MP_ANNO_PKG[[type]])
  if (is.null(ann)) warning("[methylome] no manifest available; SNP masking disabled")
  ann
}

## X and Y probes. Read from the manifest rather than from pgx$genes, whose
## chr column is a cytoband ("Xp22.2") in every playbase-built pgx and would
## need parsing; the manifest already says "chrX".
##
## Why this is a mask at all: in a mixed-sex cohort every X probe separates
## the sexes by dosage compensation alone, so any phenotype even mildly
## correlated with sex picks up a block of spurious X hits. Adjusting for sex
## does not fix it - the effect is not additive on the beta scale.
mp_sexchr_probes <- function(pgx) {
  ann <- mp_manifest(pgx)
  if (is.null(ann) || !"chr" %in% colnames(ann)) return(character(0))
  rownames(ann)[ann$chr %in% c("chrX", "chrY")]
}

## Which probes survive the selected masks.
##
## Note the default: SNP and cross-reactive are on because they are
## unconditionally wrong probes - they measure genotype or the wrong locus
## whatever the cohort. The sex-chromosome mask is conditionally right: it is
## near-universal in mixed-sex population EWAS, and plainly wrong in a
## single-sex cohort or in any study whose question is sex itself, where it
## would silently delete the signal. So it is offered, unticked.
mp_masked_probes <- function(pgx, mask = c("snp", "xreactive")) {
  out <- character(0)
  if ("snp" %in% mask) out <- c(out, mp_snp_probes(pgx))
  if ("xreactive" %in% mask) out <- c(out, mp_xreactive_probes())
  if ("sexchr" %in% mask) out <- c(out, mp_sexchr_probes(pgx))
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

## Numeric sample columns usable as a continuous EWAS outcome. A large share
## of published EWAS regress on an exposure rather than compare two groups -
## age, BMI, gestational age, pack-years, dose - and limma fits those with the
## same machinery, so refusing them would be an arbitrary restriction.
## Two-valued numerics are excluded: those are a two-group contrast wearing
## numbers, and belong on the contrast path where delta-beta is meaningful.
mp_continuous_vars <- function(pgx) {
  s <- pgx$samples
  if (is.null(s)) return(character(0))
  keep <- vapply(colnames(s), function(k) {
    v <- suppressWarnings(as.numeric(as.character(s[[k]])))
    sum(!is.na(v)) >= 6 && stats::sd(v, na.rm = TRUE) > 0 &&
      length(unique(v[!is.na(v)])) > 2
  }, logical(1))
  sort(colnames(s)[keep])
}

## Categorical sample columns usable as an EWAS outcome directly.
##
## playbase builds every contrast as a two-group comparison, so a three-level
## phenotype (disease stage, treatment arm, tissue) is not merely refused by
## the model - it is unreachable from the picker. Offering the column itself is
## the only way into the F-test path.
##
## Two-level columns are included even though a contrast usually exists for
## them: it keeps one code path (the k = 2 branch is what a two-level column
## takes), and a sample column with no matching contrast would otherwise be
## untestable. At least three samples per level, or the group mean the F-test
## compares is noise.
mp_categorical_vars <- function(pgx, max_levels = 8, min_per_level = 3) {
  s <- pgx$samples
  if (is.null(s)) return(character(0))
  cont <- mp_continuous_vars(pgx)
  keep <- vapply(colnames(s), function(k) {
    if (k %in% cont) return(FALSE)
    v <- as.character(s[[k]])
    tb <- table(v[!is.na(v) & v != ""])
    length(tb) >= 2 && length(tb) <= max_levels && all(tb >= min_per_level)
  }, logical(1))
  sort(colnames(s)[keep])
}

## The outcome picker offers contrasts, continuous and categorical variables in
## one list, so which path to take is decided by what was picked. A contrast
## wins on a name clash - it is the more specific object.
mp_is_continuous <- function(pgx, sel) {
  if (is.null(sel) || !nzchar(sel)) return(FALSE)
  if (sel %in% colnames(pgx$contrasts)) return(FALSE)
  sel %in% mp_continuous_vars(pgx)
}

mp_is_categorical <- function(pgx, sel) {
  if (is.null(sel) || !nzchar(sel)) return(FALSE)
  if (sel %in% colnames(pgx$contrasts)) return(FALSE)
  sel %in% mp_categorical_vars(pgx)
}

#' Fit the differential-methylation model.
#'
#' @param contrast a contrast name, or the name of a continuous or categorical
#'   sample column.
#' @param n_sv surrogate variables to append to the design: 0 off, -1 estimate
#'   with sva::num.sv, a positive number to force that many.
#' @return list(table, formula, n, groups, strata, dropped_covars, masked, ...)
mp_fit_ewas <- function(pgx, contrast, covars = character(0),
                        cellfracs = NULL, mask = character(0), n_sv = 0) {
  mp_require_methylomics(pgx)
  mp_require_array(pgx)
  shiny::validate(shiny::need(!is.null(contrast) && nzchar(contrast),
                              "No outcome selected."))

  beta <- mp_beta(pgx)
  continuous <- mp_is_continuous(pgx, contrast)

  if (continuous) {
    x <- suppressWarnings(as.numeric(as.character(pgx$samples[[contrast]])))
    names(x) <- rownames(pgx$samples)
    x <- x[!is.na(x)]
    ss <- intersect(names(x), colnames(beta))
    shiny::validate(shiny::need(length(ss) >= 6,
      "Too few samples carry a value for this variable."))
    dat <- data.frame(.x = as.numeric(x[ss]), row.names = ss)
  } else {
    grp <- if (mp_is_categorical(pgx, contrast)) {
      g <- as.character(pgx$samples[[contrast]])
      names(g) <- rownames(pgx$samples)
      g[!is.na(g) & g != ""]
    } else {
      mp_contrast_groups(pgx, contrast)
    }
    shiny::validate(shiny::need(!is.null(grp), "This contrast has no group labels."))
    ss <- intersect(names(grp), colnames(beta))
    grp <- droplevels(factor(grp[ss]))
    shiny::validate(shiny::need(nlevels(grp) >= 2,
      "This outcome has only one group. Pick a contrast or a variable that varies."))
    dat <- data.frame(.group = grp, row.names = ss, stringsAsFactors = FALSE)
  }

  ## Covariate frame, aligned to the outcome's samples. The outcome itself is
  ## never also a covariate.
  used <- character(0); dropped <- character(0)
  for (k in setdiff(covars, contrast)) {
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
  if (continuous) {
    shiny::validate(shiny::need(stats::sd(dat$.x) > 0,
      "This variable is constant once incomplete covariates are dropped."))
  } else {
    dat$.group <- droplevels(dat$.group)
    shiny::validate(shiny::need(nlevels(dat$.group) >= 2,
      "Only one group has samples left once incomplete covariates are dropped."))
  }
  ## Three or more groups is a moderated F-test, not a difference of two means:
  ## no signed effect and no single t exists, which several panels care about.
  anova_fit <- !continuous && nlevels(dat$.group) >= 3

  ## A continuous outcome keeps the intercept and tests its own slope; a
  ## contrast drops it and tests the difference of the two group means.
  rhs <- c(if (continuous) ".x" else "0 + .group", used)
  design <- stats::model.matrix(
    stats::as.formula(paste("~", paste(rhs, collapse = " + "))), data = dat)
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
  M <- playbase.epigenetics::betaToM(B)

  ## ------------------------------------------------------- latent factors --
  ## Surrogate variables: the only confounding control left when no
  ## deconvolution reference exists for the tissue. The outcome is protected by
  ## sitting in mod and not in mod0, so sva estimates structure orthogonal to
  ## the question rather than removing the answer.
  sv_note <- NULL
  if (n_sv != 0) {
    shiny::validate(shiny::need(requireNamespace("sva", quietly = TRUE),
      "The sva package is not installed, so latent-factor adjustment is unavailable."))
    outcome_cols <- if (continuous) ".x" else grep("^\\.group", colnames(design), value = TRUE)
    mod0 <- design[, setdiff(colnames(design), outcome_cols), drop = FALSE]
    ## The contrast design is fitted without an intercept (~ 0 + .group), so
    ## dropping the group columns can leave a null model with no constant term.
    if (!any(apply(mod0, 2, function(z) length(unique(z)) == 1))) {
      mod0 <- cbind(Intercept = 1, mod0)
    }
    sv <- tryCatch({
      k <- if (n_sv < 0) sva::num.sv(M, design, method = "be") else as.integer(n_sv)
      ## Cap: ten SVs is already more latent structure than a typical cohort
      ## supports, and every SV costs a residual degree of freedom.
      k <- min(k, 10L, nrow(dat) - ncol(design) - 2L)
      if (k < 1) NULL else {
        ## sva narrates its IRW iterations on stdout; the progress bar already
        ## says what is happening.
        invisible(utils::capture.output(
          obj <- sva::sva(M, mod = design, mod0 = mod0, n.sv = k)))
        obj$sv
      }
    }, error = function(e) e)
    shiny::validate(shiny::need(!inherits(sv, "error"),
      paste0("Surrogate-variable estimation failed on this design: ",
             if (inherits(sv, "error")) conditionMessage(sv) else "",
             ". Untick the latent-factor adjustment, or drop a covariate.")))
    if (is.null(sv)) {
      ## Degrade to the unadjusted model, but say so - a silent fall-back is
      ## exactly the fail-open this board exists to prevent.
      sv_note <- "no latent factors detected"
    } else {
      colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
      design <- cbind(design, sv)
      used <- c(used, colnames(sv))
      shiny::validate(shiny::need(nrow(dat) - ncol(design) >= 2,
        "Not enough residual degrees of freedom once the surrogate variables are added - drop a covariate."))
    }
  }

  if (continuous) {
    fit2 <- limma::eBayes(limma::lmFit(M, design), trend = TRUE)
    tt <- limma::topTable(fit2, coef = ".x", number = Inf, sort.by = "none")
    ## Effect size on the beta scale is the slope of beta on the outcome, so
    ## the same design is refitted on B. Reported per unit of the outcome -
    ## a delta-beta of 0.004 for age means 0.4% methylation per year.
    dbeta <- limma::lmFit(B, design)$coefficients[, ".x"]
    head <- paste0("~ ", contrast, " (continuous)")
    ## Tertiles, so the per-CpG and per-region panels have one code path
    ## whether the outcome is a contrast or a continuous exposure.
    br <- unique(c(-Inf, stats::quantile(dat$.x, c(1 / 3, 2 / 3), na.rm = TRUE), Inf))
    strata <- factor(cut(dat$.x, br,
                         labels = if (length(br) == 4) c("low", "mid", "high") else NULL))
    groups <- NULL
    desc <- sprintf("%s from %.4g to %.4g", contrast, min(dat$.x), max(dat$.x))
  } else {
    lv <- levels(dat$.group)
    gcol <- make.names(paste0(".group", lv))
    fit <- limma::lmFit(M, design)
    ## Every level against the first. With one contrast limma returns the
    ## moderated t; with k-1 it returns the moderated F of "any difference
    ## among the groups", which is the multi-level test. Testing the k group
    ## coefficients directly would instead ask whether every group mean is
    ## zero - true of no M-value anywhere.
    cm <- limma::makeContrasts(contrasts = paste0(gcol[-1], "-", gcol[1]),
                               levels = design)
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm), trend = TRUE)
    tt <- limma::topTable(fit2, number = Inf, sort.by = "none")
    if (anova_fit) {
      ## An F-test has no direction, so the effect size is the spread of the
      ## group means on the beta scale: unsigned, still a methylation
      ## difference, and still meaningful to the |delta-beta| filter.
      mu <- lapply(lv, function(l)
        rowMeans(B[, dat$.group == l, drop = FALSE], na.rm = TRUE))
      dbeta <- do.call(pmax, mu) - do.call(pmin, mu)
      head <- sprintf("~ %s (F-test across %d groups)", contrast, length(lv))
    } else {
      dbeta <- rowMeans(B[, dat$.group == lv[2], drop = FALSE], na.rm = TRUE) -
               rowMeans(B[, dat$.group == lv[1], drop = FALSE], na.rm = TRUE)
      head <- paste0("~ ", lv[2], " vs ", lv[1])
    }
    strata <- dat$.group
    groups <- table(dat$.group)
    desc <- NULL
  }
  names(strata) <- rownames(dat)

  tab <- data.frame(
    probe = rownames(tt),
    logFC_M = if (anova_fit) NA_real_ else tt$logFC,
    dbeta = round(dbeta[rownames(tt)], 4),
    p = tt$P.Value, q = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
  ## Only on the anova branch: the F is what was tested there, and there is no
  ## single log fold change across k-1 contrasts to report in its place.
  if (anova_fit) tab$Fstat <- tt$F

  list(
    table = tab,
    formula = paste0(head,
                     if (length(used)) paste0(" + ", paste(used, collapse = " + ")) else
                       "   (unadjusted)",
                     if (!is.null(sv_note)) paste0("   [", sv_note, "]")),
    ## dmrff needs the standard error, which limma does not return directly;
    ## se = |logFC / t|, so keep the moderated t alongside. NULL on the anova
    ## branch so nothing downstream mistakes an F for a signed t.
    t = if (anova_fit) NULL else tt$t,
    anova = anova_fit,
    n_sv = if (is.null(sv_note)) sum(grepl("^SV[0-9]+$", used)) else 0L,
    n = nrow(dat), groups = groups, desc = desc,
    continuous = continuous,
    x = if (continuous) stats::setNames(dat$.x, rownames(dat)) else NULL,
    strata = strata, samples = rownames(dat),
    dropped_covars = dropped, masked = length(masked),
    covars = used
  )
}
