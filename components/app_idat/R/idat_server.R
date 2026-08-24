##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## IDAT -> beta table. The reading itself lives in
## playbase.epigenetics::read_idats(); everything here is the app side that
## commit's plan named as out of scope for the package: stage the upload back
## into a directory minfi can resolve, run the read off the main thread, show
## the per-sample failures before the user commits, and hand back a bundle that
## drops straight into the upload board.

## T0 measured ~1.4 GB fixed + ~60 MB per sample inside one worker. Refusing a
## cohort is a message; OOM-ing the worker takes every other session with it.
## ponytail: hard cap. Chunked reading or an HDF5Array backend is the upgrade
## path if real cohorts start hitting it (playbase.epigenetics plan, D7).
IDAT_MAX_SAMPLES <- 96

## A 450k _Grn/_Red pair is ~16 MB; shiny.maxRequestSize is 999 MB
## (components/app/R/global.R). Stating the arithmetic beats an unexplained 502.
IDAT_MAX_UPLOAD_MB <- 999

## ---------------------------------------------------------------------------
## Help text. Single source: the app renders these behind its info buttons and
## presentation/build-idat-guide.py parses this same block into the printed
## guide, so the manual cannot drift from the product. Editing the guide HTML
## instead of this list is the one thing that breaks the arrangement.
##
## Every entry that describes something *computed* names the package and the
## function actually called, per presentation/GUIDE-GUIDE.md.
## ---------------------------------------------------------------------------
IDAT_HELP <- list(
  upload = list(
    title = "Uploading IDAT files",
    text = paste(
      "An Illumina scanner writes two files per sample, not one: a green",
      "channel and a red channel, named <b>&lt;prefix&gt;_Grn.idat</b> and",
      "<b>&lt;prefix&gt;_Red.idat</b>. Both are needed. Drop the files in",
      "loose, or a ZIP of them - a ZIP keeps subdirectories, which is how GEO",
      "ships them, and loose files land flat. Either resolves.",
      "<br><br>Samples are paired by filename: two files belong to the same",
      "sample when everything before the channel suffix is identical. Any",
      "prefix works - GSM accessions, Sentrix chip positions, your own names -",
      "so long as the two names differ only in <b>_Grn</b> versus <b>_Red</b>.",
      "Both <b>.idat</b> and <b>.idat.gz</b> are read.",
      "<br><br>A sample sheet CSV is optional. Include one and its",
      "<b>Sample_Name</b> column names your columns; leave it out and one is",
      "generated, with each sample named after its file prefix. A sample",
      "missing one of its two channels is listed as incomplete and skipped,",
      "never silently averaged in."
    ),
    methods = paste(
      "The sheet is parsed by <b>minfi::read.metharray.sheet()</b>, which knows",
      "the Illumina sheet dialects, then validated by",
      "<b>playbase.epigenetics</b>: Sentrix_ID/Slide and",
      "Sentrix_Position/Array naming variants are normalised, relative paths",
      "are resolved against the sheet rather than the working directory, and",
      "every row is checked for an existing Grn/Red pair before any file is",
      "read. Rows that do not resolve are reported, not fatal.",
      "A generated sheet deliberately carries no Slide or Array columns, only",
      "an absolute file prefix, because minfi rebuilds the prefix by pattern",
      "matching whenever those columns are present - omitting them is what",
      "makes arbitrary filenames work."
    )
  ),
  preprocessing = list(
    title = "Preprocessing",
    text = paste(
      "This is the only correction applied to your data. <b>Noob</b>, the",
      "default, removes background signal and the dye bias between the two",
      "colour channels using the array's own control probes; it works at any",
      "cohort size, including a single sample. <b>Funnorm</b> additionally",
      "fits control-probe components across samples to remove unwanted",
      "technical variation - it wants roughly ten samples or more, and warns",
      "below that, because under a cohort those components track noise rather",
      "than chemistry. <b>Raw</b> applies nothing.",
      "<br><br>Beyond this step no further normalization is applied to the",
      "beta values. Probe-design correction such as BMIQ is a separate choice",
      "you make when you load the table into the platform, so it is never",
      "applied twice."
    ),
    methods = paste(
      "<b>minfi::preprocessNoob()</b>, <b>minfi::preprocessFunnorm()</b> or",
      "<b>minfi::preprocessRaw()</b>, all called through",
      "<b>playbase.epigenetics::read_idats()</b>. Beta values are then taken",
      "with <b>minfi::getBeta()</b>. All three operate on raw intensities,",
      "before beta values exist, which is why they cannot be applied later to",
      "a beta table."
    )
  ),
  array = list(
    title = "Array type",
    text = paste(
      "Auto-detect reads the array type out of the files themselves and is",
      "the right setting almost always. Choosing 450K or EPIC explicitly does",
      "not override detection - it is checked against the files, and a",
      "mismatch is refused rather than annotated against the wrong manifest.",
      "<br><br>Array detection is the one place where a wrong answer is",
      "invisible: the wrong manifest gives the wrong probes, which gives beta",
      "values that look entirely normal and are wrong. So anything that is",
      "not 450K or EPIC v1 is refused by name rather than guessed at.",
      "<b>EPIC v2 is not supported</b> and is rejected with a message saying",
      "so."
    ),
    methods = paste(
      "The identifier comes from <b>minfi::annotation()</b> on the loaded",
      "object and is mapped by <b>playbase.epigenetics::array_type()</b>.",
      "Probe annotation is then",
      "<b>playbase.epigenetics::annotate_methylomics()</b> against the hg19",
      "manifests - <b>IlluminaHumanMethylation450kanno.ilmn12.hg19</b> or",
      "<b>IlluminaHumanMethylationEPICanno.ilm10b4.hg19</b>."
    )
  ),
  detect_p = list(
    title = "Detection p-value cut",
    text = paste(
      "A detection p-value asks whether a probe's signal can be told apart",
      "from background on that particular array. Probes above this cut are",
      "set to <b>NA</b> for the sample that failed them, rather than being",
      "given a number that is really noise.",
      "<br><br>Masking rather than deleting is deliberate: it keeps the table",
      "rectangular, and one bad probe on one array does not cost you that",
      "probe in every other sample. The default of 0.01 is the standard cut.",
      "Setting it to 1 masks nothing and produces a table with no NA values,",
      "at the price of keeping background noise as though it were signal."
    ),
    methods = paste(
      "<b>minfi::detectionP()</b> on the raw object, aligned to the beta",
      "matrix. The per-sample share of failing probes is kept and reported as",
      "the <b>Detection fail</b> column of the QC table - on a healthy sample",
      "it sits well under one percent."
    )
  ),
  max_fail = list(
    title = "Dropping cohort-wide failures",
    text = paste(
      "A probe that fails in more than this fraction of your samples is not",
      "masked but removed entirely - it is dead across the cohort rather than",
      "unlucky on one array. Deleting is reserved for that case because every",
      "probe removed costs coverage of the epigenetic clocks and of any",
      "region you later test.",
      "<br><br>Sex-chromosome probes are the one exception and are never",
      "dropped by this rule. In an all-female cohort most of chromosome Y",
      "fails by definition, and dropping it would silently break the sex",
      "check that the QC table depends on. So the NA values in a finished",
      "table are mostly chrX and chrY in samples of the other sex, not",
      "scattered noise."
    ),
    methods = paste(
      "Applied inside <b>playbase.epigenetics::read_idats()</b>. Sex-probe",
      "identity comes from the chromosome column of the array's annotation",
      "package. On the six-sample 450K reference set at defaults: 18,027 of",
      "485,512 probes dropped, 467,485 kept, and 1,989 masked cells - 0.071",
      "percent of the table, every one of them in the 817 sex-chromosome",
      "probes held back by this exemption."
    )
  ),
  convert = list(
    title = "Running the conversion",
    text = paste(
      "Conversion runs off the main thread, so the application stays",
      "responsive and other users are unaffected. The progress bar reports",
      "the phase the reader is actually in - reading pairs, raw QC, detection",
      "p-values and normalization, masking, annotation - alongside an elapsed",
      "clock.",
      "<br><br>The percentage is an estimate weighted by measured phase costs;",
      "the clock beside it is the exact number. Annotation is the long step",
      "and barely grows with cohort size, while reading grows per sample.",
      "Three 450K samples take about 100 seconds end to end, most of it",
      "annotation.",
      "<br><br>Conversions are capped at 96 samples. The reader needs roughly",
      "1.4 GB of memory before it starts plus about 60 MB per sample, so a",
      "cap is a message where an over-large cohort would be a failure. Uploads",
      "are separately capped at 999 MB, and an uncompressed 450K pair is",
      "about 16 MB."
    ),
    methods = paste(
      "<b>playbase.epigenetics::read_idats()</b> inside a",
      "<b>promises::future_promise()</b> worker, driven by a Shiny",
      "<b>ExtendedTask</b>. The progress feed is the reader's own phase",
      "messages, streamed to the browser as they happen rather than",
      "reconstructed from a separate timer."
    )
  ),
  beta_tab = list(
    title = "Beta values",
    text = paste(
      "The converted table: one row per CpG probe, one column per sample,",
      "each value the fraction of DNA methylated at that site on a 0 to 1",
      "scale. The caption gives the true dimensions; the view shows the first",
      "5,000 rows because a browser cannot usefully scroll half a million.",
      "The full table is in the download.",
      "<br><br>NA entries are probes masked by the detection cut - see the",
      "Advanced settings. They are carried through rather than filled in,",
      "because an invented value is worse than a missing one."
    ),
    methods = NULL
  ),
  qc_tab = list(
    title = "Quality control",
    text = paste(
      "Two views that answer different questions. The scatter plots median",
      "methylated against median unmethylated intensity, one point per",
      "sample. This is the check only raw files can give you: a dead sample",
      "and a healthy one can produce beta distributions that look alike, but",
      "they cannot produce alike intensities. Samples below the standard",
      "diagonal are drawn in red and labelled.",
      "<br><br>The table below is the same per-sample ledger the Methylome",
      "Profiler shows, so a sample judged here is judged the same way there.",
      "Bimodality is the fraction of probes outside the intermediate band and",
      "should sit near 0.9. Bisulfite conversion below roughly 80 percent",
      "means incomplete conversion. Detection fail above about one percent",
      "means a poor sample. Verdict is relative to your cohort, not absolute.",
      "<br><br>The four raw-intensity columns exist only before beta values",
      "are computed. Convert the files elsewhere and they are gone; this is",
      "the only point at which they can be captured."
    ),
    methods = paste(
      "Median intensities from <b>minfi::getQC()</b>, bisulfite conversion",
      "from <b>wateRmelon::bscon()</b>, both on the raw object. The scatter is",
      "<b>playbase.epigenetics::plotIntensityQC()</b>, with the standard",
      "cutoff at mean log2 intensity 10.5. The remaining columns are",
      "<b>playbase.epigenetics::sample_qc()</b>, the same function the",
      "Methylome Profiler and the platform's upload step both call."
    )
  ),
  log_tab = list(
    title = "Log and failures",
    text = paste(
      "Which samples did not make it, and why. One corrupt file costs its own",
      "sample and nothing else - the rest of the chip converts, and the",
      "failure is named here with its reason. Check this before using the",
      "output: a cohort quietly three samples short is the failure mode this",
      "tab exists to prevent, and the beta table alone will not show it.",
      "<br><br>Mixed array types in one upload are refused outright rather",
      "than truncated to their shared probes, naming both groups. A table",
      "silently reduced to the intersection of two arrays is not something",
      "anyone can interpret later.",
      "<br><br>Below the table is the reader's own log: array type detected,",
      "probes masked and dropped, annotation coverage."
    ),
    methods = NULL
  ),
  send = list(
    title = "Sending it straight to the upload wizard",
    text = paste(
      "The alternative to downloading and re-uploading. <b>Send to upload</b>",
      "opens the platform's upload wizard with the beta table and the sample",
      "table already in it, and the datatype and array already set - so the",
      "conversion never leaves the browser as a file at all.",
      "<br><br>Everything the wizard normally asks for is still asked for. You",
      "name the dataset, fill in the groups you want to compare, and choose",
      "the normalization - including the probe-design correction this",
      "application deliberately leaves to that step. Nothing is computed",
      "behind you.",
      "<br><br>The download is still there and produces the same tables. Use",
      "it when the beta values are going somewhere other than this platform,",
      "or when you want a copy before committing to an import."
    ),
    methods = paste(
      "The button writes the conversion into the same channel the upload",
      "board already watches to pre-fill itself, so there is no second route",
      "into the wizard to keep in step with the first. The array type sent",
      "with it is the one <b>playbase.epigenetics::read_idats()</b> detected",
      "in the files, not a default - the platform's probe annotation is",
      "manifest-specific, so 450K settings over EPIC data would be wrong in a",
      "way nothing downstream would show you."
    )
  ),
  download = list(
    title = "The download",
    text = paste(
      "One ZIP with three files, named for where they go next.",
      "<b>counts.csv</b> is the beta table - it is what the platform's upload",
      "step expects for a methylomics dataset, and beta values are detected",
      "and passed through unchanged. <b>samples.csv</b> lists the sample",
      "names as a starting point for your phenotype table; add your groups to",
      "it. <b>qc.csv</b> is the ledger from the QC tab, including the",
      "raw-intensity columns that cannot be recovered later.",
      "<br><br>Load counts.csv and samples.csv through the normal upload",
      "screen as datatype methylomics, or skip the round trip entirely with",
      "<b>Send to upload</b>. Probe-design correction is offered there, which",
      "is why it is not applied here.",
      "<br><br>The table is large: six 450K samples are about 57 MB of CSV,",
      "25 MB zipped."
    ),
    methods = NULL
  )
)

## Info button for one help entry, rendered as a popover.
idat_info <- function(key) {
  h <- IDAT_HELP[[key]]
  body <- h$text
  if (!is.null(h$methods)) {
    body <- paste0(body, "<br><br><b>How it is computed.</b> ", h$methods)
  }
  bslib::popover(
    shiny::icon("circle-info", style = "color: #ccc; cursor: pointer;"),
    shiny::HTML(body),
    title = h$title,
    options = list(customClass = "idat-help")
  )
}

#' Array platform choices; "auto" lets read_idats() detect it from the files.
idat_platform_choices <- function() {
  c("Auto-detect" = "auto", "450K array" = "450K array", "EPIC array" = "EPIC array")
}

#' Intensity-level preprocessing choices, as read_idats() names them.
idat_normalization_choices <- function() {
  c(
    "Noob (background + dye)" = "noob",
    "Funnorm (needs ~10+ samples)" = "funnorm",
    "Raw (none)" = "raw"
  )
}

## An empty inventory, in the one shape every caller expects back.
idat_no_samples <- function() {
  data.frame(
    sample = character(0), grn = character(0), red = character(0),
    complete = logical(0), stringsAsFactors = FALSE
  )
}

## Pair the .idat files under `dir` into samples. minfi resolves IDATs by
## filename convention, so this works off names on disk, never off upload order.
idat_inventory <- function(dir) {
  paths <- list.files(dir,
    pattern = "\\.idat(\\.gz)?$", recursive = TRUE,
    ignore.case = TRUE, full.names = TRUE
  )
  if (!length(paths)) {
    return(idat_no_samples())
  }
  base <- sub("_(Grn|Red)\\.idat(\\.gz)?$", "", basename(paths), ignore.case = TRUE)
  chan <- ifelse(grepl("_Grn\\.idat(\\.gz)?$", paths, ignore.case = TRUE), "grn", "red")
  samples <- sort(unique(base))
  pick <- function(s, ch) {
    hit <- paths[base == s & chan == ch]
    if (length(hit)) hit[1] else NA_character_
  }
  out <- data.frame(
    sample = samples,
    grn = vapply(samples, pick, character(1), "grn"),
    red = vapply(samples, pick, character(1), "red"),
    stringsAsFactors = FALSE
  )
  out$complete <- !is.na(out$grn) & !is.na(out$red)
  rownames(out) <- NULL
  out
}

#' Stage a Shiny upload into a directory minfi can read.
#'
#' fileInput hands back \code{data.frame(name, datapath)} where datapath is a
#' meaningless temp name, so the real filenames have to be restored on disk
#' before minfi sees them. ZIPs are unpacked with their structure intact - GEO
#' ships \code{<Slide>/<Slide>_<Array>_Grn.idat}, and read_idat_targets()
#' already resolves idats parked in a subdirectory.
#'
#' @param files What \code{input$idat} gives, or NULL.
#' @param exdir Directory to stage into; created if absent.
#'
#' @return \code{list(dir, inventory, sheet)}. \code{sheet} is the user's own
#'   sample sheet if the upload carried one (a name matching sheet|sample wins
#'   over any other CSV), else NA.
idat_stage_upload <- function(files, exdir = tempfile("idat")) {
  if (is.null(files) || nrow(files) == 0) {
    return(list(dir = exdir, inventory = idat_no_samples(), sheet = NA_character_))
  }
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(files))) {
    if (grepl("\\.zip$", files$name[i], ignore.case = TRUE)) {
      utils::unzip(files$datapath[i], exdir = exdir)
    } else {
      ## basename(): a browser can send a relative path, and a crafted one
      ## ("../../x") must not escape the staging directory.
      dest <- file.path(exdir, basename(files$name[i]))
      file.copy(files$datapath[i], dest, overwrite = TRUE)
    }
  }

  ## The user's own sheet, if any. A name saying "sheet" or "sample" beats a
  ## stray readme.csv; ours is written later and never lives here yet.
  csv <- list.files(exdir,
    pattern = "\\.csv$", recursive = TRUE,
    ignore.case = TRUE, full.names = TRUE
  )
  ## Match on the basename but keep the full path - read_idats() resolves
  ## everything relative to the sheet, so a bare filename resolves to nothing.
  hit <- grep("sheet|sample", basename(csv), ignore.case = TRUE)
  sheet <- if (length(hit)) csv[hit[1]] else if (length(csv)) csv[1] else NA_character_

  list(dir = exdir, inventory = idat_inventory(exdir), sheet = sheet)
}

#' Write a minimal sample sheet for an upload that carried none.
#'
#' Deliberately no Slide/Array columns: minfi::read.metharray.sheet() rebuilds
#' Basename by grepping "<Slide>_<Array>_Grn.idat" whenever an Array column
#' exists, but skips that block entirely when df$Array is NULL - so the
#' absolute Basename written here survives, and any filename convention works.
#'
#' @param inventory The \code{inventory} frame from \code{idat_stage_upload()}.
#' @param dir Directory to write the sheet into.
#' @return Path to the written sheet.
idat_write_sheet <- function(inventory, dir) {
  ok <- inventory[inventory$complete, , drop = FALSE]
  if (!nrow(ok)) {
    stop("idat_write_sheet: no complete IDAT pairs to write a sheet for.", call. = FALSE)
  }
  path <- file.path(dir, "_playground_sheet.csv")
  utils::write.csv(
    data.frame(
      Sample_Name = ok$sample,
      Basename = normalizePath(sub("_Grn\\.idat(\\.gz)?$", "", ok$grn, ignore.case = TRUE),
        mustWork = FALSE
      ),
      stringsAsFactors = FALSE
    ),
    path,
    row.names = FALSE
  )
  path
}

## A future worker loads the whole Bioconductor stack on first use, so the raw
## message stream is ~110 lines of "Loading required package:" and masked-object
## banners around ~15 that say anything. Every message the playbase/minfi layers
## emit on purpose is tagged "[name] ...", so keep those and the problems.
## Two warnings are expected and would read as problems: minfi complains that
## our synthesized sheet has no Slide/Array columns (it deliberately does not -
## see idat_write_sheet), and playbase reports its own dplyr/data.table import
## masking on load.
IDAT_LOG_EXPECTED <- "Could not infer (array|slide) name|replacing previous import"

idat_log_lines <- function(x) {
  keep <- grepl("^\\[|error|warning", x, ignore.case = TRUE)
  x[keep & !grepl(IDAT_LOG_EXPECTED, x)]
}

## read_idats() announces each phase as it starts, so the progress feed is its
## own message stream read back - no second source of truth to drift.
##
## `pct` is where the bar sits once that phase has been announced, timed
## against a real 3-sample GSE68777 run (100 s wall clock): startup 5 s,
## read 6 s, raw QC 6 s, detection-p + noob 18 s, mask 1 s, annotate 60 s,
## ledger 4 s. Annotation is the long tail, so it deliberately parks the bar
## at 60 rather than 80 - a bar sitting near-done for a minute reads as a hang.
##
## ponytail: a calibration knob, not a measurement. Phase costs scale
## differently with cohort size (reading grows per sample, annotation barely
## does), so the bar is an estimate and the elapsed clock beside it is the
## honest number. Re-time it if the shape looks wrong.
##
## Note the fourth row covers two steps: preprocessNoob/funnorm runs straight
## after detectionP and emits nothing, so it has no phase line of its own.
IDAT_PHASES <- data.frame(
  pattern = c(
    "read_idat_targets\\] Parsed",
    "read_idats\\] Reading",
    "read_idats\\] Raw QC",
    "read_idats\\] Detection p-values",
    "mask_by_detection\\]",
    "read_idats\\] Annotating",
    "annotate_methylomics\\] Annotation completed"
  ),
  label = c(
    "Parsing the sample sheet",
    "Reading IDAT pairs",
    "Raw intensity QC (intensities, bisulfite conversion)",
    "Detection p-values, then normalization",
    "Masking failed probes",
    "Annotating probes (the long step)",
    "Building the QC ledger"
  ),
  pct = c(5, 20, 35, 50, 55, 60, 95),
  stringsAsFactors = FALSE
)

#' Read the worker's progress file back into a phase, a percentage and a clock.
#'
#' @param path File the worker appends its phase lines to; may not exist yet.
#' @param elapsed Seconds since the task was invoked.
#' @return \code{list(pct, label, elapsed)}.
idat_progress <- function(path, elapsed = 0) {
  lines <- if (file.exists(path)) readLines(path, warn = FALSE) else character(0)
  hit <- 0L
  for (i in seq_len(nrow(IDAT_PHASES))) {
    if (any(grepl(IDAT_PHASES$pattern[i], lines))) hit <- i
  }
  list(
    ## Before the first phase line the worker is loading the Bioconductor
    ## stack, which on a cold worker is the slowest silent stretch of the run.
    pct = if (hit) IDAT_PHASES$pct[hit] else 2,
    label = if (hit) IDAT_PHASES$label[hit] else "Starting up (loading the methylation stack)",
    elapsed = elapsed
  )
}

#' Run read_idats() and capture everything it says on the way.
#'
#' Called inside a future worker, where an error is a dead promise and a
#' message() is lost entirely - so both are caught and returned as data.
#' read_idats() phrases its stops for a user (EPIC v2, mixed array types,
#' nothing readable), so `msg` is safe to show verbatim.
#'
#' @return \code{list(ok, res, log, msg)}.
idat_convert <- function(sheet, method = "noob", detect_p = 0.01,
                         max_fail = 0.05, array = "auto", progress_file = NULL,
                         sheet_csv = NA_character_) {
  log <- character(0)
  msg <- NA_character_

  ## Every phase line is appended to progress_file as it happens, so the UI can
  ## tail it while the worker is still busy - the log returned at the end is
  ## the same stream, just too late to watch.
  note <- function(x) {
    line <- sub("\n$", "", conditionMessage(x))
    log <<- c(log, line)
    if (!is.null(progress_file) && length(idat_log_lines(line))) {
      cat(line, "\n", sep = "", file = progress_file, append = TRUE)
    }
  }

  ## Handle the conditions rather than sinking them: future's multisession
  ## worker installs its own calling handler and muffles message(), so a
  ## capture.output(type = "message") there comes back empty. A calling
  ## handler runs first and works in both the worker and a direct call.
  res <- withCallingHandlers(
    {
      out <- tryCatch(
        playbase.epigenetics::read_idats(
          input = sheet,
          method = method,
          detect_p = detect_p,
          max_fail = max_fail,
          array = if (identical(array, "auto")) NULL else array
        ),
        error = function(e) {
          msg <<- conditionMessage(e)
          NULL
        }
      )
      ## The QC ledger the package documents as the call to make. A failure
      ## here must not lose the betas we just spent minutes computing.
      if (!is.null(out)) {
        out$array <- idat_array_type(log)
        out$samples <- idat_samples(out, sheet_csv)
        out$ledger <- tryCatch(
          cbind(playbase.epigenetics::sample_qc(out$beta, out$annot), out$qc),
          error = function(e) {
            log <<- c(log, paste("sample_qc failed:", conditionMessage(e)))
            out$qc
          }
        )
      }
      out
    },
    message = function(m) {
      note(m)
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      log <<- c(log, paste("Warning:", conditionMessage(w)))
      tryInvokeRestart("muffleWarning")
    }
  )

  list(ok = !is.null(res), res = res, log = idat_log_lines(log), msg = msg)
}

## read_idats() announces the array it detected exactly once. Taking it from
## there rather than re-deriving it keeps one answer in play: the platform's
## methylomics annotation is manifest-specific, so a 450K default over EPIC
## data is the same silent wrong answer array_type() refuses to make.
idat_array_type <- function(log) {
  hit <- regmatches(log, regexpr("Array type: [^;]+", log))
  if (!length(hit)) {
    return(NA_character_)
  }
  trimws(sub("Array type: ", "", hit[1]))
}

#' The sample table for a conversion: the user's sheet where there was one.
#'
#' Rows are aligned to the beta columns and named by them, which is the shape
#' both the upload board (\code{intersect(colnames(X), rownames(Y))}) and the
#' downloaded samples.csv need. Columns that only locate files on disk are
#' dropped - they are not phenotypes and mean nothing on another machine.
#'
#' @return data.frame with \code{rownames == colnames(beta)}; a single
#'   \code{sample} column when the upload carried no usable sheet.
idat_samples <- function(res, sheet = NA_character_) {
  ids <- colnames(res$beta)
  bare <- data.frame(sample = ids, row.names = ids, stringsAsFactors = FALSE)
  if (is.na(sheet) || !file.exists(sheet)) {
    return(bare)
  }
  df <- tryCatch(utils::read.csv(sheet, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || !nrow(df)) {
    return(bare)
  }

  ## read_idats() names the beta columns from Sample_Name when the sheet has
  ## one, so that is the join key; without it there is nothing to join on.
  key <- grep("^sample_name$", colnames(df), ignore.case = TRUE)
  if (!length(key)) {
    return(bare)
  }
  rownames(df) <- make.unique(as.character(df[[key[1]]]))
  if (!any(ids %in% rownames(df))) {
    return(bare)
  }

  drop <- grepl("^(basename|slide|array|sentrix|sample_name|sample_well|pool_id|filenames?|path)$",
    colnames(df), ignore.case = TRUE)
  df <- df[, !drop, drop = FALSE]
  if (!ncol(df)) {
    return(bare)
  }
  out <- df[ids, , drop = FALSE]
  rownames(out) <- ids
  out
}

#' Bundle a conversion into the three CSVs the upload board wants, zipped.
#'
#' counts.csv holds beta values: mToBeta() detects a 0-1 matrix and returns it
#' unchanged, so betas are what belongs in a methylomics counts.csv.
idat_bundle <- function(res, file) {
  dir <- tempfile("bundle")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  utils::write.csv(res$beta, file.path(dir, "counts.csv"))
  utils::write.csv(res$ledger, file.path(dir, "qc.csv"))
  ## The same table the Send to upload button hands over, so the two routes
  ## into the platform cannot disagree. Row names are the sample ids, which is
  ## what the upload board matches counts columns against.
  utils::write.csv(res$samples, file.path(dir, "samples.csv"))

  ## zip() is path-relative: run it from the bundle dir so the archive holds
  ## three files, not three nested directories.
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE, after = FALSE)
  utils::zip(file, c("counts.csv", "samples.csv", "qc.csv"), flags = "-q")
  invisible(file)
}

## One DT config for all four tables here, so they scroll and size alike.
idat_dt <- function(df, rownames = FALSE, ...) {
  DT::datatable(
    data = df,
    rownames = rownames,
    extensions = c("Buttons", "Scroller"),
    plugins = "scrollResize",
    selection = "none",
    options = list(
      dom = "lfrtip", scroller = TRUE, scrollX = TRUE,
      scrollY = "calc(100vh - 260px)", scrollResize = TRUE
    ),
    ...
  ) %>% DT::formatStyle(columns = 0, target = "row", lineHeight = "90%")
}

#' The application server-side logic
#'
#' @param id Shiny module id.
#' @param recompute_pgx The app-wide reactiveVal the upload board watches. Set
#'   it and the board opens pre-filled - see the Send to upload button below.
#'   NULL (the standalone test app) hides that button.
#' @export
idat_server <- function(id, recompute_pgx = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    ## An IDAT upload is GB-scale; tempdir() outlives the session otherwise.
    staged_dirs <- shiny::reactiveVal(character(0))
    session$onSessionEnded(function() {
      unlink(shiny::isolate(staged_dirs()), recursive = TRUE, force = TRUE)
    })

    staged <- shiny::reactive({
      shiny::req(input$idat)
      st <- idat_stage_upload(input$idat)
      staged_dirs(c(shiny::isolate(staged_dirs()), st$dir))
      st
    })

    result <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$idat, result(NULL))

    ## Where the worker appends its phase lines; re-made per conversion so a
    ## second run does not read the first one's phases back.
    progress_file <- shiny::reactiveVal(tempfile("idat_progress"))
    started_at <- shiny::reactiveVal(NULL)

    convert_task <- shiny::ExtendedTask$new(
      function(sheet, method, detect_p, max_fail, array, progress_file) {
        promises::future_promise(
          idat_convert(sheet, method, detect_p, max_fail, array, progress_file,
            sheet_csv = sheet)
        )
      }
    ) |> bslib::bind_task_button("convert")

    ## Live feed. read_idats() takes minutes and says nothing to the browser
    ## otherwise, so tail the phase file once a second while the task runs.
    output$progress <- shiny::renderUI({
      if (!identical(convert_task$status(), "running")) {
        return(NULL)
      }
      shiny::invalidateLater(1000)
      elapsed <- round(as.numeric(difftime(Sys.time(), started_at(), units = "secs")))
      p <- idat_progress(progress_file(), elapsed)
      shiny::div(
        class = "mb-2",
        shiny::div(
          class = "progress", style = "height: 6px;",
          shiny::div(
            class = "progress-bar progress-bar-striped progress-bar-animated",
            style = sprintf("width: %d%%;", p$pct),
            role = "progressbar"
          )
        ),
        shiny::div(
          style = "font-size: 11px; color: #888; margin-top: 4px;",
          sprintf("%d%% - %s", p$pct, p$label),
          shiny::div(sprintf(
            "%d:%02d elapsed", p$elapsed %/% 60, p$elapsed %% 60
          ))
        )
      )
    })

    output$files_summary <- shiny::renderUI({
      st <- staged()
      inv <- st$inventory
      if (nrow(inv) == 0) {
        return(shiny::div(
          style = "font-size: 12px; color: #c66; margin-top: -8px; margin-bottom: 12px;",
          "No .idat files found in this upload."
        ))
      }
      n.ok <- sum(inv$complete)
      shiny::div(
        style = "font-size: 12px; color: #888; margin-top: -8px; margin-bottom: 12px;",
        sprintf("%d sample%s ready", n.ok, ifelse(n.ok == 1, "", "s")),
        if (any(!inv$complete)) {
          shiny::div(style = "color: #c66;", sprintf(
            "%d incomplete (missing Grn or Red)", sum(!inv$complete)
          ))
        },
        if (n.ok > IDAT_MAX_SAMPLES) {
          shiny::div(style = "color: #c66;", sprintf(
            "Over the %d-sample limit; split the cohort.", IDAT_MAX_SAMPLES
          ))
        },
        shiny::div(if (is.na(st$sheet)) {
          "No sample sheet - one will be generated."
        } else {
          paste("Sample sheet:", basename(st$sheet))
        })
      )
    })

    shiny::observeEvent(input$convert, {
      st <- staged()
      inv <- st$inventory
      n.ok <- sum(inv$complete)
      shiny::validate(
        shiny::need(n.ok > 0, "Please upload at least one complete pair of IDAT files."),
        shiny::need(n.ok <= IDAT_MAX_SAMPLES, sprintf(
          "%d samples exceeds the %d-sample limit for one conversion.",
          n.ok, IDAT_MAX_SAMPLES
        ))
      )
      ## Always hand read_idats() a CSV path, never the directory: it then
      ## narrows to that one file instead of rbinding every CSV in the upload.
      sheet <- if (is.na(st$sheet)) idat_write_sheet(inv, st$dir) else st$sheet

      result(NULL)
      progress_file(tempfile("idat_progress"))
      started_at(Sys.time())
      convert_task$invoke(
        sheet, input$method, input$detect_p, input$max_fail, input$platform,
        progress_file()
      )
    })

    shiny::observeEvent(convert_task$result(), {
      out <- convert_task$result()
      if (!out$ok) {
        shiny::showNotification(
          if (is.na(out$msg)) "IDAT conversion failed." else out$msg,
          type = "error", duration = NULL
        )
      }
      result(out)
    })

    ## ------------------------------- results -------------------------------

    output$table_area <- shiny::renderUI({
      out <- result()
      if (is.null(out) || !out$ok) {
        return(shiny::div(
          style = paste(
            "display: flex; flex-direction: column; align-items: center;",
            "justify-content: center; height: 100%; color: #aaa;",
            "font-size: 14px; text-align: center;"
          ),
          shiny::icon("table", class = "fa-2x", style = "margin-bottom: 10px;"),
          "Upload IDAT files and click Convert to see beta values here."
        ))
      }
      n.failed <- if (is.null(out$res$failed)) 0 else nrow(out$res$failed)
      ## The info buttons sit inside each panel, not in its tab title: a
      ## popover trigger in a tab trigger swallows the click meant to switch
      ## tabs, so the two controls fight over the same pixels.
      bslib::navset_card_tab(
        height = "100%",
        bslib::nav_panel(
          "Beta values",
          shiny::div(
            style = "font-size: 12px; color: #888; padding: 4px 0;",
            shiny::textOutput(session$ns("beta_caption"), inline = TRUE),
            idat_info("beta_tab")
          ),
          DT::dataTableOutput(session$ns("beta_table"))
        ),
        bslib::nav_panel(
          "QC",
          shiny::div(
            style = "font-size: 12px; color: #888; padding: 4px 0;",
            "Raw intensity and per-sample checks ", idat_info("qc_tab")
          ),
          shiny::plotOutput(session$ns("qc_plot"), height = "300px"),
          DT::dataTableOutput(session$ns("qc_table"))
        ),
        bslib::nav_panel(
          if (n.failed) sprintf("Log & failures (%d)", n.failed) else "Log & failures",
          shiny::div(
            style = "font-size: 12px; color: #888; padding: 4px 0;",
            "Samples that did not convert, and the reader's own log ",
            idat_info("log_tab")
          ),
          DT::dataTableOutput(session$ns("failed_table")),
          shiny::verbatimTextOutput(session$ns("log"))
        )
      )
    })

    ## A 485k-row table is a hung browser, not a feature. Preview, then download.
    beta_preview <- shiny::reactive({
      beta <- shiny::req(result()$res$beta)
      n <- min(nrow(beta), 5000)
      data.frame(
        probe = rownames(beta)[seq_len(n)],
        round(as.data.frame(beta[seq_len(n), , drop = FALSE]), 4),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    })

    output$beta_caption <- shiny::renderText({
      beta <- shiny::req(result()$res$beta)
      sprintf(
        "%d probes x %d samples. Showing the first %d - download for the full table.",
        nrow(beta), ncol(beta), min(nrow(beta), 5000)
      )
    })

    output$beta_table <- DT::renderDataTable(idat_dt(beta_preview()))

    output$qc_table <- DT::renderDataTable(
      idat_dt(shiny::req(result()$res$ledger), rownames = TRUE)
    )

    output$qc_plot <- shiny::renderPlot({
      qc <- shiny::req(result()$res$qc)
      ## plotIntensityQC carries @export in playbase.epigenetics R/plotting.R
      ## but is missing from that package's NAMESPACE - document() was not
      ## re-run in 122f67d. get() from the namespace reaches it either way, so
      ## this keeps working once the NAMESPACE is regenerated.
      get("plotIntensityQC", envir = asNamespace("playbase.epigenetics"))(qc)
    })

    output$failed_table <- DT::renderDataTable({
      failed <- result()$res$failed
      if (is.null(failed) || !nrow(failed)) {
        failed <- data.frame(sample = "-", reason = "All samples read successfully.")
      }
      idat_dt(failed)
    })

    output$log <- shiny::renderText(paste(shiny::req(result()$log), collapse = "\n"))

    ## ------------------------------ download -------------------------------

    ## Rendered only once there is something to download. downloadButton is an
    ## <a>, so an `enabled`/`disabled` attribute on it is inert - a button that
    ## looks clickable and serves an empty file is worse than no button.
    output$download_ui <- shiny::renderUI({
      out <- result()
      if (is.null(out) || !isTRUE(out$ok)) {
        return(NULL)
      }
      shiny::tagList(
        shiny::downloadButton(session$ns("download"), "Download ZIP",
          style = "width: 100%;"
        ),
        shiny::div(idat_info("download"), style = "text-align: center; margin-top: 4px;")
      )
    })

    output$download <- shiny::downloadHandler(
      filename = function() "idat_beta.zip",
      content = function(file) idat_bundle(result()$res, file)
    )

    ## ---------------------------- send to upload ---------------------------
    ## recompute_pgx is the channel the upload board already watches to
    ## pre-fill itself (upload_server.R: it sets uploaded$counts.csv,
    ## $samples.csv and the datatype from whatever is put here). Reusing it
    ## means no new plumbing through UploadBoard, and no second way for a
    ## dataset to enter the wizard.
    output$send_ui <- shiny::renderUI({
      out <- result()
      if (is.null(recompute_pgx) || is.null(out) || !isTRUE(out$ok)) {
        return(NULL)
      }
      shiny::tagList(
        shiny::actionButton(session$ns("send"), "Send to upload",
          icon = shiny::icon("right-to-bracket"),
          class = "btn-outline-primary mt-2", width = "100%"
        ),
        shiny::div(idat_info("send"), style = "text-align: center; margin-top: 4px;")
      )
    })

    shiny::observeEvent(input$send, {
      res <- shiny::req(result()$res)

      ## Set the channel and nothing else - the same thing the Library's
      ## recompute button does (loading_table_datasets.R, where the two
      ## navigation lines next to it are commented out on purpose). The board
      ## answers by switching to Upload and showing its wizard itself; a
      ## nav_select from here lands a second panel render on top of the wizard
      ## and it never appears.
      payload <- list(
        counts = res$beta,
        samples = res$samples,
        contrast = NULL,
        organism = "Human", ## the only organism the 450K/EPIC manifests cover
        datatype = "methylomics",
        meth_type = res$array,
        name = "",
        description = ""
      )

      ## Handed over from a later flush, not this button's. The Library's
      ## recompute button reaches the same channel from inside a shinyalert
      ## callbackR, which is a client-originated flush of its own - and that
      ## is the path where the board's wizard actually opens.
      shinyjs::delay(300, recompute_pgx(payload))

      shiny::showNotification(
        sprintf("Sent %d samples to the upload wizard.", ncol(res$beta)),
        type = "message"
      )
    })
  })
}
