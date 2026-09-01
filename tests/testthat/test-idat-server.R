## test-idat-server.R

.board_dir <- if (dir.exists("components/app_idat/R")) {
  "components/app_idat/R"
} else {
  "../../components/app_idat/R"
}

source(file.path(.board_dir, "idat_server.R"), local = TRUE)

## Fake what shiny::fileInput hands the server: original name + temp datapath.
.upload <- function(names, dir = tempfile("src")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(dir, paste0("upload", seq_along(names)))
  for (p in paths) writeLines("x", p)
  data.frame(name = names, datapath = paths, stringsAsFactors = FALSE)
}

## IDAT pairs on disk. idat_samples() reads a sheet the way read_idats() did -
## resolved against the files that actually exist - so a sheet fixture only
## means anything when its files are there.
.idats <- function(prefixes, dir = tempfile("stage")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (p in prefixes) {
    for (ch in c("Grn", "Red")) {
      writeLines("x", file.path(dir, paste0(p, "_", ch, ".idat")))
    }
  }
  dir
}

## A ZIP holding entry names verbatim. utils::zip() shells out to the zip
## binary, which strips "../" on the way in, so a traversing archive cannot be
## built with it - and that is the archive the staging code has to survive.
## Stored, zero-length members: CRC and both sizes are then 0.
.raw_zip <- function(names, file = tempfile(fileext = ".zip")) {
  u16 <- function(x) writeBin(as.integer(x), raw(), size = 2, endian = "little")
  u32 <- function(x) writeBin(as.integer(x), raw(), size = 4, endian = "little")
  local <- central <- raw(0)
  for (nm in names) {
    n <- charToRaw(nm)
    central <- c(
      central, charToRaw("PK\001\002"), u16(20), u16(10), u16(0), u16(0),
      u16(0), u16(0), u32(0), u32(0), u32(0), u16(length(n)), u16(0),
      u16(0), u16(0), u16(0), u32(0), u32(length(local)), n
    )
    local <- c(
      local, charToRaw("PK\003\004"), u16(10), u16(0), u16(0), u16(0),
      u16(0), u32(0), u32(0), u32(0), u16(length(n)), u16(0), n
    )
  }
  writeBin(c(
    local, central, charToRaw("PK\005\006"), u16(0), u16(0),
    u16(length(names)), u16(length(names)), u32(length(central)),
    u32(length(local)), u16(0)
  ), file)
  file
}

## ---- staging ---------------------------------------------------------------

test_that("no upload yields an empty inventory with the right columns", {
  st <- idat_stage_upload(NULL)
  expect_equal(nrow(st$inventory), 0)
  expect_named(st$inventory, c("sample", "grn", "red", "complete"))
  expect_true(is.na(st$sheet))
})

test_that("Grn/Red pairs are grouped per sample and flagged complete", {
  files <- .upload(c(
    "200123_R01C01_Grn.idat", "200123_R01C01_Red.idat",
    "200123_R02C01_Grn.idat", "200123_R02C01_Red.idat"
  ))

  inv <- idat_stage_upload(files, exdir = tempfile("out"))$inventory
  expect_equal(inv$sample, c("200123_R01C01", "200123_R02C01"))
  expect_true(all(inv$complete))
  expect_true(all(grepl("_Grn\\.idat$", inv$grn)))
  expect_true(all(grepl("_Red\\.idat$", inv$red)))
})

test_that("a sample missing one channel is kept but marked incomplete", {
  files <- .upload(c("S1_Grn.idat", "S1_Red.idat", "S2_Grn.idat"))

  inv <- idat_stage_upload(files, exdir = tempfile("out"))$inventory
  expect_equal(inv$complete, c(TRUE, FALSE))
  expect_true(is.na(inv$red[inv$sample == "S2"]))
})

test_that("a ZIP upload is unpacked, subdirectories and all", {
  ## GEO ships <Slide>/<Slide>_<Array>_Grn.idat; the structure must survive.
  src <- tempfile("src")
  dir.create(file.path(src, "5723646052"), recursive = TRUE)
  for (f in c("5723646052_R02C02_Grn.idat", "5723646052_R02C02_Red.idat")) {
    writeLines("x", file.path(src, "5723646052", f))
  }
  zip.file <- file.path(src, "idats.zip")
  old <- setwd(src)
  utils::zip(zip.file, "5723646052", flags = "-qr")
  setwd(old)
  skip_if_not(file.exists(zip.file), "zip binary not available")

  files <- data.frame(name = "idats.zip", datapath = zip.file, stringsAsFactors = FALSE)
  st <- idat_stage_upload(files, exdir = tempfile("out"))
  expect_equal(st$inventory$sample, "5723646052_R02C02")
  expect_true(st$inventory$complete)
  expect_true(grepl("5723646052/5723646052_R02C02_Grn", st$inventory$grn))
})

test_that("non-IDAT uploads contribute no samples", {
  files <- .upload(c("readme.txt", "notes.md"))
  expect_equal(nrow(idat_stage_upload(files, exdir = tempfile("out"))$inventory), 0)
})

test_that("a path-traversing upload name cannot escape the staging directory", {
  ## A crafted filename must land inside exdir, never above it.
  exdir <- tempfile("out")
  files <- .upload(c("../../evil_Grn.idat", "../../evil_Red.idat"))
  st <- idat_stage_upload(files, exdir = exdir)
  expect_true(all(startsWith(normalizePath(st$inventory$grn), normalizePath(exdir))))
})

test_that("a traversing ZIP entry is refused and named, not written above exdir", {
  ## unzip() honours "../" in an entry name and overwrites what it lands on, so
  ## a crafted archive could write anywhere the app user can. The good entries
  ## in the same archive must still stage.
  root <- tempfile("root")
  exdir <- file.path(root, "a", "b", "exdir")
  zip.file <- .raw_zip(c(
    "../../ESCAPED_Grn.idat", "../../ESCAPED_Red.idat",
    "ok_Grn.idat", "ok_Red.idat"
  ))
  files <- data.frame(name = "idats.zip", datapath = zip.file, stringsAsFactors = FALSE)

  st <- idat_stage_upload(files, exdir = exdir)
  expect_equal(basename(list.files(root, recursive = TRUE)), c("ok_Grn.idat", "ok_Red.idat"))
  expect_equal(st$inventory$sample, "ok")
  ## Refused, not silently dropped: the summary says which entries went.
  expect_equal(st$rejected, c("../../ESCAPED_Grn.idat", "../../ESCAPED_Red.idat"))
})

test_that("the same pair name in two subdirectories stays two samples", {
  ## Per-plate subdirectories routinely repeat a position. Keying on the
  ## basename folded them into one row that kept plate1 and lost plate2 -
  ## half the cohort gone, reported as "1 sample ready".
  src <- tempfile("src")
  for (p in c("plate1", "plate2")) {
    dir.create(file.path(src, p), recursive = TRUE)
    for (ch in c("Grn", "Red")) {
      writeLines("x", file.path(src, p, sprintf("A01_%s.idat", ch)))
    }
  }
  inv <- idat_inventory(src)
  expect_equal(nrow(inv), 2)
  expect_true(all(inv$complete))
  expect_equal(inv$sample, c("A01", "A01.1"))
  expect_true(grepl("plate1/", inv$grn[1]) && grepl("plate2/", inv$grn[2]))
})

## ---- sheet detection and synthesis -----------------------------------------

test_that("a sample-sheet CSV in the upload is preferred over another CSV", {
  files <- .upload(c("S1_Grn.idat", "S1_Red.idat", "readme.csv", "SampleSheet.csv"))
  st <- idat_stage_upload(files, exdir = tempfile("out"))
  expect_equal(basename(st$sheet), "SampleSheet.csv")
  ## The full path, not the basename: read_idats() resolves relative to the
  ## sheet, so a bare filename fails with "neither a directory nor a CSV file".
  expect_true(file.exists(st$sheet))
})

test_that("a lone unnamed CSV is still taken as the sheet", {
  files <- .upload(c("S1_Grn.idat", "S1_Red.idat", "targets.csv"))
  st <- idat_stage_upload(files, exdir = tempfile("out"))
  expect_equal(basename(st$sheet), "targets.csv")
  expect_true(file.exists(st$sheet))
})

test_that("the synthesized sheet carries an absolute Basename and no Slide/Array", {
  ## minfi::read.metharray.sheet() rebuilds Basename by grepping
  ## "<Slide>_<Array>_Grn.idat" whenever an Array column exists. Omitting
  ## Slide/Array is what makes our absolute path survive intact.
  files <- .upload(c("weird.name_Grn.idat", "weird.name_Red.idat"))
  st <- idat_stage_upload(files, exdir = tempfile("out"))
  sheet <- idat_write_sheet(st$inventory, st$dir)

  df <- read.csv(sheet, stringsAsFactors = FALSE)
  expect_named(df, c("Sample_Name", "Basename"))
  expect_false(any(c("Slide", "Array") %in% colnames(df)))
  expect_equal(df$Sample_Name, "weird.name")
  expect_true(startsWith(df$Basename, "/"))
  expect_true(file.exists(paste0(df$Basename, "_Grn.idat")))
  expect_true(file.exists(paste0(df$Basename, "_Red.idat")))
})

test_that("only complete pairs reach the synthesized sheet", {
  files <- .upload(c("S1_Grn.idat", "S1_Red.idat", "S2_Grn.idat"))
  st <- idat_stage_upload(files, exdir = tempfile("out"))
  df <- read.csv(idat_write_sheet(st$inventory, st$dir), stringsAsFactors = FALSE)
  expect_equal(df$Sample_Name, "S1")
})

test_that("writing a sheet with nothing complete is an error, not an empty sheet", {
  files <- .upload("S2_Grn.idat")
  st <- idat_stage_upload(files, exdir = tempfile("out"))
  expect_error(idat_write_sheet(st$inventory, st$dir), "no complete IDAT pairs")
})

## ---- conversion wrapper ----------------------------------------------------

test_that("a read_idats() failure comes back as data, not a thrown error", {
  ## The future worker cannot raise: an error has to survive as list(ok, msg).
  out <- idat_convert(tempfile("nope"))
  expect_false(out$ok)
  expect_true(nzchar(out$msg))
  expect_null(out$res)
})

test_that("the log filter keeps tagged messages and problems, drops load noise", {
  ## A worker loads the whole Bioc stack on first use; without this the Log tab
  ## is ~110 lines of banner around the ~15 that say anything.
  x <- c(
    "[playbase.epigenetics::read_idats] Array type: 450K array",
    "Loading required package: minfi",
    "Attaching package: 'BiocGenerics'",
    "    IQR, mad, sd, var, xtabs",
    "",
    "WARNING: Low matching ratio",
    "Error in `$<-.data.frame`: boom"
  )
  expect_equal(idat_log_lines(x), x[c(1, 6, 7)])
})

test_that("the log filter drops the two expected-but-alarming warnings", {
  ## minfi complains about the Slide/Array columns idat_write_sheet omits on
  ## purpose, and playbase reports its own import masking. Neither is a problem.
  x <- c(
    "Warning: Could not infer array name for file: /tmp/x/_playground_sheet.csv",
    "Warning: Could not infer slide name for file: /tmp/x/_playground_sheet.csv",
    "Warning: replacing previous import 'data.table::first' by 'dplyr::first'",
    "Warning: something that actually matters"
  )
  expect_equal(idat_log_lines(x), x[4])
})

## ---- progress feed ---------------------------------------------------------

test_that("progress starts at the load phase when the file does not exist yet", {
  ## The worker's slowest silent stretch is loading Bioconductor, before it
  ## has written a single line - that must not read as a stalled bar at 0.
  p <- idat_progress(tempfile("never"), elapsed = 3)
  expect_equal(p$pct, 2)
  expect_match(p$label, "Starting up")
  expect_equal(p$elapsed, 3)
})

test_that("progress tracks the furthest phase the log has reached", {
  f <- tempfile("prog")
  writeLines(c(
    "[playbase.epigenetics::read_idat_targets] Parsed 6 sample row(s) from the sheet.",
    "[playbase.epigenetics::read_idats] Reading 6 IDAT pair(s)...",
    "[playbase.epigenetics::read_idats] Raw QC (median intensities, bisulfite conversion)..."
  ), f)
  p <- idat_progress(f)
  expect_equal(p$pct, IDAT_PHASES$pct[3])
  expect_match(p$label, "Raw intensity QC")
})

test_that("progress never goes backwards when a later phase logs more lines", {
  ## mask_by_detection can emit two lines, one line or none depending on the
  ## cohort; the bar tracks the furthest phase seen, not the last line.
  f <- tempfile("prog")
  writeLines(c(
    "[playbase.epigenetics::read_idats] Annotating 467485 probes...",
    "[playbase.epigenetics::mask_by_detection] dropping 18027 probe(s)"
  ), f)
  expect_equal(idat_progress(f)$pct, IDAT_PHASES$pct[6])
})

test_that("every phase pattern matches a real read_idats() message", {
  ## Pinned against the actual message text: a reworded message upstream
  ## silently freezes the bar, and this is what catches it.
  real <- c(
    "[playbase.epigenetics::read_idat_targets] Parsed 6 sample row(s) from the sheet.",
    "[playbase.epigenetics::read_idats] Reading 6 IDAT pair(s)...",
    "[playbase.epigenetics::read_idats] Raw QC (median intensities, bisulfite conversion)...",
    "[playbase.epigenetics::read_idats] Detection p-values...",
    "[playbase.epigenetics::mask_by_detection] dropping 18027 probe(s)",
    "[playbase.epigenetics::read_idats] Annotating 467485 probes...",
    "[playbase.epigenetics::annotate_methylomics] Annotation completed"
  )
  for (i in seq_len(nrow(IDAT_PHASES))) {
    expect_true(
      any(grepl(IDAT_PHASES$pattern[i], real)),
      info = paste("phase", i, IDAT_PHASES$label[i], "matches nothing")
    )
  }
  expect_true(all(diff(IDAT_PHASES$pct) > 0))
})

## ---- array type and sample table -------------------------------------------

test_that("the detected array type is read back out of the log", {
  ## Pinned against read_idats()'s actual message: get this wrong and the
  ## upload board silently annotates EPIC data against the 450K manifest.
  log <- c(
    "[playbase.epigenetics::read_idats] Reading 6 IDAT pair(s)...",
    "[playbase.epigenetics::read_idats] Array type: EPIC array; method: noob",
    "[playbase.epigenetics::read_idats] Detection p-values..."
  )
  expect_equal(idat_array_type(log), "EPIC array")
  expect_true(is.na(idat_array_type("nothing here")))
})

test_that("with no sheet the sample table is the ids, named by themselves", {
  res <- list(beta = matrix(0, 2, 2, dimnames = list(NULL, c("S1", "S2"))))
  s <- idat_samples(res, NA_character_)
  expect_equal(rownames(s), c("S1", "S2"))
  expect_equal(s$sample, c("S1", "S2"))
})

test_that("a sheet's phenotype columns are carried over, aligned to the betas", {
  ## Note the sheet's row order differs from the beta columns on purpose.
  d <- .idats(c("a", "b"))
  f <- file.path(d, "SampleSheet.csv")
  write.csv(data.frame(
    Sample_Name = c("S2", "S1"), Slide = c("x", "y"),
    Basename = c("a", "b"), group = c("case", "ctrl"),
    age = c(61, 44), stringsAsFactors = FALSE
  ), f, row.names = FALSE)

  res <- list(beta = matrix(0, 2, 2, dimnames = list(NULL, c("S1", "S2"))))
  s <- idat_samples(res, f)
  expect_equal(rownames(s), c("S1", "S2"))
  expect_equal(s$group, c("ctrl", "case"))
  expect_equal(s$age, c(44, 61))
  ## Basename locates a file on one machine; Sample_Name is already the row
  ## names. Slide stays - chip is a batch covariate, not bookkeeping.
  expect_false(any(c("Basename", "Sample_Name") %in% colnames(s)))
  expect_true("Slide" %in% colnames(s))
})

test_that("an Illumina [Header]/[Data] sheet is read from the right row", {
  ## Real Illumina sheets open with a [Header] block, and a reader that takes
  ## "[Header]" as the column names finds no Sample_Name and drops every
  ## phenotype the sheet carried - silently.
  d <- .idats(c("5723646052_R01C01", "5723646052_R02C01"))
  f <- file.path(d, "SampleSheet.csv")
  writeLines(c(
    "[Header],,,,",
    "Investigator Name,MrNoName,,,",
    ",,,,",
    "[Data],,,,",
    "Sample_Name,Sentrix_ID,Sentrix_Position,age,status",
    "S1,5723646052,R01C01,58,normal",
    "S2,5723646052,R02C01,75,cancer"
  ), f)

  res <- list(beta = matrix(0, 2, 2, dimnames = list(NULL, c("S1", "S2"))))
  s <- idat_samples(res, f)
  expect_equal(rownames(s), c("S1", "S2"))
  expect_equal(s$age, c(58, 75))
  expect_equal(s$status, c("normal", "cancer"))
  ## Sentrix_ID/Sentrix_Position come back as Slide/Array: the reader
  ## normalises the naming variants, and the chip stays as a batch covariate.
  expect_true(all(c("Slide", "Array") %in% colnames(s)))
})

test_that("a duplicated Sample_Name attaches to the row that was actually read", {
  ## Technical replicates share a Sample_Name. Uniquifying over every CSV row -
  ## including the one whose IDATs were never staged - names the dropped row
  ## "S1" and the survivor "S1.1", so the survivor's phenotypes come back as
  ## the dropped row's. Silently wrong ages, and a wrong Slide batch covariate.
  d <- .idats(c("S1b", "S2")) ## note: no S1a on disk
  f <- file.path(d, "SampleSheet.csv")
  write.csv(data.frame(
    Sample_Name = c("S1", "S1", "S2"),
    Basename = c("S1a", "S1b", "S2"),
    Slide = c("plate1", "plate2", "plate2"),
    group = c("never_read", "survivor", "ctrl"),
    age = c(99, 44, 61), stringsAsFactors = FALSE
  ), f, row.names = FALSE)

  res <- list(beta = matrix(0, 2, 2, dimnames = list(NULL, c("S1", "S2"))))
  s <- idat_samples(res, f)
  expect_equal(rownames(s), c("S1", "S2"))
  expect_equal(s$group, c("survivor", "ctrl"))
  expect_equal(s$age, c(44, 61))
  expect_equal(s$Slide, c("plate2", "plate2"))
})

test_that("a sheet naming only some of the samples is refused, not half-joined", {
  ## One match used to accept the whole sheet, and every other sample came
  ## back all-NA - phenotypes that look filled in and are not.
  d <- .idats(c("S1", "S2"))
  f <- file.path(d, "SampleSheet.csv")
  write.csv(data.frame(
    Sample_Name = "S1", Basename = "S1", group = "case",
    stringsAsFactors = FALSE
  ), f, row.names = FALSE)

  res <- list(beta = matrix(0, 2, 2, dimnames = list(NULL, c("S1", "S2"))))
  expect_equal(colnames(idat_samples(res, f)), "sample")
})

test_that("a sheet that cannot be joined falls back rather than misaligning", {
  ## No Sample_Name means no key; guessing the order would silently attach the
  ## wrong phenotype to the wrong sample.
  f <- tempfile(fileext = ".csv")
  write.csv(data.frame(Slide = c("a", "b"), group = c("case", "ctrl")), f,
    row.names = FALSE)
  res <- list(beta = matrix(0, 2, 2, dimnames = list(NULL, c("S1", "S2"))))
  expect_equal(colnames(idat_samples(res, f)), "sample")

  ## A sheet naming entirely different samples is no better a join.
  f2 <- tempfile(fileext = ".csv")
  write.csv(data.frame(Sample_Name = c("X9", "X8"), group = c("case", "ctrl")),
    f2, row.names = FALSE)
  expect_equal(colnames(idat_samples(res, f2)), "sample")
})

## ---- reconciling the inventory against the conversion ----------------------

test_that("staged pairs the sheet never listed are reported, not lost", {
  ## Four complete pairs plus a sheet naming two: the reader returns two
  ## samples and an empty failures frame, so the Log tab used to say "All
  ## samples read successfully" with half the cohort gone.
  inv <- data.frame(
    sample = c("S1", "S2", "S3", "S4"), complete = TRUE,
    stringsAsFactors = FALSE
  )
  res <- list(
    beta = matrix(0, 2, 2, dimnames = list(NULL, c("S1", "S2"))),
    failed = data.frame(sample = character(0), reason = character(0))
  )
  out <- idat_reconcile(res, inv)
  expect_equal(out$sample, c("S3", "S4"))
  expect_match(out$reason[1], "never read")
})

test_that("a shortfall whose names cannot be matched is still reported by count", {
  ## A user sheet names samples whatever it likes, so the inventory's file
  ## prefixes need not be the beta column names. Better a count than a guess.
  inv <- data.frame(
    sample = c("GSM1_R01C01", "GSM2_R01C02"), complete = TRUE,
    stringsAsFactors = FALSE
  )
  res <- list(
    beta = matrix(0, 2, 1, dimnames = list(NULL, "Control_1")),
    failed = NULL
  )
  out <- idat_reconcile(res, inv)
  expect_equal(nrow(out), 1)
  expect_equal(out$sample, "1 sample(s)")
})

test_that("a fully accounted-for conversion adds nothing", {
  ## Incomplete pairs are already reported by the summary, and a sample the
  ## reader named as failed must not be counted twice.
  inv <- data.frame(
    sample = c("S1", "S2", "S3"), complete = c(TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  res <- list(
    beta = matrix(0, 2, 1, dimnames = list(NULL, "S1")),
    failed = data.frame(sample = "S2", reason = "unreadable IDAT",
      stringsAsFactors = FALSE)
  )
  expect_equal(idat_reconcile(res, inv), res$failed)
})

## ---- download bundle -------------------------------------------------------

test_that("the bundle holds counts/samples/qc and round-trips the betas", {
  beta <- matrix(c(0.1, 0.9, 0.5, 0.4), nrow = 2,
    dimnames = list(c("cg001", "cg002"), c("S1", "S2")))
  res <- list(
    beta = beta,
    samples = idat_samples(list(beta = beta)),
    ledger = data.frame(`Median meth` = c(11.2, 11.4), row.names = c("S1", "S2"),
      check.names = FALSE)
  )
  zip.file <- tempfile(fileext = ".zip")
  idat_bundle(res, zip.file)
  skip_if_not(file.exists(zip.file), "zip binary not available")

  expect_setequal(utils::unzip(zip.file, list = TRUE)$Name,
    c("counts.csv", "samples.csv", "qc.csv"))

  out <- tempfile("unz")
  utils::unzip(zip.file, exdir = out)
  back <- as.matrix(read.csv(file.path(out, "counts.csv"), row.names = 1))
  expect_equal(back, beta)
  ## row.names = 1: the upload board matches counts columns against these.
  samples <- read.csv(file.path(out, "samples.csv"), row.names = 1,
    stringsAsFactors = FALSE)
  expect_equal(rownames(samples), c("S1", "S2"))
})
