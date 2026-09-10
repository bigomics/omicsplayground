#!/usr/bin/env Rscript
##
## Execute a params.RData, seeded. Used for BOTH paths -- the caller decides
## which playbase and which omicsplayground tree via R_LIBS_USER and <opg-root>:
##
##   A: R_LIBS_USER=libs/pb-main Rscript run_pgx.R runs/<sc>/A <master-worktree> create
##   B: R_LIBS_USER=libs/pb-edgy Rscript run_pgx.R runs/<sc>/B <edgy-worktree>   create
##
## stage=create -> pgx.createPGX only, saved as createpgx.rds        (gate 2)
## stage=full   -> the real bin/pgxcreate_op.R, producing <name>.pgx (gate 3)
##

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 3) stop("usage: run_pgx.R <params-dir> <opg-root> <create|full> [seed]")
params_dir <- normalizePath(argv[1])
OPG <- normalizePath(argv[2])
stage <- argv[3]
seed <- if (length(argv) >= 4) as.integer(argv[4]) else 42L

set.seed(seed)
suppressMessages(library(playbase))
message("[run_pgx] stage=", stage, " playbase=", as.character(packageVersion("playbase")),
  " opg=", basename(OPG), " seed=", seed)

if (stage == "full") {
  ## bin/pgxcreate_op.R reads its own args from commandArgs(trailingOnly=TRUE),
  ## so re-point them at (params_dir, OPG) and source it in this seeded session.
  local({
    real_commandArgs <- base::commandArgs
    assign("commandArgs", function(trailingOnly = FALSE) {
      if (trailingOnly) c(params_dir, OPG) else real_commandArgs(FALSE)
    }, envir = globalenv())
  })
  source(file.path(OPG, "bin", "pgxcreate_op.R"), local = FALSE)

  ## A captured params still carries the APP's pgx.save.folder (auth$user_dir), so
  ## Path A's .pgx lands in data/ rather than the run dir. Collect it either way so
  ## gate 3 has both files side by side.
  p <- readRDS(file.path(params_dir, "params.RData"))
  want <- paste0(p$name, ".pgx")
  for (d in unique(c(p$pgx.save.folder, params_dir))) {
    src <- file.path(d, want)
    if (file.exists(src) && !file.exists(file.path(params_dir, want))) {
      file.copy(src, file.path(params_dir, want))
      message("[run_pgx] collected ", src, " -> ", params_dir)
    }
  }
  message("[run_pgx] full pipeline done")
  quit(status = 0)
}

if (stage != "create") stop("unknown stage: ", stage)

params <- readRDS(file.path(params_dir, "params.RData"))

## Same argument list bin/pgxcreate_op.R uses, plus `preprocess`. playbase main
## has no `preprocess` formal, so drop anything this version does not accept --
## that keeps ONE runner working against both branches.
args <- list(
  organism = params$organism,
  counts = params$counts,
  X = params$countsX,
  preprocess = params$preprocess,
  norm_method = params$norm_method,
  samples = params$samples,
  contrasts = params$contrasts,
  azimuth_ref = params$azimuth_ref,
  dotimeseries = params$dotimeseries,
  name = params$name,
  datatype = params$datatype,
  datatype_subtype = params$datatype_subtype,
  probe_type = params$probe_type,
  description = params$description,
  metadata = params$metadata,
  creator = params$creator,
  batch.correct.method = params$batch.correct.method,
  batch.pars = params$batch.pars,
  covariates = params$covariates,
  dma = params$dma,
  remove.xy.probes = params$remove.xy.probes,
  meth_type = params$meth_type,
  prune.samples = params$prune.samples,
  filter.genes = params$filter.genes,
  exclude.genes = params$exclude.genes,
  only.known = params$only.known,
  average.duplicated = params$average.duplicated,
  only.proteincoding = params$only.proteincoding,
  only.hugo = params$only.hugo,
  convert.hugo = params$convert.hugo,
  custom.geneset = params$custom.geneset,
  max.genesets = params$max.genesets,
  annot_table = params$annot_table,
  settings = params$settings,
  sc_compute_settings = params$sc_compute_settings
)

known <- names(formals(playbase::pgx.createPGX))
dropped <- setdiff(names(args), known)
if (length(dropped)) {
  message("[run_pgx] this playbase has no formal(s): ", paste(dropped, collapse = ", "),
    " -- dropping")
}
args <- args[intersect(names(args), known)]

if (!is.null(args$preprocess) && !is.null(args$X)) {
  ## createPGX silently prefers X, so a run with both would prove nothing.
  stop("[run_pgx] params has BOTH X (countsX) and preprocess -- preprocess would be ignored")
}

pgx <- do.call(playbase::pgx.createPGX, args)

saveRDS(pgx, file = file.path(params_dir, "createpgx.rds"))
message("[run_pgx] createPGX done: ", nrow(pgx$counts), " x ", ncol(pgx$counts),
  " counts, X ", nrow(pgx$X), " x ", ncol(pgx$X))
message("[run_pgx] wrote ", file.path(params_dir, "createpgx.rds"))
