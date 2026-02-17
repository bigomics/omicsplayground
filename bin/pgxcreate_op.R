#!/usr/bin/env Rscript

##
## Make PGX file from CSV files
##
## (c) 2023 BigOmics Analytics
##

message("[create PGX process] : starting process")
args <- commandArgs(trailingOnly = TRUE)

temp_dir <- args[1]
OPG <- args[2]
if (!exists("temp_dir")) temp_dir <- getwd()

params_from_op <- file.path(temp_dir, "params.RData")
if (file.exists(params_from_op)) {
  params <- readRDS(params_from_op)
} else {
  yaml::yaml.load_file(file.path(temp_dir, "PARAMS.yml"))
}

# Call create_pgx function
pgx <- playbase::pgx.createPGX(
  organism = params$organism,
  counts = params$counts,
  X = params$countsX,
  norm_method = params$norm_method,
  samples = params$samples,
  contrasts = params$contrasts,
  azimuth_ref = params$azimuth_ref,
  dotimeseries = params$dotimeseries,
  name = params$name,
  datatype = params$datatype,
  probe_type = params$probe_type,
  description = params$description,
  metadata = params$metadata,
  creator = params$creator,
  batch.correct.method = params$batch.correct.method,
  batch.pars = params$batch.pars,
  covariates = params$covariates, ## new
  dma = params$dma, ## new
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

message("[create PGX process] : PGX created successfully")
message("[compute PGX process] : starting process")

pgx <- playbase::pgx.computePGX(
  pgx = pgx,
  max.genes = params$max.genes,
  gx.methods = params$gx.methods,
  gset.methods = params$gset.methods,
  custom.geneset = pgx$custom.geneset,
  custom_fc = params$custom_fc,
  extra.methods = params$extra.methods,
  use.design = params$use.design,
  prune.samples = params$prune.samples,
  do.clustergenes = params$do.cluster,
  do.clustergenesets = params$do.cluster,
  cluster.contrasts = params$cluster.contrasts,
  pgx.dir = params$pgx.save.folder,
  libx.dir = params$libx.dir,
  user_input_dir = temp_dir
)

# embed opg version
pgx$versions$omicsplayground_version <- scan(file.path(OPG, "VERSION"), character())[1]

# annotate pgx

message("[ComputePgxServer:@compute] initialize object\n")

# Save output to a PGX file
pgx_name <- paste0(params$name, ".pgx")
# if pgx.save folder exists, save pgx file to it, otherwise save in temp_dir
if (dir.exists(params$pgx.save.folder)) {
  pgx_name <- file.path(params$pgx.save.folder, pgx_name)
} else {
  pgx_name <- file.path(temp_dir, pgx_name)
}
save(pgx, file = pgx_name)

# We now send success message to user from the shiny app, because the process only runs if the app is alive, so there is no need to send the message here.

# ds_name <- paste0("<b>", params$name, "</b>")
# gmail_creds <- file.path(params$ETC, "gmail_creds")

# params$sendSuccessMessageToUser(
#   user_email = params$email,
#   pgx_name = ds_name,
#   path_to_creds = gmail_creds
# )

message("[compute PGX process] : process finished, pgx is saved as", pgx_name, "\n")
