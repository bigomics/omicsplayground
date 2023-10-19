#!/usr/bin/env Rscript

##
## Make PGX file from CSV files
##
## (c) 2023 BigOmics Analytics
##

message("[compute PGX process] : starting process\n")
args = commandArgs(trailingOnly=TRUE)

temp_dir <- args[1]

if(!exists("temp_dir")) {
  temp_dir <- getwd()
}

params_from_op <- file.path(temp_dir,"params.RData")

if (file.exists(params_from_op)) {
  params <- readRDS(params_from_op)
} else {
  yaml::yaml.load_file(file.path(temp_dir,"PARAMS.yml"))
}

# Call create_pgx function
pgx <- playbase::pgx.createPGX(
  organism = params$organism,
  counts = params$counts,
  samples = params$samples,
  contrasts = params$contrasts,
  name <- params$name,
  datatype <- params$datatype,
  description <- params$description,
  creator <- params$creator,
  X = NULL,
  batch.correct = params$batch.correct,
  prune.samples = params$prune.samples,
  filter.genes = params$filter.genes,
  only.known = params$only.known,
  only.proteincoding = params$only.proteincoding,
  only.hugo = params$only.hugo,
  convert.hugo = params$convert.hugo,
  do.cluster = params$do.cluster,
  cluster.contrasts = params$cluster.contrasts
)

pgx <- playbase::pgx.computePGX(
  pgx = pgx,
  max.genes = params$max.genes,
  max.genesets = params$max.genesets,
  gx.methods = params$gx.methods,
  gset.methods = params$gset.methods,
  custom.geneset = params$custom.geneset,
  extra.methods = params$extra.methods,
  use.design = params$use.design,        ## no.design+prune are combined
  prune.samples = params$prune.samples,  ##
  do.cluster = params$do.cluster,
  pgx.dir = params$pgx.save.folder,
  libx.dir = params$libx.dir
  )

# annotate pgx

message("[ComputePgxServer:@compute] initialize object\n")

# Save output to a PGX file

pgx_name <- paste0(params$name,".pgx")
# if pgx.save folder exists, save pgx file to it, otherwise save in temp_dir
if (dir.exists(params$pgx.save.folder)) {
  pgx_name <- file.path(params$pgx.save.folder,pgx_name)
} else {
  pgx_name <- file.path(temp_dir,pgx_name)
}
save(pgx, file = pgx_name)

message("[compute PGX process] : process finished, pgx is saved as", pgx_name,"\n")
