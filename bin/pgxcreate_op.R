#!/usr/bin/env Rscript

##
## Make PGX file from CSV files
##
## (c) 2023 BigOmics Analytics 
##

message("[compute PGX process] : starting process")
args = commandArgs(trailingOnly=TRUE)

temp_dir <- args[1]

params  <- readRDS(file.path(temp_dir,"params.RData"))

# Call create_pgx function
pgx <- playbase::create_pgx(
  counts = params$counts,
  samples = params$samples,
  contrasts = params$contrasts,
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

pgx <- playbase::compute_pgx(
  pgx = pgx,	
  max.genes = params$max.genes,
  max.genesets = params$max.genesets,
  gx.methods = params$gx.methods,
  gset.methods = params$gset.methods,
  extra.methods = params$extra.methods,
  use.design = params$use.design,        ## no.design+prune are combined
  prune.samples = params$prune.samples,  ##
  do.cluster = params$do.cluster,
  )

# annotate pgx

pgx$name <- params$name
pgx$datatype <- params$datatype
pgx$description <- params$description
pgx$creator <- params$creator
pgx$date <- params$date

message("[ComputePgxServer:@compute] initialize object")

# Save output to a CSV file
save(pgx, file = file.path(temp_dir,"my.pgx"))

message("[compute PGX process] : process finished, pgx is saved to my.pgx")