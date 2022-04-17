## Generated automatically: do not edit by hand

## find package root folder
PKG <<- pkgload::pkg_path()

message('package root = ',PKG)

source(file.path(PKG,'R/app_config.R'),local=TRUE)
source(file.path(PKG,'R/auth.R'),local=TRUE)
source(file.path(PKG,'R/compute2-extra.R'),local=TRUE)
source(file.path(PKG,'R/compute2-genes.R'),local=TRUE)
source(file.path(PKG,'R/compute2-genesets.R'),local=TRUE)
source(file.path(PKG,'R/ggplot-theme.R'),local=TRUE)
source(file.path(PKG,'R/gset-fisher.r'),local=TRUE)
source(file.path(PKG,'R/gset-gsea.r'),local=TRUE)
source(file.path(PKG,'R/gset-meta.r'),local=TRUE)
source(file.path(PKG,'R/gset-rankcor.r'),local=TRUE)
source(file.path(PKG,'R/gx-combat.r'),local=TRUE)
source(file.path(PKG,'R/gx-heatmap.r'),local=TRUE)
source(file.path(PKG,'R/gx-limma.r'),local=TRUE)
source(file.path(PKG,'R/gx-plot.r'),local=TRUE)
source(file.path(PKG,'R/gx-util.r'),local=TRUE)
source(file.path(PKG,'R/gx-volcano.r'),local=TRUE)
source(file.path(PKG,'R/ngs-cook.r'),local=TRUE)
source(file.path(PKG,'R/ngs-fit.r'),local=TRUE)
source(file.path(PKG,'R/ngs-functions.R'),local=TRUE)
source(file.path(PKG,'R/ngs-hiveplot.R'),local=TRUE)
source(file.path(PKG,'R/ngs-salmon.r'),local=TRUE)
source(file.path(PKG,'R/pgx-api.R'),local=TRUE)
source(file.path(PKG,'R/pgx-cluster.R'),local=TRUE)
source(file.path(PKG,'R/pgx-cna.R'),local=TRUE)
source(file.path(PKG,'R/pgx-compute.R'),local=TRUE)
source(file.path(PKG,'R/pgx-contrasts.R'),local=TRUE)
source(file.path(PKG,'R/pgx-correct.R'),local=TRUE)
source(file.path(PKG,'R/pgx-correlation.R'),local=TRUE)
source(file.path(PKG,'R/pgx-deconv.R'),local=TRUE)
source(file.path(PKG,'R/pgx-drugs.R'),local=TRUE)
source(file.path(PKG,'R/pgx-files.R'),local=TRUE)
source(file.path(PKG,'R/pgx-functions.R'),local=TRUE)
source(file.path(PKG,'R/pgx-getgeo.R'),local=TRUE)
source(file.path(PKG,'R/pgx-graph.R'),local=TRUE)
source(file.path(PKG,'R/pgx-grep2.R'),local=TRUE)
source(file.path(PKG,'R/pgx-init.R'),local=TRUE)
source(file.path(PKG,'R/pgx-links.R'),local=TRUE)
source(file.path(PKG,'R/pgx-modules.R'),local=TRUE)
source(file.path(PKG,'R/pgx-plotting.R'),local=TRUE)
source(file.path(PKG,'R/pgx-predict.R'),local=TRUE)
source(file.path(PKG,'R/pgx-proteomics.R'),local=TRUE)
source(file.path(PKG,'R/pgx-pubmed.R'),local=TRUE)
source(file.path(PKG,'R/pgx-purifyExpression.R'),local=TRUE)
source(file.path(PKG,'R/pgx-signature.R'),local=TRUE)
source(file.path(PKG,'R/pgx-singlecell.R'),local=TRUE)
source(file.path(PKG,'R/pgx-tcga.R'),local=TRUE)
source(file.path(PKG,'R/pgx-ui.R'),local=TRUE)
source(file.path(PKG,'R/pgx-vizpanels.R'),local=TRUE)
source(file.path(PKG,'R/pgx-wordcloud.R'),local=TRUE)
source(file.path(PKG,'R/setdirs.R'),local=TRUE)
source(file.path(PKG,'R/ui-code.R'),local=TRUE)
source(file.path(PKG,'R/xcr-graph.r'),local=TRUE)
source(file.path(PKG,'R/xcr-math.r'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/biomarker_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/clustering_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/compare_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/connectivity_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/correlation_server.R'),local=TRUE)

#' DataView module server function
#'
#' @description A shiny Module (server code).
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @param inputData Reactive expression that provides the input ngs/pgx data object 
#'
#' @export 
DataViewBoard <- function(){}
source(file.path(PKG,'shiny/boards/server/dataview_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/drugconnectivity_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/enrichment_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/expression_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/featuremap_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/functional_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/intersection_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/loading_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/signature_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/singlecell_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/tcga_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/user_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/wgcna_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/server/wordcloud_server.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/biomarker_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/clustering_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/compare_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/connectivity_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/correlation_ui.R'),local=TRUE)

#' DataView module UI Input function
#'
#' @description A shiny Module. Renders the input parts (sidebar contents) for the module.
#'
#' @param id Internal parameters for {shiny}.
#' #'
#' @export 
DataViewInputs <- function(){}
#' DataView module UI output function
#'
#' @description Renders the output part for the module as tabsetPanel object
#'
#' @param id Internal parameters for {shiny}.
#' #'
#' @export 
DataViewUI <- function(){}
source(file.path(PKG,'shiny/boards/ui/dataview_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/drugconnectivity_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/enrichment_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/expression_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/featuremap_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/functional_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/intersection_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/loading_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/signature_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/singlecell_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/tcga_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/user_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/wgcna_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/boards/ui/wordcloud_ui.R'),local=TRUE)
source(file.path(PKG,'shiny/modules/AuthenticationModule.R'),local=TRUE)
source(file.path(PKG,'shiny/modules/BatchCorrectModule.R'),local=TRUE)
source(file.path(PKG,'shiny/modules/ComputePgxModule.R'),local=TRUE)
source(file.path(PKG,'shiny/modules/MakeContrastModule.R'),local=TRUE)
source(file.path(PKG,'shiny/modules/NormalizeCountsModule.R'),local=TRUE)
source(file.path(PKG,'shiny/modules/PlotModule.R'),local=TRUE)
source(file.path(PKG,'shiny/modules/plotModules/dataviewTSNEPlotModule.R'),local=TRUE)
source(file.path(PKG,'shiny/modules/plotModules/examplePlotModule.R'),local=TRUE)
source(file.path(PKG,'shiny/modules/plotModules/PlotModule.R'),local=TRUE)
source(file.path(PKG,'shiny/modules/QuestionModule.R'),local=TRUE)
source(file.path(PKG,'shiny/modules/TimerModule.R'),local=TRUE)
source(file.path(PKG,'shiny/modules/UploadModule.R'),local=TRUE)

#' The main application Server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
app_server <- function(){}
source(file.path(PKG,'shiny/server.R'),local=TRUE)

#' The main application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @export
app_ui <- function(){}
source(file.path(PKG,'shiny/ui.R'),local=TRUE)
source(file.path(PKG,'shiny/utils/utils.R'),local=TRUE)
