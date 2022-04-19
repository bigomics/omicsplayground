## Generated automatically: do not edit by hand

## find package root folder
PKG <<- pkgload::pkg_path()

message('package root = ',PKG)

source(file.path(PKG,'R/app_config.R'),encoding='UTF-8')
source(file.path(PKG,'R/auth.R'),encoding='UTF-8')
source(file.path(PKG,'R/compute2-extra.R'),encoding='UTF-8')
source(file.path(PKG,'R/compute2-genes.R'),encoding='UTF-8')
source(file.path(PKG,'R/compute2-genesets.R'),encoding='UTF-8')
source(file.path(PKG,'R/ggplot-theme.R'),encoding='UTF-8')
source(file.path(PKG,'R/gset-fisher.r'),encoding='UTF-8')
source(file.path(PKG,'R/gset-gsea.r'),encoding='UTF-8')
source(file.path(PKG,'R/gset-meta.r'),encoding='UTF-8')
source(file.path(PKG,'R/gset-rankcor.r'),encoding='UTF-8')
source(file.path(PKG,'R/gx-combat.r'),encoding='UTF-8')
source(file.path(PKG,'R/gx-heatmap.r'),encoding='UTF-8')
source(file.path(PKG,'R/gx-limma.r'),encoding='UTF-8')
source(file.path(PKG,'R/gx-plot.r'),encoding='UTF-8')
source(file.path(PKG,'R/gx-util.r'),encoding='UTF-8')
source(file.path(PKG,'R/gx-volcano.r'),encoding='UTF-8')
source(file.path(PKG,'R/ngs-cook.r'),encoding='UTF-8')
source(file.path(PKG,'R/ngs-fit.r'),encoding='UTF-8')
source(file.path(PKG,'R/ngs-functions.R'),encoding='UTF-8')
source(file.path(PKG,'R/ngs-hiveplot.R'),encoding='UTF-8')
source(file.path(PKG,'R/ngs-salmon.r'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-api.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-cluster.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-cna.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-compute.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-contrasts.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-correct.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-correlation.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-deconv.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-drugs.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-files.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-functions.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-getgeo.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-graph.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-grep2.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-init.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-links.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-modules.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-plotting.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-predict.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-proteomics.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-pubmed.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-purifyExpression.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-signature.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-singlecell.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-tcga.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-ui.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-vizpanels.R'),encoding='UTF-8')
source(file.path(PKG,'R/pgx-wordcloud.R'),encoding='UTF-8')
source(file.path(PKG,'R/setdirs.R'),encoding='UTF-8')
source(file.path(PKG,'R/ui-code.R'),encoding='UTF-8')
source(file.path(PKG,'R/xcr-graph.r'),encoding='UTF-8')
source(file.path(PKG,'R/xcr-math.r'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/biomarker_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/clustering_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/compare_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/connectivity_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/correlation_server.R'),encoding='UTF-8')

#' DataView module server function
#'
#' @description A shiny Module (server code).
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @param inputData Reactive expression that provides the input ngs/pgx data object 
#'
#' @export 
DataViewBoard <- function(){}
source(file.path(PKG,'shiny/boards/server/dataview_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/drugconnectivity_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/enrichment_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/expression_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/featuremap_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/functional_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/intersection_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/loading_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/signature_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/singlecell_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/tcga_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/user_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/wgcna_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/server/wordcloud_server.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/biomarker_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/clustering_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/compare_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/connectivity_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/correlation_ui.R'),encoding='UTF-8')

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
source(file.path(PKG,'shiny/boards/ui/dataview_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/drugconnectivity_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/enrichment_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/expression_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/featuremap_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/functional_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/intersection_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/loading_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/signature_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/singlecell_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/tcga_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/user_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/wgcna_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/boards/ui/wordcloud_ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/modules/AuthenticationModule.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/modules/BatchCorrectModule.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/modules/ComputePgxModule.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/modules/MakeContrastModule.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/modules/NormalizeCountsModule.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/modules/PlotModule.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/modules/plotModules/dataviewTSNEPlotModule.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/modules/plotModules/examplePlotModule.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/modules/plotModules/PlotModule.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/modules/QuestionModule.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/modules/TimerModule.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/modules/UploadModule.R'),encoding='UTF-8')

#' The main application Server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
app_server <- function(){}
source(file.path(PKG,'shiny/server.R'),encoding='UTF-8')

#' The main application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @export
app_ui <- function(){}
source(file.path(PKG,'shiny/ui.R'),encoding='UTF-8')
source(file.path(PKG,'shiny/utils/utils.R'),encoding='UTF-8')
