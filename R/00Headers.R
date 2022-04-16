# Generated automatically: do not edit by hand!

#' Access files in the current app
#'
#' NOTE: If you manually change your package name in the DESCRIPTION,
#' don't forget to change it here too, and in the config file.
#' For a safer name change mechanism, use the `golem::set_golem_name()` function.
#'
#' @param ... character vectors, specifying subdirectory and file(s)
#' within your package. The default, none, returns the root of the app.
#'
#' @noRd
app_sys <- function(){}
#' Read App Config
#'
#' @param value Value to retrieve from the config file.
#' @param config GOLEM_CONFIG_ACTIVE value. If unset, R_CONFIG_ACTIVE.
#' If unset, "default".
#' @param use_parent Logical, scan the parent directory for config file.
#' @param file Location of the config file
#'
#' @noRd
get_golem_config <- function(){}
source('R/app_config.R',local=TRUE)


source('R/auth.R',local=TRUE)


source('R/compute2-extra.R',local=TRUE)


source('R/compute2-genes.R',local=TRUE)


source('R/compute2-genesets.R',local=TRUE)

#' The plot theme to be used for figures in the OmicsPlayground app.
#'
#' @param style (string) Overall color style of text labels.
#' Either "default" or "light".
#' @param base_size (integer) Base point size
#' @param grid (string) Grid lines. Options include  "none" or
#' any combination of "X", "Y", "x" and "y".
#' @param axistitle (string) Axis titles. Options include "none" or
#' any combination of "X", "Y", "x" and "y".
#' @param axistext (string) Axis text labels for values or groups.
#' Options include "none" or any combination of "X", "Y", "x" and "y".
#' @param axis_num (string) Should axis text be formatted as monospaced? Set 
#' to  x and y, respectively, in case numeric values are displayed. Options 
#' include "none" or any combination of "X", "Y", "x" and "y". 
#' @param legend_num (logical) Should legend text be formatted as monospaced?
#' Defaults to FALSE (no monospace font). Set to TRUE in case of numeric values.
#' #' @param panelborder (logical) Should a panel border be drawn?
#' Defaults to FALSE (no border). If TRUE it also adds tick marks to both axes.
#' @param margin (numeric) Should a margin of x be added to the plot?
#' Defaults to 0 (no margin by default).
#' @param ... Other arguments passed to ggplot methods.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mpg, aes(class)) + geom_bar() + theme_omics()
#' ggplot(mpg, aes(class)) + geom_bar() +
#'   theme_omics(style = "light", grid = "xy", margin = 20)
#' ggplot(mpg, aes(class)) + geom_bar() +
#'   theme_omics(grid = "none", axistext = "x", 
#'               axistitle = "none", panelborder = TRUE)
#' }
#'
#' @export
theme_omics <- function(){}
#' A nicely styled legend for categorical aesthetics (color, fill, shape,
#' size and alpha)
#'
#' @param aes (string) Aesthetic of the legend that should be modified.
#' Options include "color", "fill", "shape", "size" and "alpha".
#' @param reverse (boolean) If TRUE legend items are shown in reversed order.
#' Defaults to FALSE.
#' @param ... Other arguments passed to ggplot methods.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mpg, aes(cty, displ)) + geom_point(aes(color = manufacturer)) +
#'   theme_omics() + guide_discrete()
#' ggplot(mpg, aes(class)) + geom_bar(aes(fill = manufacturer)) +
#'   theme_omics() + guide_discrete(aes = "fill")
#' }
#'
#' @export
guide_discrete <- function(){}
#' A nicely styled colorbar for continuous aesthetics (color, fill and alpha)
#'
#' @param aes (string) Aesthetic of the legend that should be modified.
#' Options include "color", "fill", "shape" and "size".
#' @param aes (type) Type of the guide. Options include "bar" and "steps".
#' @param width (number) Width of the color bar.
#' Options include "color", "fill", "shape" and "size".
#' @param ... Other arguments passed to ggplot methods.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mpg, aes(cty, hwy)) + geom_point(aes(color = displ)) +
#'   theme_omics() + guide_continuous()
#' ggplot(mpg, aes(class)) + geom_bar(aes(fill = manufacturer)) +
#'   theme_omics() + guide_continuous(aes = "fill")
#' }
#'
#' @export
guide_continuous <- function(){}
source('R/ggplot-theme.R',local=TRUE)


source('R/gset-fisher.r',local=TRUE)


source('R/gset-gsea.r',local=TRUE)


source('R/gset-meta.r',local=TRUE)


source('R/gset-rankcor.r',local=TRUE)


source('R/gx-combat.r',local=TRUE)


source('R/gx-heatmap.r',local=TRUE)


source('R/gx-limma.r',local=TRUE)


source('R/gx-plot.r',local=TRUE)


source('R/gx-util.r',local=TRUE)


source('R/gx-volcano.r',local=TRUE)


source('R/ngs-cook.r',local=TRUE)


source('R/ngs-fit.r',local=TRUE)


source('R/ngs-functions.R',local=TRUE)


source('R/ngs-hiveplot.R',local=TRUE)


source('R/ngs-salmon.r',local=TRUE)


source('R/pgx-api.R',local=TRUE)


source('R/pgx-cluster.R',local=TRUE)


source('R/pgx-cna.R',local=TRUE)


source('R/pgx-compute.R',local=TRUE)


source('R/pgx-contrasts.R',local=TRUE)


source('R/pgx-correct.R',local=TRUE)


source('R/pgx-correlation.R',local=TRUE)


source('R/pgx-deconv.R',local=TRUE)


source('R/pgx-drugs.R',local=TRUE)


source('R/pgx-files.R',local=TRUE)

#' Given a Set of Points and Box sizes,
#' https://github.com/slowkow/ggrepel/issues/24
util.findboxes <- function(){}
source('R/pgx-functions.R',local=TRUE)


source('R/pgx-getgeo.R',local=TRUE)


source('R/pgx-graph.R',local=TRUE)

#' An R-based wrapper for Trim Galore!
#'
#' @description Run the Trim Galore! tool
#'
#' @details This script runs the Trim Galore! tool and requires installation
#' of both Cutadapt and Trim Galore!  It is essential that Cutadapt is in the
#' executable path otherwise this tool will not work.
#'
#' @param fastq1 a character vector indicating the read files to be trimmed.
#' @param fastq2 (optional) a character vector indicating read files to be
#' trimmmed.  If specified, it is assumed the reads are paired, and this vector
#' MUST be in the same order as those listed in \code{fastq1}.  If \code{NULL}
#' then it is assumed the reads are single-end.
#' @param adapter1 a character string specifying the adapter sequence to be
#' trimmed. If not specified explicitly, Trim Galore will try to auto-detect
#' whether the Illumina universal, Nextera transposase or Illumina small RNA
#' adapter sequence was used. Also see \code{illumina}, \code{nextera} and
#' \code{small_rna} options. If no adapter can be detected within the first 1
#' million sequences of the first file specified Trim Galore defaults to
#' \code{illumina}.
#' @param adapter2 a character string specifying an optional adapter sequence to
#' be trimmed off read 2 of paired-end files. This option requires paired-end
#' reads.
#' @param illumina a logical specifying that the adapter sequence to be trimmed
#' is the first 13bp of the Illumina universal adapter AGATCGGAAGAGC instead of
#' the default auto-detection of adapter sequence.  Default: \code{FALSE}
#' @param nextera adapter sequence to be trimmed is the first 12bp of the
#' Nextera adapter CTGTCTCTTATA instead of the default auto-detection of adapter
#' sequence.
#' @param small_rna a logical specifying that the adapter sequence to be trimmed
#' is the first 12bp of the Illumina Small RNA 3' Adapter TGGAATTCTCGG instead
#' of the default auto-detection of adapter sequence.  Selecting to trim
#' smallRNA adapters will also lower the \code{length} value to 18bp. If the
#' smallRNA libraries are paired-end then \code{adapter2} will be set to the
#' Illumina small RNA 5' adapter automatically (GATCGTCGGACT) unless
#' \code{adapter2} had been defined explicitly.
#' @param minlength an integer value; reads that become shorter than this length
#' as a result of either quality or adapter trimming are discarded. A value of 0
#' effectively disables this behaviour.  Default: 20 bp.  For paired-end files,
#' both reads of a read-pair need to be longer than bp to be printed out to
#' validated paired-end files. If only one read became too short there is the
#' possibility of keeping such unpaired single-end reads (see
#' \code{retain_unpaired}). Default pair-cutoff: 20 bp.
#' @param minqual an integer value specifying the quality threshold below which
#' to trim low-quality ends from reads in addition to adapter removal. Default
#' Phred score: 20.
#' @param trimN a logical specifying whether to remove Ns from the end of reads.
#' @param retainUnpaired a logical.  If only one of the two paired-end reads
#' become too short, the longer read will be written to either .unpaired_1.fq or
#' .unpaired_2.fq output files. The length cutoff for unpaired single-end reads
#' is governed by the parameters \code{retain1length} and \code{retain2length}.
#' Default: ON.
#' @param retain1length an integer.  Unpaired single-end read length cutoff
#' needed for read 1 to be written to .unpaired_1.fq output file. These reads
#' may then be mapped in single-end mode. Default: 35 bp.
#' @param retain2length an integer.  Unpaired single-end read length cutoff
#' needed for read 2 to be written to .unpaired_1.fq output file. These reads
#' may then be mapped in single-end mode. Default: 35 bp
#' @param clipR1 an integer instructing Trim Galore to remove the specified
#' number of bp from the 5' end of read 1 (or single-end reads). This may be
#' useful if the qualities were very poor, or if there is some sort of unwanted
#' bias at the 5' end. Default: 0
#' @param clipR2 an integer instructing Trim Galore to remove the specified
#' number of bp from the 5' end of read 2 (paired-end reads only). This may be
#' useful if the qualities were very poor, or if there is some sort of unwanted
#' bias at the 5' end. Default: 0
#' @param clip3primeR1 an integer instructing Trim Galore to remove the
#' specified number of bp from the 3' end of read 1 (or single-end reads) AFTER
#' adapter/quality trimming has been performed. This may remove some unwanted
#' bias from the 3' end that is not directly related to adapter sequence or
#' basecall quality. Default: 0.
#' @param clip3primeR2 an integer instructing Trim Galore to remove the
#' specified number of bp from the 3' end of read 1 (or single-end reads) AFTER
#' adapter/quality trimming has been performed. This may remove some unwanted
#' bias from the 3' end that is not directly related to adapter sequence or
#' basecall quality. Default: 0.
#' @param robust_check a logical indicating whether to check that the paired
#' files specified are matching and have equal numbers of reads.  Default:
#' \code{FALSE}
#' @param trimgalore a character string specifying the path to the trimgalore executable.
#' On Unix systems, if the executable is in \code{$PATH}, then it may be left as
#' the default. If it is not in \code{$PATH}, then the absolute path should be given.
#` If using the WSL on Windows 10, then the path must be the absolute path in WSL,
#' unless the system has been set up as described in the vignette.
#' @param dest.dir a character string specifying the output directory.  If NULL
#' a directory named "TRIMMED_FASTQC" is created in the current working directory
#' [DEFAULT = NULL].
#' @param threads an integer value indicating the number of parallel threads to
#' be used by FastQC. [DEFAULT = maximum number of available threads - 1].
#'
#' @export
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel

source('R/pgx-grep2.R',local=TRUE)


source('R/pgx-links.R',local=TRUE)


source('R/pgx-modules.R',local=TRUE)


source('R/pgx-plotting.R',local=TRUE)


source('R/pgx-predict.R',local=TRUE)


source('R/pgx-proteomics.R',local=TRUE)


source('R/pgx-pubmed.R',local=TRUE)


source('R/pgx-purifyExpression.R',local=TRUE)


source('R/pgx-signature.R',local=TRUE)


source('R/pgx-singlecell.R',local=TRUE)


source('R/pgx-tcga.R',local=TRUE)


source('R/pgx-ui.R',local=TRUE)


source('R/pgx-vizpanels.R',local=TRUE)


source('R/pgx-wordcloud.R',local=TRUE)


source('R/setdirs.R',local=TRUE)


source('R/ui-code.R',local=TRUE)


source('R/xcr-graph.r',local=TRUE)


source('R/xcr-math.r',local=TRUE)


source('shiny/boards/server/biomarker_server.R',local=TRUE)


source('shiny/boards/server/clustering_server.R',local=TRUE)


source('shiny/boards/server/compare_server.R',local=TRUE)


source('shiny/boards/server/connectivity_server.R',local=TRUE)


source('shiny/boards/server/correlation_server.R',local=TRUE)

#' DataView module server function
#'
#' @description A shiny Module (server code).
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @param inputData Reactive expression that provides the input ngs/pgx data object 
#'
#' @export 
DataViewBoard <- function(){}
source('shiny/boards/server/dataview_server.R',local=TRUE)


source('shiny/boards/server/drugconnectivity_server.R',local=TRUE)


source('shiny/boards/server/enrichment_server.R',local=TRUE)


source('shiny/boards/server/expression_server.R',local=TRUE)


source('shiny/boards/server/featuremap_server.R',local=TRUE)


source('shiny/boards/server/functional_server.R',local=TRUE)


source('shiny/boards/server/intersection_server.R',local=TRUE)


source('shiny/boards/server/loading_server.R',local=TRUE)


source('shiny/boards/server/signature_server.R',local=TRUE)


source('shiny/boards/server/singlecell_server.R',local=TRUE)


source('shiny/boards/server/tcga_server.R',local=TRUE)


source('shiny/boards/server/user_server.R',local=TRUE)


source('shiny/boards/server/wgcna_server.R',local=TRUE)


source('shiny/boards/server/wordcloud_server.R',local=TRUE)


source('shiny/boards/ui/biomarker_ui.R',local=TRUE)


source('shiny/boards/ui/clustering_ui.R',local=TRUE)


source('shiny/boards/ui/compare_ui.R',local=TRUE)


source('shiny/boards/ui/connectivity_ui.R',local=TRUE)


source('shiny/boards/ui/correlation_ui.R',local=TRUE)

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
source('shiny/boards/ui/dataview_ui.R',local=TRUE)


source('shiny/boards/ui/drugconnectivity_ui.R',local=TRUE)


source('shiny/boards/ui/enrichment_ui.R',local=TRUE)


source('shiny/boards/ui/expression_ui.R',local=TRUE)


source('shiny/boards/ui/featuremap_ui.R',local=TRUE)


source('shiny/boards/ui/functional_ui.R',local=TRUE)


source('shiny/boards/ui/intersection_ui.R',local=TRUE)


source('shiny/boards/ui/loading_ui.R',local=TRUE)


source('shiny/boards/ui/signature_ui.R',local=TRUE)


source('shiny/boards/ui/singlecell_ui.R',local=TRUE)


source('shiny/boards/ui/tcga_ui.R',local=TRUE)


source('shiny/boards/ui/user_ui.R',local=TRUE)


source('shiny/boards/ui/wgcna_ui.R',local=TRUE)


source('shiny/boards/ui/wordcloud_ui.R',local=TRUE)


source('shiny/modules/AuthenticationModule.R',local=TRUE)


source('shiny/modules/BatchCorrectModule.R',local=TRUE)


source('shiny/modules/ComputePgxModule.R',local=TRUE)


source('shiny/modules/MakeContrastModule.R',local=TRUE)


source('shiny/modules/NormalizeCountsModule.R',local=TRUE)


source('shiny/modules/PlotModule.R',local=TRUE)


source('shiny/modules/plotModules/dataviewTSNEPlotModule.R',local=TRUE)


source('shiny/modules/plotModules/examplePlotModule.R',local=TRUE)


source('shiny/modules/plotModules/PlotModule.R',local=TRUE)


source('shiny/modules/QuestionModule.R',local=TRUE)


source('shiny/modules/TimerModule.R',local=TRUE)


source('shiny/modules/UploadModule.R',local=TRUE)

#' The main application Server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
app_server <- function(){}
source('shiny/server.R',local=TRUE)

#' The main application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @export
app_ui <- function(){}
source('shiny/ui.R',local=TRUE)


source('shiny/utils/utils.R',local=TRUE)

