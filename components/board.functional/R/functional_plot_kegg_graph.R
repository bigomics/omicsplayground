##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' Importance plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
functional_plot_kegg_graph_ui <- function(id,
                                          label = "",
                                          rowH = 660) {
  ns <- shiny::NS(id)
  info_text <- strwrap("<strong>KEGG pathways</strong> are a collection of
    manually curated pathways representing the current knowledge of molecular
    interactions, reactions and relation networks as pathway maps. In the
    pathway map, genes are colored according to their upregulation (red) or
    downregulation (blue) in the contrast profile. Each pathway is scored for
    the selected contrast profile and reported in the table below.")

  PlotModuleUI(ns("plot"),
    title = "Kegg pathway map",
    label = label,
    plotlib = "image",
    info.text = info_text,
    info.width = "350px",
    options = NULL,
    download.fmt = "png",
    height = c(0.53 * rowH, 700),
    width = c("100%", 1280),
  )
}

#' Importance plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
functional_plot_kegg_graph_server <- function(id,
                                              pgx,
                                              getFilteredKeggTable,
                                              kegg_table,
                                              fa_contrast,
                                              watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
        res <- list(
          pgx = pgx,
          df = getFilteredKeggTable(),
          kegg_table = kegg_table,
          fa_contrast = fa_contrast()
        )
        return(res)
      })

      plot_RENDER <- shiny::reactive({
        res <- plot_data()
        pgx <- res$pgx
        df <- res$df
        fa_contrast <- res$fa_contrast
        kegg_table <- res$kegg_table

        ###############

        NULL.IMG <- list(src = "", contentType = "image/png")
        if (is.null(pgx)) {
          return(NULL.IMG)
        }

        comparison <- fa_contrast
        if (is.null(comparison) || length(comparison) == 0) {
          return(NULL.IMG)
        }
        if (comparison == "") {
          return(NULL.IMG)
        }

        ## get fold-change vector
        fc <- pgx$gx.meta$meta[[comparison]]$meta.fx
        pp <- rownames(pgx$gx.meta$meta[[comparison]])

        if ("hgnc_symbol" %in% colnames(pgx$genes)) {
          names(fc) <- pgx$genes[pp, "hgnc_symbol"]
        } else {
          names(fc) <- toupper(pgx$genes[pp, "gene_name"])
        }
        fc <- fc[order(-abs(fc))]
        fc <- fc[which(!duplicated(names(fc)) & names(fc) != "")]

        ## get selected KEGG id
        # df <- getFilteredKeggTable()
        if (is.null(df)) {
          return(NULL.IMG)
        }

        sel.row <- kegg_table$rows_selected()
        if (is.null(sel.row) || length(sel.row) == 0) {
          return(NULL.IMG)
        }
        sel.row <- as.integer(sel.row)

        pathway.id <- "04110" ## CELL CYCLE
        pathway.name <- pw.genes <- "x"
        if (is.null(sel.row) || length(sel.row) == 0) {
          return(NULL.IMG)
        }

        if (!is.null(sel.row) && length(sel.row) > 0) {
          pathway.id <- df[sel.row, "kegg.id"]
          pathway.name <- df[sel.row, "pathway"]
          pw.genes <- unlist(getGSETS(as.character(pathway.name)))
        }

        ## folder with predownloaded XML files
        xml.dir <- file.path(FILES, "kegg-xml")
        xml.dir <- normalizePath(xml.dir) ## absolute path

        ## We temporarily switch the working directory to always readable
        ## TMP folder
        curwd <- getwd()
        tmpdir <- tempdir()
        setwd(tmpdir)
        pv.out <- pathview::pathview(
          gene.data = fc, pathway.id = pathway.id, gene.idtype = "SYMBOL",
          gene.annotpkg = "org.Hs.eg.db", species = "hsa",
          out.suffix = "pathview", limit = list(gene = 2, cpd = 1),
          low = list(gene = "dodgerblue2", cpd = "purple"),
          high = list(gene = "firebrick2", cpd = "yellow"),
          kegg.dir = xml.dir, kegg.native = TRUE, same.layer = FALSE
        )
        Sys.sleep(0.2) ## wait for graph

        ## back to previous working folder
        setwd(curwd)

        outfile <- file.path(tmpdir, paste0("hsa", pathway.id, ".pathview.png"))
        if (!file.exists(outfile)) {
          return(NULL.IMG)
        }

        list(
          src = outfile,
          contentType = "image/png",
          width = "100%", height = "100%", ## actual size: 1040x800
          alt = "pathview image"
        )
      })

      PlotModuleServer(
        "plot",
        plotlib = "image",
        func = plot_RENDER,
        func2 = plot_RENDER,
        csvFunc = plot_data,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
