##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##



## Annotate clusters ##########


clustannot_plot_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)

  a_MSigDB <- "<a href='http://software.broadinstitute.org/gsea/msigdb'> GO</a>"
  a_KEGG <- "<a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102409/'> KEGG</a>"
  a_GO <- "<a href='http://geneontology.org/'>Gene Ontology</a>"

  info_text <- paste0("The top features of the heatmap in the <code>Heatmap</code> panel are divided into gene (or gene set) clusters based on their expression profile patterns. For each cluster, the platform provides a functional annotation in the <code>Annotate cluster</code> panel by correlating annotation features from more than 42 published reference databases, including well-known databases such as ", a_MSigDB, ", ", a_KEGG, " and ", a_GO, ". In the plot settings, users can specify the level and reference set to be used under the <code>Reference level</code> and <code>Reference set</code> settings, respectively.")

  plots_opts <- shiny::tagList(
    withTooltip(
      shiny::selectInput(ns("xann_level"), "Reference level:",
        choices = c("gene", "geneset", "phenotype"),
        selected = "geneset", width = "80%"
      ),
      "Select the level of an anotation analysis.",
      placement = "left", options = list(container = "body")
    ),
    shiny::conditionalPanel(
      "input.xann_level == 'geneset'",
      ns = ns,
      withTooltip(shiny::checkboxInput(ns("xann_odds_weighting"), "Fisher test weighting"),
        "Enable weighting with Fisher test probability for gene sets. This will effectively penalize small clusters and increase robustness.",
        placement = "left", options = list(container = "body")
      )
    ),
    withTooltip(shiny::selectInput(ns("xann_refset"), "Reference set:", choices = "", width = "80%"),
      "Specify a reference set to be used in the annotation.",
      placement = "left", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = "Functional annotation of clusters",
    label = label,
    outputFunc = plotly::plotlyOutput,
    outputFunc2 = plotly::plotlyOutput,
    info.text = info_text,
    options = plots_opts,
    download.fmt = c("png", "pdf", "csv"),
    width = c("auto", "100%"),
    height = height
  )
}

clustannot_table_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)
  plotWidget(ns("clustannot_table"))
}

clustannot_server <- function(id,
                              pgx,
                              top_matrix = reactive(NULL),
                              hm_level = reactive("gene"),
                              hm_topmode = reactive("sd"),
                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    shiny::observe({
      shiny::req(pgx$X, pgx$gsetX, pgx$families)

      if (is.null(input$xann_level)) {
        return(NULL)
      }
      ann.types <- sel <- NULL
      if (input$xann_level != "phenotype") {
        if (input$xann_level == "geneset") {
          ann.types <- names(COLLECTIONS)
          cc <- sapply(COLLECTIONS, function(s) length(intersect(s, rownames(pgx$gsetX))))
          ann.types <- ann.types[cc >= 3]
        }
        if (input$xann_level == "gene") {
          ann.types <- names(pgx$families)
          cc <- sapply(pgx$families, function(g) length(intersect(g, rownames(pgx$X))))
          ann.types <- ann.types[cc >= 3]
        }
        ann.types <- setdiff(ann.types, "<all>") ## avoid slow...
        ann.types <- grep("^<", ann.types, invert = TRUE, value = TRUE) ## remove special groups
        sel <- ann.types[1]
        if ("H" %in% ann.types) sel <- "H"
        j <- grep("^transcription", ann.types, ignore.case = TRUE)
        if (input$xann_level == "geneset") j <- grep("hallmark", ann.types, ignore.case = TRUE)
        if (length(j) > 0) sel <- ann.types[j[1]]
        ann.types <- sort(ann.types)
      } else {
        ann.types <- sel <- "<all>"
      }
      shiny::updateSelectInput(session, "xann_refset", choices = ann.types, selected = sel)
    })

    ## This is used both for plot and table
    get_annot_correlation <- shiny::reactive({
      shiny::req(pgx$X, pgx$Y, pgx$gsetX, pgx$families)

      ## filt <- getTopMatrix()
      filt <- top_matrix()
      shiny::req(filt)

      zx <- filt$mat
      idx <- filt$idx
      samples <- filt$samples

      if (nrow(zx) <= 1) {
        return(NULL)
      }

      ann.level <- "geneset"
      ann.refset <- "Hallmark collection"
      ann.level <- input$xann_level
      ## if(is.null(ann.level)) return(NULL)
      ann.refset <- input$xann_refset
      ## if(is.null(ann.refset)) return(NULL)
      shiny::req(input$xann_level, input$xann_refset)

      ref <- NULL
      ref <- pgx$gsetX[, , drop = FALSE]
      ref <- pgx$X[, , drop = FALSE]
      if (ann.level == "gene" && ann.refset %in% names(pgx$families)) {
        gg <- pgx$families[[ann.refset]]
        jj <- match(toupper(gg), toupper(pgx$genes$gene_name))
        jj <- setdiff(jj, NA)
        pp <- rownames(pgx$genes)[jj]
        ref <- pgx$X[intersect(pp, rownames(pgx$X)), , drop = FALSE]
      }
      if (ann.level == "geneset" && ann.refset %in% names(COLLECTIONS)) {
        ss <- COLLECTIONS[[ann.refset]]
        ss <- intersect(ss, rownames(pgx$gsetX))
        length(ss)
        ref <- pgx$gsetX[ss, ]
      }
      if (ann.level == "phenotype") {
        ref <- t(playbase::expandAnnotationMatrix(pgx$Y))
      }
      if (is.null(ref)) {
        cat("<clustering:get_annot_correlation> WARNING:: ref error\n")
        return(NULL)
      }

      ## -----------  restrict to top??
      dim(ref)
      if (nrow(ref) > 1000) {
        ref <- head(ref[order(-apply(ref, 1, sd)), ], 1000)
      }

      ## -----------  get original data level
      X <- pgx$X
      if (hm_level() == "geneset") X <- pgx$gsetX

      ## ----------- for each gene cluster compute average correlation
      idxx <- setdiff(idx, c(NA, " ", "   "))
      rho <- matrix(NA, nrow(ref), length(idxx))
      colnames(rho) <- idxx
      rownames(rho) <- rownames(ref)

      i <- 1
      if (nrow(ref) > 0) {
        for (i in 1:length(idxx)) {
          gg <- rownames(zx)[which(idx == idxx[i])]
          aa <- t(X[gg, samples, drop = FALSE])
          bb <- t(ref[, samples, drop = FALSE])
          ## rr = cor(aa , bb, use="pairwise", method="spearman")
          rr <- cor(apply(aa, 2, rank), apply(bb, 2, rank), use = "pairwise")
          if (hm_topmode() == "pca") rr <- abs(rr)
          rho[, i] <- colMeans(rr, na.rm = TRUE)
        }
      }

      if (input$hm_level == "gene" && ann.level == "geneset" && input$xann_odds_weighting) {
        table(idx)
        grp <- tapply(toupper(rownames(zx)), idx, list) ## toupper for mouse!!
        ## gmt <- GSETS[rownames(rho)]
        gmt <- getGSETS(rownames(rho))
        bg.genes <- toupper(rownames(X))
        P <- c()
        for (i in 1:ncol(rho)) {
          k <- colnames(rho)[i]
          res <- playbase::gset.fisher(
            grp[[k]], gmt,
            fdr = 1, min.genes = 0, max.genes = Inf,
            background = bg.genes
          )
          res <- res[rownames(rho), ]
          r <- res[, "odd.ratio"]
          odd.prob <- r / (1 + r)
          ## odd.1mpv <- 1 - res[,"p.value"]
          ## P <- cbind(P,odd.1mpv)
          P <- cbind(P, odd.prob)
        }
        colnames(P) <- colnames(rho)
        rownames(P) <- rownames(rho)
        rho <- rho * (P / max(P))
      }

      ## rho = round(rho, digits=3)
      dim(rho)
      return(rho)
    })


    ## Plot ##########

    plot_data <- shiny::reactive({
      get_annot_correlation()
    })

    plot.RENDER <- function() {
      rho <- plot_data()
      shiny::req(rho)

      NTERMS <- 6
      NTERMS <- 12
      slen <- 40
      if (ncol(rho) >= 5) {
        slen <- 20
      }
      if (ncol(rho) > 6) {
        NTERMS <- 6
      }
      if (ncol(rho) <= 2) {
        NTERMS <- 22
      }

      klrpal <- rep(RColorBrewer::brewer.pal(8, "Set2"), 2)
      ## klrpal = paste0(klrpal,"88")
      col.addalpha <- function(clr, a = 100) {
        paste0("rgba(", paste(col2rgb(clr)[, 1], collapse = ","), ",", a, ")")
      }
      ## klrpal = as.character(sapply(klrpal, col.addalpha, a=50))
      klrpal <- paste0(klrpal, "55")

      plot_list <- list()
      i <- 1
      for (i in 1:min(9, ncol(rho))) {
        x <- rev(head(sort(rho[, i], decreasing = TRUE), NTERMS))
        names(x) <- sub(".*:", "", names(x))
        names(x) <- gsub(GSET.PREFIX.REGEX, "", names(x))

        y <- names(x)
        y <- factor(y, levels = y)
        anntitle <- function(tt) {
          list(
            text = tt, font = list(size = 13),
            xref = "paper", yref = "paper",
            yanchor = "bottom", xanchor = "center",
            align = "center", x = 0.5, y = 1.02, showarrow = FALSE
          )
        }

        ## NOTE: The same plotly code (originally) as in `clustering_server.R`
        ##       -> Seems it uses the function from that file, not this one
        ## TODO: clean-up; we should stick to the general setup of individual
        ##       scripts for the plotting functions
        plot_list[[i]] <- plotly::plot_ly(
          x = x, y = y, type = "bar", orientation = "h",
          ## text=y,
          hoverinfo = "text",
          hovertemplate = paste0("%{y}<extra>", colnames(rho)[i], "</extra>"),
          ## hovertemplate = "%{y}",
          marker = list(color = klrpal[i])
        ) %>%
          plotly::layout(
            showlegend = FALSE,
            annotations = anntitle(colnames(rho)[i]),
            ## annotations = list(text="TITLE"),
            ## margin = c(0, 0.0, 0.05, 0.05),
            margin = list(l = 5, r = 0, t = 25, b = 15),
            xaxis = list(
              range = c(0, 0.9),
              titlefont = list(size = 11),
              tickfont = list(size = 10),
              showgrid = FALSE,
              title = "\ncorrelation (R)"
            ),
            yaxis = list(
              title = "",
              showgrid = FALSE,
              showline = FALSE,
              showticklabels = FALSE,
              showgrid = FALSE,
              zeroline = FALSE
            )
          ) %>%
          ## labeling the y-axis inside bars
          plotly::add_annotations(
            xref = "paper", yref = "y",
            x = 0.01, y = y, xanchor = "left",
            text = playbase::shortstring(y, slen),
            font = list(size = 10),
            showarrow = FALSE, align = "right"
          )
        ## layout(margin = c(0,0,0,0))
      }

      if (length(plot_list) <= 4) {
        nrows <- ceiling(length(plot_list) / 2)
      } else {
        nrows <- ceiling(length(plot_list) / 3)
      }

      plotly::subplot(plot_list,
        nrows = nrows, shareX = TRUE,
        ## template = "plotly_dark",
        margin = c(0, 0.0, 0.05, 0.05)
      ) %>%
        plotly::config(displayModeBar = FALSE)
    }

    modal_plot.RENDER <- function() {
      plot.RENDER()
    }

    PlotModuleServer(
      "pltmod",
      ## plotlib = "plotly",
      ## plotlib2 = "plotly",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      renderFunc = plotly::renderPlotly,
      renderFunc2 = plotly::renderPlotly,
      res = c(90, 170), ## resolution of plots
      pdf.width = 8, pdf.height = 5,
      add.watermark = watermark
    )


    ## Table ##########

    table.RENDER <- shiny::reactive({
      rho <- get_annot_correlation()
      if (is.null(rho)) {
        return(NULL)
      }

      ## rownames(rho) = playbase::shortstring(rownames(rho),50)
      rho.name <- playbase::shortstring(sub(".*:", "", rownames(rho)), 60)
      ## rho = data.frame(cbind( name=rho.name, rho))
      df <- data.frame(feature = rho.name, round(as.matrix(rho), digits = 3))
      rownames(df) <- rownames(rho)
      if (input$xann_level == "geneset") {
        df$feature <- playbase::wrapHyperLink(df$feature, rownames(df))
      }

      DT::datatable(
        df,
        rownames = FALSE, escape = c(-1, -2),
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = c(1)),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip", buttons = c("copy", "csv", "pdf"),
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE, ## scrollY = TRUE,
          ## scrollY = 170,
          scrollY = "70vh",
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    table_info_text <- "In this table, users can check mean correlation values of features in the clusters with respect to the annotation references database selected in the settings."

    ## clustannot_table_module <- tableModule(
    clustannot_table_module <- shiny::callModule(
      tableModule,
      id = "clustannot_table",
      func = table.RENDER,
      ## options = clustannot_table_opts,
      info.text = table_info_text,
      title = "Annotation scores", label = "b",
      height = c(240, 700), width = c("auto", 1000),
      ## caption = clustannot_caption
    )
  })
}
