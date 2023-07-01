##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

upload_plot_phenostats_ui <- function(
        id,
        label = "",
        height,
        width,
        title,
        caption,
        info.text) {
    ns <- shiny::NS(id)

    PlotModuleUI(
        ns("pltmod"),
        title = title,
        label = label,
        plotlib = "base",
        info.text = info.text,
        caption = caption,
        options = NULL,
        download.fmt = c("png", "pdf", "csv"),
        width = width,
        height = height
    )
}

upload_plot_phenostats_server <- function(id, checkTables, uploaded, watermark = FALSE) {
    moduleServer(id, function(input, output, session) {

        ## extract data from pgx object
        plot_data <- shiny::reactive({

            check <- checkTables()
            status.ok <- check["samples.csv", "status"]
            if (status.ok != "OK") {
                frame()
                status.ds <- check["samples.csv", "description"]
                msg <- paste(
                    toupper(status.ok), "\n", "(Required) Upload 'samples.csv'",
                    tolower(status.ds)
                )
                graphics::text(0.5, 0.5, paste(strwrap(msg, 30), collapse = "\n"), col = "grey25")
                graphics::box(lty = 1, col = "grey60")
                return(NULL)
            } else {
                pheno <- uploaded[["samples.csv"]]
                return(pheno)
            }

        })

        plot.RENDER <- function() {
            pheno <- plot_data()

            shiny::validate(
                shiny::need(!is.null(pheno),
                    "Please upload a samples.csv file (REQUIRED)."
                )
            )

            px <- head(colnames(pheno), 20) ## show maximum??

            df <- type.convert(pheno[, px, drop = FALSE])
            vt <- df %>% inspectdf::inspect_types()
            vt

            ## discretized continuous variable into 10 bins
            ii <- unlist(vt$col_name[c("numeric", "integer")])
            ii
            if (!is.null(ii) && length(ii)) {
                cat("[UploadModule::phenoStats] discretizing variables:", ii, "\n")
                df[, ii] <- apply(df[, ii, drop = FALSE], 2, function(x) {
                    if (any(is.infinite(x))) x[which(is.infinite(x))] <- NA
                    cut(x, breaks = 10)
                })
            }

            p1 <- df %>%
                inspectdf::inspect_cat() %>%
                inspectdf::show_plot()
            tt2 <- paste(nrow(pheno), "samples x", ncol(pheno), "phenotypes")

            p1 <- p1 + ggplot2::ggtitle("PHENOTYPES", subtitle = tt2) +
                ggplot2::theme(

                    axis.text.y = ggplot2::element_text(
                        size = 12,
                        margin = ggplot2::margin(0, 0, 0, 25),
                        hjust = 1
                    )
                )
            return(p1)
        }

        modal_plot.RENDER <- function() {
            plot.RENDER()
        }

        PlotModuleServer(
            "pltmod",
            plotlib = "base",
            func = plot.RENDER,
            func2 = modal_plot.RENDER,
            csvFunc = plot_data, ##  *** downloadable data as CSV
            res = c(90, 90), ## resolution of plots
            pdf.width = 4, pdf.height = 4,
            add.watermark = watermark
        )
    }) ## end of moduleServer
}
