##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

upload_plot_countstats_ui <- function(
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

upload_plot_countstats_server <- function(id, checkTables, uploaded, watermark = FALSE) {
    moduleServer(id, function(input, output, session) {

        ## extract data from pgx object
        plot_data <- shiny::reactive({
            check <- checkTables()
            req(check)
            status.ok <- check["counts.csv", "status"]
            dbg("[countStats] status.ok = ", status.ok)

            if (status.ok != "OK") {
                frame()
                status.ds <- check["counts.csv", "description"]
                msg <- paste(
                    toupper(status.ok), "\n", "(Required) Upload 'counts.csv'",
                    tolower(status.ds)
                )
                graphics::text(0.5, 0.5, paste(strwrap(msg, 30), collapse = "\n"), col = "grey25")
                graphics::box(lty = 1, col = "grey60")
                return(NULL)
            } else {
                counts <- uploaded[["counts.csv"]]
                return(counts)
            }
        })

        plot.RENDER <- function() {
            counts <- plot_data()

            shiny::validate(
                shiny::need(
                    !is.null(counts),
                    "Please upload a counts.csv file (REQUIRED)."
                )
            )

            xx <- log2(1 + counts)
            if (nrow(xx) > 1000) xx <- xx[sample(1:nrow(xx), 1000), , drop = FALSE]
            suppressWarnings(dc <- data.table::melt(xx))
            dc$value[dc$value == 0] <- NA
            tt2 <- paste(nrow(counts), "genes x", ncol(counts), "samples")
            ggplot2::ggplot(dc, ggplot2::aes(x = value, color = Var2)) +
                ggplot2::geom_density() +
                ggplot2::xlab("log2(1+counts)") +
                ggplot2::theme(legend.position = "none") +
                ggplot2::ggtitle("COUNTS", subtitle = tt2)
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
