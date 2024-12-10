##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


PlotModuleUI <- function(id,
                         info.text = "Figure",
                         info.methods = NULL,
                         info.references = NULL,
                         info.extra_link = NULL,
                         title = "",
                         options = NULL,
                         label = "",
                         caption = "",
                         caption2 = caption,
                         plotlib = "base",
                         plotlib2 = NULL,
                         outputFunc = NULL,
                         outputFunc2 = outputFunc,
                         no.download = FALSE,
                         download.fmt = c("png", "pdf"),
                         just.info = FALSE,
                         info.width = "300px",
                         show.maximize = TRUE,
                         height = c(400, 800),
                         card_footer_height = "3rem",
                         width = c("auto", "100%"),
                         pdf.width = 8,
                         pdf.height = 8,
                         cards = FALSE,
                         card_names = NULL,
                         header_buttons = NULL,
                         translate = TRUE,
                         translate_js = TRUE) {
  ns <- shiny::NS(id)

  if (is.null(plotlib2)) plotlib2 <- plotlib
  if (length(height) == 1) height <- c(height, 800)
  if (length(width) == 1) width <- c(width, "100%")

  ifnotchar.int <- function(s) {
    suppressWarnings(
      ifelse(!is.na(as.integer(s)), paste0(as.integer(s), "px"), s)
    )
  }
  width.1 <- ifnotchar.int(width[1])
  width.2 <- "100%"
  height.1 <- ifnotchar.int(height[1])
  height.2 <- ifnotchar.int(height[2])

  if (translate) {
    info.text <- tspan(info.text, js = translate_js)
    info.methods <- tspan(info.methods, js = translate_js)
    title <- tspan(title, js = translate_js)
    caption2 <- tspan(caption2, js = translate_js)
    caption <- tspan(caption, js = translate_js)
  }

  getOutputFunc <- function(plotlib) {
    FUN <- switch(plotlib,
      generic = NULL,
      htmlwidget = NULL,
      plotly = plotly::plotlyOutput,
      pairsD3 = pairsD3::pairsD3Output,
      visnetwork = visNetwork::visNetworkOutput,
      ggplot = shiny::plotOutput,
      grid = shiny::plotOutput,
      iheatmapr = iheatmapr::iheatmaprOutput,
      image = shiny::imageOutput,
      base = shiny::plotOutput,
      svgPanZoom = svgPanZoom::svgPanZoomOutput,
      renderUI = shiny::htmlOutput,
      shiny::plotOutput
    )
    FUN
  }

  if (is.null(plotlib2)) plotlib2 <- plotlib
  if (cards) {
    if (length(plotlib) != length(card_names)) {
      plotlib <- rep(plotlib[1], length(card_names))
    }
    if (length(outputFunc) == 1) {
      outputFunc <- rep(list(outputFunc), length(card_names))
    }
    if (length(outputFunc2) == 1) {
      outputFunc2 <- rep(list(outputFunc2), length(card_names))
    }
    if (is.null(outputFunc)) outputFunc <- lapply(plotlib, getOutputFunc)
    if (is.null(outputFunc2)) outputFunc2 <- lapply(plotlib2, getOutputFunc)
  } else {
    if (is.null(outputFunc)) outputFunc <- getOutputFunc(plotlib)
    if (is.null(outputFunc2)) outputFunc2 <- getOutputFunc(plotlib2)
  }

  ## --------------------------------------------------------------------------------
  ## ------------------------ BUTTONS -----------------------------------------------
  ## --------------------------------------------------------------------------------
  options.button <- ""

  if (!just.info && !is.null(options) && length(options) > 0) {
    options.button <- DropdownMenu(
      options,
      size = "xs",
      width = "auto",
      icon = shiny::icon("bars"),
      status = "default"
    )
  }

  if (cards == FALSE) {
    download_buttons <- div(
      shiny::downloadButton(
        outputId = ns("download"),
        label = "Download",
        class = "btn-outline-primary"
      )
    )
  } else {
    button_list <- lapply(seq_along(card_names), function(x) {
      shiny::conditionalPanel(
        condition = paste0(
          "input.card_selector == '", card_names[x], "'"
        ),
        ns = ns,
        div(
          shiny::downloadButton(
            outputId = ns(paste0(
              "download", x
            )),
            label = "Download",
            class = "btn-outline-primary"
          )
        )
      )
    })
    download_buttons <- do.call(div, button_list)
  }

  pdf_size_ui <- shiny::tagList(
    shiny::fillRow(
      shiny::numericInput(ns("pdf_width"),
        "Width",
        pdf.width,
        1, 20, 1,
        width = "95%"
      ),
      shiny::numericInput(ns("pdf_height"),
        "Height",
        pdf.height,
        1, 20, 1,
        width = "100%"
      )
    ),
    shiny::br(), shiny::br(), shiny::br()
  )

  if ("csv" %in% download.fmt) download.fmt <- c(download.fmt, "excel")

  dload.button <- DropdownMenu(
    div(
      style = "width: 150px;",
      shiny::selectInput(
        inputId = ns("downloadOption"),
        label = "Format",
        choices = download.fmt
      ),
      div(
        id = ns("pdf_size_panel"),
        shiny::div(
          pdf_size_ui,
          shiny::br()
        )
      ),
      shiny::conditionalPanel(
        condition = "input.downloadOption == 'pdf'",
        ns = ns,
        shiny::checkboxInput(
          inputId = ns("get_pdf_settings"),
          label = "Include plot settings",
          TRUE
        )
      ),
      download_buttons,
    ),
    size = "xs",
    icon = shiny::icon("download"),
    status = "default"
  )

  if (no.download || length(download.fmt) == 0) dload.button <- ""

  zoom.button <- NULL
  if (show.maximize) {
    zoom.button <- modalTrigger(
      ns("zoombutton"),
      ns("plotPopup"),
      icon("up-right-and-down-left-from-center"),
      class = "btn-circle-xs"
    )
  }

  # Build cards or single plot
  if (cards) {
    tabs <- lapply(1:length(card_names), function(x) {
      bslib::nav_panel(
        card_names[x],
        outputFunc[[x]](ns(paste0("renderfigure", x))) %>%
          bigLoaders::useSpinner()
      )
    })
    tabs <- c(tabs, title = "", id = ns("card_selector"))
    plot_cards <- do.call(
      bslib::navset_card_pill,
      tabs
    )
  } else {
    plot_cards <- outputFunc(ns("renderfigure")) %>%
      bigLoaders::useSpinner()
  }

  if (is.null(header_buttons)) {
    header_buttons <- div()
  }

  header <- shiny::fillRow(
    flex = c(1, NA, NA, NA, NA, NA, NA),
    class = "plotmodule-header",
    shiny::div(
      class = "plotmodule-title",
      style = "white-space: nowrap; overflow: hidden; text-overflow: clip;",
      title
    ),
    if (cards) {
      nav_bar <- gsub("nav nav-pills shiny-tab-input card-header-pills", "nav navbar-nav shiny-tab-input header-nav", plot_cards$children[[1]])
      nav_bar <- gsub("card-header bslib-navs-card-title", "bslib-navs-card-title", nav_bar) |> shiny::HTML()
      nav_bar
    } else {
      shiny::div()
    },
    header_buttons,
    DropdownMenu(
      shiny::div(
        class = "plotmodule-info",
        shiny::HTML("<b>Plot info</b><br>"),
        shiny::HTML(as.character(info.text))
      ),
      if (!is.null(info.methods)) {
        shiny::div(
          class = "plotmodule-info",
          shiny::HTML("<b>Methods</b><br>"),
          shiny::HTML(info.methods)
        )
      } else {
        NULL
      },
      if (!is.null(info.references)) {
        html_code <- ""
        for (i in seq_along(info.references)) {
          ref <- info.references[[i]]
          name <- ref[[1]]
          link <- ref[[2]]

          # Create the formatted HTML string
          formatted_ref <- paste0("[", i, "] ", name, " <a href='", link, "' target='_blank'>", link, "</a><br>")

          # Append the formatted string to the HTML code
          html_code <- paste0(html_code, formatted_ref)
        }
        shiny::div(
          class = "plotmodule-info",
          shiny::HTML("<b>References</b>"),
          shiny::div(
            class = "plotmodule-info plotmodule-references",
            shiny::HTML(html_code)
          )
        )
      } else {
        NULL
      },
      if (!is.null(info.extra_link)) {
        shiny::div(
          class = "plotmodule-info",
          shiny::HTML(
            paste0(
              "<b><a href='",
              info.extra_link,
              "' target='_blank'>Further information...</a></b>"
            )
          )
        )
      } else {
        NULL
      },
      shiny::HTML("<br>"),
      shiny::actionButton(
        ns("copy_info"),
        "Copy text",
        icon = shiny::icon("clipboard"),
        class = "btn-outline-dark btn-sm",
        onclick = "copyPlotModuleInfo();"
      ),
      size = "xs",
      icon = shiny::icon("info"),
      status = "default",
      width = "300px"
    ),
    options.button,
    shiny::div(class = "download-button", title = "download", dload.button),
    shiny::div(class = "zoom-button", title = "zoom", zoom.button)
  )

  ## ------------------------------------------------------------------------
  ## --------------- modal UI (former output$popupfig) ----------------------
  ## ------------------------------------------------------------------------

  if (cards) {
    tabs_modal <- lapply(1:length(card_names), function(x) {
      bslib::nav_panel(
        card_names[x],
        id = card_names[x],
        bslib::card_body(
          outputFunc2[[x]](ns(paste0("renderpopup", x)),
            width = width.2, height = height.2
          ) %>%
            bigLoaders::useSpinner()
        )
      )
    })
    tabs_modal <- c(tabs_modal, id = "card_selector_modal", bg = "transparent", inverse = FALSE)
    plot_cards_modal <- do.call(
      bslib::navset_bar,
      tabs_modal
    )
    plot_cards_modal[[1]] <- gsub("nav navbar-nav nav-underline", "nav navbar-nav", plot_cards_modal[[1]]) |> shiny::HTML()
    plot_cards_modal[[1]] <- gsub("navbar navbar-default navbar-static-top", "navbar navbar-default navbar-static-top navbar-custom", plot_cards_modal[[1]]) |> shiny::HTML()
  } else {
    plot_cards_modal <- outputFunc2(ns("renderpopup"), width = width.2, height = height.2) %>%
      bigLoaders::useSpinner()
  }


  popupfigUI <- function() {
    ## render caption2 (for modal)
    if (any(class(caption2) == "reactive")) {
      caption2 <- caption2()
    }
    caption2 <- shiny::div(
      class = "caption2 popup-plot-caption",
      shiny::HTML(paste0(
        "<b>", as.character(title), ".</b>&nbsp;&nbsp;",
        as.character(caption2)
      ))
    )
    shiny::div(
      class = "popup-plot-body",
      shiny::div(
        class = "popup-plot",
        plot_cards_modal
      ),
      caption2
    )
  }

  popupfigUI_editor <- function(card = NULL) {
    if (!is.null(card)) {
      htmlOutput(ns(paste0("editor_frame", card)))
    } else {
      htmlOutput(ns("editor_frame"))
    }
  }

  ## inline styles (should be in CSS...)
  modaldialog.style <- paste0("#", ns("plotPopup"), " .modal-dialog {width:", width.2, ";}")
  modalbody.style <- paste0("#", ns("plotPopup"), " .modal-body {min-height:", height.2, "; padding:30px 150px;}")
  modalcontent.style <- paste0("#", ns("plotPopup"), " .modal-content {width:100vw;}")
  modalfooter.none <- paste0("#", ns("plotPopup"), " .modal-footer{display:none;}")

  if (any(class(caption) == "reactive")) {
    caption <- caption()
  }

  e <- bslib::card(
    class = "plotmodule",
    full_screen = FALSE,
    style = paste0("height:", height.1, ";overflow: visible;"),
    bslib::as.card_item(div(header)),
    bslib::card_body(
      gap = "0px",
      if (cards) {
        plot_cards$children[[2]]
      } else {
        plot_cards
      },
      shiny::div(
        class = "popup-modal",
        modalUI(
          id = ns("plotPopup"),
          title = title,
          size = "fullscreen",
          footer = NULL,
          popupfigUI()
        )
      ),
      if (cards) {
        div(
          lapply(1:length(card_names), function(x) {
            shiny::div(
              class = "popup-modal",
              modalUI(
                id = ns(paste0("plotPopup_editor", x)),
                title = "Editor",
                size = "fullscreen",
                footer = NULL,
                popupfigUI_editor(x)
              )
            )
          })
        )
      } else {
        shiny::div(
          class = "popup-modal",
          modalUI(
            id = ns("plotPopup_editor"),
            title = "Editor",
            size = "fullscreen",
            footer = NULL,
            popupfigUI_editor()
          )
        )
      },
      shiny::tagList(
        shiny::tags$head(shiny::tags$style(modaldialog.style)),
        shiny::tags$head(shiny::tags$style(modalbody.style)),
        shiny::tags$head(shiny::tags$style(modalcontent.style)),
        shiny::tags$head(shiny::tags$style(modalfooter.none)),
        shiny::tags$script(src = "custom/dropdown-helper.js")
      )
    ),
    bslib::card_body(
      class = "card-footer", # center the content horizontally and vertically
      style = paste0("height:", card_footer_height, ";"), # add left and top margin of 2 pixels
      div(
        class = "caption",
        shiny::HTML(paste0(
          "<b>", as.character(title), ".</b>&nbsp;",
          as.character(caption)
        ))
      )
    )
  ) # end of card
  return(e)
}


PlotModuleServer <- function(id,
                             func,
                             func2 = NULL,
                             plotlib = "base",
                             plotlib2 = plotlib,
                             renderFunc = NULL,
                             renderFunc2 = renderFunc,
                             csvFunc = NULL,
                             download.fmt = c("png", "pdf"),
                             height = c(640, 800),
                             width = c("auto", 1400),
                             res = c(100, 170),
                             download.pdf = NULL,
                             download.png = NULL,
                             download.html = NULL,
                             download.csv = NULL,
                             download.excel = NULL,
                             download.obj = NULL,
                             pdf.width = 8,
                             pdf.height = 6,
                             pdf.pointsize = 12,
                             add.watermark = FALSE,
                             remove_margins = FALSE,
                             vis.delay = 0,
                             card = NULL) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      filename <- sub("-$", "", ns("")) ## filename root

      observeEvent(input$downloadOption,
        {
          if (!input$downloadOption %in% c("pdf", "png")) {
            shinyjs::hide("pdf_size_panel")
          } else {
            shinyjs::show("pdf_size_panel")
          }
        },
        ignoreInit = TRUE
      )

      ## --------------------------------------------------------------------------------
      ## ------------------------ Plotly editor------------------------------------------
      ## --------------------------------------------------------------------------------
      octocat <- list(
        name = "Editor",
        icon = list(
          path = "M410.052,46.398c-0.812-10.885-5.509-21.129-13.226-28.845c-16.089-16.089-41.044-17.965-59.34-4.459l-7.427,5.487C257.281,72.291,191.872,135.46,135.647,206.336c-14.115,17.795-27.792,36.287-40.715,55.015c-0.928-0.042-1.859-0.068-2.795-0.068c-16.279,0-31.583,6.339-43.094,17.851C28.607,299.57,27.77,319.906,26.96,339.572c-0.745,18.1-1.449,35.196-16.99,54.271L0,406.081h15.785c37.145,0,96.119-17.431,119.447-40.759c11.511-11.511,17.85-26.815,17.85-43.094c0-0.941-0.026-1.877-0.068-2.81c18.766-12.941,37.258-26.614,55.01-40.704C278.873,222.52,342.046,157.11,395.787,84.302l5.479-7.419C407.747,68.111,410.867,57.284,410.052,46.398z M124.625,354.715c-16.334,16.334-58.786,31.89-94.095,35.555c10.098-18.012,10.791-34.866,11.417-50.082c0.754-18.326,1.406-34.152,17.702-50.449c8.678-8.678,20.216-13.457,32.488-13.457s23.81,4.779,32.488,13.457s13.457,20.215,13.457,32.487C138.082,334.5,133.303,346.037,124.625,354.715z M135.232,279.133c-6.875-6.875-15.11-11.889-24.091-14.825c10.801-15.429,22.107-30.656,33.724-45.426c12.79,1.717,24.7,7.567,33.871,16.737c9.174,9.174,15.027,21.087,16.745,33.875c-14.743,11.601-29.97,22.905-45.427,33.719C147.116,294.236,142.104,286.006,135.232,279.133z M389.2,67.971l-5.481,7.421c-50.415,68.302-109.268,129.976-175.037,183.518c-3.279-12.747-9.915-24.473-19.34-33.897c-9.421-9.421-21.145-16.055-33.893-19.333C209.017,139.887,270.692,81.036,338.97,30.649l7.427-5.488c12.275-9.062,29.023-7.801,39.823,3c5.177,5.177,8.329,12.05,8.874,19.355C395.641,54.822,393.548,62.086,389.2,67.971z",
          transform = "scale(2.2)",
          width = "1000",
          height = "1000"
        ),
        click = htmlwidgets::JS(
          if (!is.null(card)) {
            paste0("function(gd){$('#", ns(paste0("plotPopup_editor", card)), "').modal('show')}")
          } else {
            paste0("function(gd){$('#", ns("plotPopup_editor"), "').modal('show')}")
          }
        )
      )

      getEditorUrl <- function(session, path_object, path_object2) {
        cd <- session$clientData
        sprintf(
          "%s/?plotURL=%s//%s:%s%s%s&plotDS=%s//%s:%s%s%s",
          "custom/editor/index.html",
          cd$url_protocol,
          cd$url_hostname,
          cd$url_port,
          cd$url_pathname,
          path_object,
          cd$url_protocol,
          cd$url_hostname,
          cd$url_port,
          cd$url_pathname,
          path_object2
        )
      }


      if (!is.null(card)) {
        output[[paste0("editor_frame", card)]] <- renderUI({
          plot <- func()
          if (exists("csvFunc") && is.function(csvFunc)) {
            plot_data_csv <- csvFunc()
            if (inherits(plot_data_csv, "list")) {
              plot_data_csv <- plot_data_csv[[1]]
            }
          } else {
            plot_data_csv <- NULL
          }
          for (i in 1:length(plot$x$data)) {
            plot$x$data[[i]]$hovertemplate <- NULL
          }
          for (i in 1:length(plot$x$attrs)) {
            plot$x$attrs[[i]]$hovertemplate <- NULL
          }
          json <- plotly::plotly_json(plot, TRUE) # requires `listviewer` to work properly
          res <- session$registerDataObj(
            "plotly_graph", json$x$data,
            function(data, req) {
              httpResponse(
                status = 200,
                content_type = "application/json",
                content = data
              )
            }
          )
          res2 <- session$registerDataObj(
            "plotly_graph2", jsonlite::toJSON(list(data = as.list(plot_data_csv))),
            function(data, req) {
              httpResponse(
                status = 200,
                content_type = "application/json",
                content = data
              )
            }
          )
          url <- getEditorUrl(session, res, res2)
          tags$iframe(src = url, style = "height: 85vh; width: 100%;")
        })
      } else {
        output$editor_frame <- renderUI({
          plot <- func()
          if (exists("csvFunc") && is.function(csvFunc)) {
            plot_data_csv <- csvFunc()
            if (inherits(plot_data_csv, "list")) {
              plot_data_csv <- plot_data_csv[[1]]
            }
          } else {
            plot_data_csv <- NULL
          }
          for (i in 1:length(plot$x$data)) {
            plot$x$data[[i]]$hovertemplate <- NULL
          }
          for (i in 1:length(plot$x$attrs)) {
            plot$x$attrs[[i]]$hovertemplate <- NULL
          }
          json <- plotly::plotly_json(plot, TRUE) # requires `listviewer` to work properly
          res <- session$registerDataObj(
            "plotly_graph", json$x$data,
            function(data, req) {
              httpResponse(
                status = 200,
                content_type = "application/json",
                content = data
              )
            }
          )
          res2 <- session$registerDataObj(
            "plotly_graph2", jsonlite::toJSON(list(data = as.list(plot_data_csv))),
            function(data, req) {
              httpResponse(
                status = 200,
                content_type = "application/json",
                content = data
              )
            }
          )
          url <- getEditorUrl(session, res, res2)
          tags$iframe(src = url, style = "height: 85vh; width: 100%;")
        })
      }


      ## --------------------------------------------------------------------------------
      ## ------------------------ FIGURE ------------------------------------------------
      ## --------------------------------------------------------------------------------

      ## these engines cannot (yet) provide html
      if (plotlib %in% c("base")) {
        download.fmt <- setdiff(download.fmt, c("html"))
      }

      do.pdf <- "pdf" %in% download.fmt
      do.png <- "png" %in% download.fmt
      do.html <- "html" %in% download.fmt
      do.obj <- "obj" %in% download.fmt
      do.csv <- !is.null(csvFunc)
      do.excel <- !is.null(csvFunc)

      PNGFILE <- PDFFILE <- HTMLFILE <- CSVFILE <- EXCELFILE <- NULL
      if (do.pdf) PDFFILE <- paste0(gsub("file", "plot", tempfile()), ".pdf")
      if (do.png) PNGFILE <- paste0(gsub("file", "plot", tempfile()), ".png")
      if (do.csv) CSVFILE <- paste0(gsub("file", "data", tempfile()), ".csv")
      if (do.excel) EXCELFILE <- paste0(gsub("file", "data", tempfile()), ".xlsx")
      HTMLFILE <- paste0(tempfile(), ".html") ## tempory for webshot
      HTMLFILE
      unlink(HTMLFILE)

      ## ============================================================
      ## =============== Download Handlers ==========================
      ## ============================================================

      if (do.png && is.null(download.png)) {
        download.png <- shiny::downloadHandler(
          filename = paste0(filename, ".png"),
          content = function(file) {
            png.width <- input$pdf_width * 80
            png.height <- input$pdf_height * 80
            resx <- 4 ## upresolution
            shiny::withProgress(
              {
                if (plotlib == "plotly") {
                  p <- func()
                  p$width <- png.width
                  p$height <- png.height
                  plotlyExport(p, PNGFILE, width = p$width, height = p$height, scale = resx)
                } else if (plotlib == "iheatmapr") {
                  p <- func()
                  iheatmapr::save_iheatmap(p, vwidth = png.width, vheight = png.height, PNGFILE)
                } else if (plotlib == "visnetwork") {
                  p <- func()
                  visPrint(p,
                    file = PNGFILE,
                    width = png.width * resx * 2,
                    height = png.height * resx * 2,
                    delay = vis.delay,
                    zoom = 1
                  )
                } else if (plotlib %in% c("htmlwidget", "pairsD3", "scatterD3")) {
                  p <- func()
                  htmlwidgets::saveWidget(p, HTMLFILE)
                  webshot2::webshot(
                    url = HTMLFILE, file = PNGFILE,
                    vwidth = png.width * resx, vheight = png.height * resx
                  )
                } else if (plotlib %in% c("ggplot", "ggplot2")) {
                  ggplot2::ggsave(PNGFILE, plot = func(), dpi = 72 * resx)
                } else if (plotlib == "grid") {
                  p <- func()
                  png(PNGFILE,
                    width = png.width * resx,
                    height = png.height * resx,
                    pointsize = 1.2 * pdf.pointsize,
                    res = 72 * resx
                  )
                  grid::grid.draw(p)
                  dev.off()
                } else if (plotlib == "image") {
                  p <- func()
                  dbg("[downloadHandler.PNG] copy image ", p$src, "to PNGFILE", PNGFILE)
                  file.copy(p$src, PNGFILE, overwrite = TRUE)
                } else if (plotlib == "generic") {
                  ## generic function should produce PNG inside plot func()
                } else if (plotlib == "base") {
                  # Save original plot parameters
                  if (remove_margins == TRUE) {
                    par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
                  }
                  png(PNGFILE,
                    width = png.width * resx,
                    height = png.height * resx,
                    pointsize = 1.2 * pdf.pointsize,
                    res = 72 * resx
                  )
                  print(func())
                  dev.off() ## important!!
                } else if (plotlib == "svgPanZoom") {
                  p <- func()
                  htmlwidgets::saveWidget(p, HTMLFILE)
                  webshot2::webshot(
                    url = HTMLFILE, file = PNGFILE,
                    vwidth = png.width, vheight = png.height, zoom = 4
                  )
                } else { ## end base

                  png(PNGFILE, pointsize = pdf.pointsize)
                  plot.new()
                  mtext("Error. PNG not available.", line = -8)
                  dev.off()
                }

                ## finally copy to final exported file
                dbg("[downloadHandler.PNG] copy PNGFILE", PNGFILE, "to download file", file)
                file.copy(PNGFILE, file, overwrite = TRUE)
                ## ImageMagick or pdftk
                if (TRUE && !add.watermark %in% c("none", FALSE)) {
                  message("[plotModule] adding watermark to PNG...")
                  markfile <- file.path(FILES, "watermark-logo.png")
                  addWatermark.PNG2(file, mark = markfile, position = add.watermark)
                }
                ## Record downloaded plot
                record_plot_download(ns("") %>% substr(1, nchar(.) - 1))
              },
              message = "Exporting to PNG",
              value = 0.8
            )
          } ## content
        ) ## PNG downloadHandler
      } ## end if do.png

      if (do.pdf && is.null(download.pdf)) {
        download.pdf <- shiny::downloadHandler(
          filename = paste0(filename, ".pdf"),
          content = function(file) {
            pdf.width <- input$pdf_width
            pdf.height <- input$pdf_height
            shiny::withProgress(
              {
                if (plotlib == "plotly") {
                  p <- func()
                  p$width <- pdf.width * 80
                  p$height <- pdf.height * 80
                  plotlyExport(p, PDFFILE, width = p$width, height = p$height)
                } else if (plotlib == "iheatmapr") {
                  p <- func()
                  iheatmapr::save_iheatmap(p, vwidth = pdf.width * 80, vheight = pdf.height * 80, PDFFILE)
                } else if (plotlib == "visnetwork") {
                  p <- func()
                  visPrint(p, file = PDFFILE, width = pdf.width, height = pdf.height, delay = vis.delay, zoom = 1)
                } else if (plotlib %in% c("htmlwidget", "pairsD3", "scatterD3")) {
                  p <- func()
                  htmlwidgets::saveWidget(p, HTMLFILE)
                  webshot2::webshot(url = HTMLFILE, file = PDFFILE, vwidth = pdf.width * 100, vheight = pdf.height * 100)
                } else if (plotlib %in% c("ggplot", "ggplot2")) {
                  p <- func()
                  pdf(PDFFILE, width = pdf.width, height = pdf.height, pointsize = pdf.pointsize)
                  print(p)
                  dev.off()
                } else if (plotlib %in% c("grid")) {
                  p <- func()
                  pdf(PDFFILE, width = pdf.width, height = pdf.height, pointsize = pdf.pointsize)
                  grid::grid.draw(p)
                  dev.off()
                } else if (plotlib == "image") {
                  p <- func()
                } else if (plotlib == "generic") {
                  ## generic function should produce PDF inside plot func()
                } else if (plotlib == "base") {
                  pdf(
                    file = PDFFILE, width = pdf.width, height = pdf.height,
                    pointsize = pdf.pointsize
                  )
                  print(func())
                  dev.off() ## important!!
                } else if (plotlib == "svgPanZoom") {
                  p <- func()
                  htmlwidgets::saveWidget(p, HTMLFILE)
                  webshot2::webshot(
                    url = HTMLFILE, file = PDFFILE,
                    vwidth = pdf.width * 100, vheight = pdf.height * 100
                  )
                } else { ## end base
                  pdf(PDFFILE, pointsize = pdf.pointsize)
                  plot.new()
                  mtext("Error. PDF not available.", line = -8)
                  dev.off()
                }

                ## finally copy to final exported file
                dbg("[downloadHandler.PDF] copy PDFFILE", PDFFILE, "to download file", file)
                file.copy(PDFFILE, file, overwrite = TRUE)

                ## ImageMagick or pdftk
                if (TRUE && !add.watermark %in% c(FALSE, "none")) {
                  message("[plotModule] adding watermark to PDF...")
                  markfile <- file.path(FILES, "watermark-logo.pdf")
                  addWatermark.PDF2(file, w = pdf.width, h = pdf.height, mark = markfile)
                }
                # Add settings
                if (input$get_pdf_settings) {
                  addSettings(ns, session, file)
                }
                ## Record downloaded plot
                record_plot_download(ns("") %>% substr(1, nchar(.) - 1))
              },
              message = "Exporting to PDF",
              value = 0.8
            )
          } ## content
        ) ## PDF downloadHandler
      } ## end if do.pdf

      saveHTML <- function() {
        if (plotlib == "plotly") {
          p <- func()
          htmlwidgets::saveWidget(p, HTMLFILE)
        } else if (plotlib %in% c("htmlwidget", "pairsD3", "scatterD3")) {
          p <- func()
          htmlwidgets::saveWidget(p, HTMLFILE)
        } else if (plotlib == "iheatmapr") {
          p <- func()
          iheatmapr::save_iheatmap(p, HTMLFILE)
        } else if (plotlib == "visnetwork") {
          p <- func()
          visNetwork::visSave(p, HTMLFILE)
        } else if (plotlib %in% c("ggplot", "ggplot2")) {
          p <- func()
          htmlwidgets::saveWidget(plotly::ggplotly(p), file = HTMLFILE)
        } else if (plotlib == "image") {
          write("<body>image cannot export to HTML</body>", HTMLFILE)
        } else if (plotlib == "generic") {
          ## generic function should produce PDF inside plot func()
        } else if (plotlib == "base") {
          write("<body>R base plots cannot export to HTML</body>", HTMLFILE)
        } else if (plotlib == "svgPanZoom") {
          p <- func()
          htmlwidgets::saveWidget(p, HTMLFILE)
        } else { ## end base
          write("<body>HTML export error</body>", file = HTMLFILE)
        }
        return(HTMLFILE)
      }

      if (do.html && is.null(download.html)) {
        download.html <- shiny::downloadHandler(
          filename = paste0(filename, ".html"),
          content = function(file) {
            shiny::withProgress(
              {
                if (plotlib == "plotly") {
                  p <- func()
                  htmlwidgets::saveWidget(p, HTMLFILE)
                } else if (plotlib %in% c("htmlwidget", "pairsD3", "scatterD3")) {
                  p <- func()
                  htmlwidgets::saveWidget(p, HTMLFILE)
                } else if (plotlib == "iheatmapr") {
                  p <- func()
                  iheatmapr::save_iheatmap(p, HTMLFILE)
                } else if (plotlib == "visnetwork") {
                  p <- func()
                  visNetwork::visSave(p, HTMLFILE)
                } else if (plotlib %in% c("ggplot", "ggplot2")) {
                  p <- func()
                  htmlwidgets::saveWidget(plotly::ggplotly(p), file = HTMLFILE)
                } else if (plotlib == "generic") {
                  ## generic function should produce PDF inside plot func()
                } else if (plotlib == "image") {
                  write("<body>image cannot be exported to HTML</body>", HTMLFILE)
                } else if (plotlib == "base") {
                  write("<body>R base plots cannot be exported to HTML</body>", HTMLFILE)
                } else { ## end base
                  write("<body>HTML export error</body>", file = HTMLFILE)
                }
                ## Record downloaded plot
                record_plot_download(ns("") %>% substr(1, nchar(.) - 1))
                ## finally copy to fina lexport file
                file.copy(HTMLFILE, file, overwrite = TRUE)
              },
              message = "Exporting to HTML",
              value = 0.8
            )
          } ## end of content
        ) ## end of HTML downloadHandler
      } ## end of do HTML

      if (do.obj) {
        if (plotlib == "plotly") {
          download.obj <- shiny::downloadHandler(
            filename = paste0(filename, ".rds"),
            content = function(file) {
              shiny::withProgress(
                {
                  p <- func()
                  ## we need to strip away unnecessary environment to prevent save bloat
                  b <- plotly::plotly_build(p)$x[c("data", "layout", "config")]
                  b <- plotly::as_widget(b) ## from JSON back to R object
                  saveRDS(b, file = file)
                  ## Record downloaded plot
                  record_plot_download(ns("") %>% substr(1, nchar(.) - 1))
                },
                message = "saving plot object",
                value = 0.2
              )
            } ## end of content
          ) ## end of object downloadHandler
        }
      } ## end of do object

      ## if(do.csv && is.null(download.csv) )  {
      if (do.csv) {
        download.csv <- shiny::downloadHandler(
          filename = paste0(filename, ".csv"),
          content = function(file) {
            shiny::withProgress(
              {
                data <- csvFunc()
                if (is.list(data) && !is.data.frame(data)) data <- data[[1]]
                write.csv(data, file = file)
                ## Record downloaded plot
                record_plot_download(ns("") %>% substr(1, nchar(.) - 1))
              },
              message = "Exporting to CSV",
              value = 0.8
            )
          } ## end of content
        ) ## end of HTML downloadHandler
      } ## end of do HTML

      # Excel download
      if (do.excel) {
        download.excel <- shiny::downloadHandler(
          filename = paste0(filename, ".xlsx"),
          content = function(file) {
            shiny::withProgress({
              data <- csvFunc()
              if (is.list(data) && !is.data.frame(data)) data <- data[[1]]
              openxlsx::write.xlsx(data, file = file)
            })
          }
        )
      }

      ## --------------------------------------------------------------------------------
      ## ------------------------ OUTPUT ------------------------------------------------
      ## --------------------------------------------------------------------------------
      if (is.null(card)) {
        observeEvent(input$downloadOption, {
          if (input$downloadOption == "png") {
            output$download <- download.png
          }
          if (input$downloadOption == "pdf") {
            output$download <- download.pdf
          }
          if (input$downloadOption == "csv") {
            output$download <- download.csv
          }
          if (input$downloadOption == "excel") {
            output$download <- download.excel
          }
          if (input$downloadOption == "html") {
            output$download <- download.html
          }
          if (input$downloadOption == "obj") {
            output$download <- download.obj
          }
        })
      } else {
        observeEvent(input$downloadOption, {
          if (input$downloadOption == "png") {
            output[[paste0(
              "download",
              card
            )]] <- download.png
          }
          if (input$downloadOption == "pdf") {
            output[[paste0(
              "download",
              card
            )]] <- download.pdf
          }
          if (input$downloadOption == "csv") {
            output[[paste0(
              "download",
              card
            )]] <- download.csv
          }
          if (input$downloadOption == "excel") {
            output[[paste0(
              "download",
              card
            )]] <- download.excel
          }
          if (input$downloadOption == "html") {
            output[[paste0(
              "download",
              card
            )]] <- download.html
          }
          if (input$downloadOption == "obj") {
            output[[paste0(
              "download",
              card
            )]] <- download.obj
          }
        })
      }

      ## --------------------------------------------------------------------------------
      ## ---------------------------- UI ------------------------------------------------
      ## --------------------------------------------------------------------------------

      if (is.null(func2)) func2 <- func
      if (is.null(plotlib2)) plotlib2 <- plotlib
      if (length(height) == 1) height <- c(height, 700)
      if (length(width) == 1) width <- c(width, 1200)
      if (length(res) == 1) res <- c(res, 1.3 * res)

      res.1 <- res[1]
      res.2 <- res[2]

      ## width and height should actually be speficied in UI, not here.
      ifnotchar.int <- function(s) ifelse(grepl("[%]$|auto|vmin|vh|vw|vmax", s), s, as.integer(s))
      width.1 <- ifnotchar.int(width[1])
      width.2 <- ifnotchar.int(width[2])
      height.1 <- ifnotchar.int(height[1])
      height.2 <- ifnotchar.int(height[2])

      ## This sets the correct render and output functions for different
      ## plotting libraries.

      getRenderFunc <- function(plotlib) {
        switch(plotlib,
          generic = NULL,
          htmlwidget = NULL,
          plotly = plotly::renderPlotly,
          pairsD3 = pairsD3::renderPairsD3,
          visnetwork = visNetwork::renderVisNetwork,
          ggplot = shiny::renderPlot,
          grid = function(x) shiny::renderPlot(grid::grid.draw(x, recording = FALSE)),
          iheatmapr = iheatmapr::renderIheatmap,
          image = shiny::renderImage,
          base = shiny::renderPlot,
          svgPanZoom = svgPanZoom::renderSvgPanZoom,
          renderUI = shiny::renderUI,
          shiny::renderPlot
        )
      }

      if (is.null(renderFunc)) {
        renderFunc <- getRenderFunc(plotlib)
      }

      if (is.null(renderFunc2)) {
        renderFunc2 <- getRenderFunc(plotlib2)
      }

      render <- render2 <- NULL

      if (!is.null(func) && plotlib == "base") {
        render <- shiny::renderPlot(
          {
            func()
          },
          res = res.1
        )
      }
      if (!is.null(func2) && plotlib2 == "base") {
        render2 <- shiny::renderPlot(
          {
            func2()
          },
          res = res.2
        )
      }
      if (!is.null(func) && plotlib == "grid") {
        render <- shiny::renderPlot(
          {
            grid::grid.draw(func(), recording = FALSE)
          },
          res = res.1
        )
      }
      if (!is.null(func2) && plotlib2 == "grid") {
        render2 <- shiny::renderPlot(
          {
            grid::grid.draw(func2(), recording = FALSE)
          },
          res = res.2
        )
      }
      if (plotlib == "image") {
        render <- shiny::renderImage(func(), deleteFile = FALSE)
      }
      if (!is.null(func2) && plotlib2 == "image") {
        render2 <- shiny::renderImage(func2(), deleteFile = FALSE)
      }

      if (grepl("cacheKeyExpr", head(renderFunc, 1))) {
        render <- shiny::renderCachedPlot(
          func(),
          cacheKeyExpr = {
            list(csvFunc())
          },
          res = res.1
        )
      }
      if (grepl("cacheKeyExpr", head(renderFunc2, 1))) {
        render2 <- shiny::renderCachedPlot(
          func2(),
          cacheKeyExpr = {
            list(csvFunc())
          },
          res = res.2
        )
      }

      if (is.null(render)) {
        if (plotlib == "plotly") {
          # If the plotting function is `plotly`, add the edit button
          render <- renderFunc({
            # By default remove plotly logo from all plots
            plot <- func() %>%
              plotly::config(
                displaylogo = FALSE,
                scrollZoom = TRUE
              ) %>%
              plotly::plotly_build()

            if (remove_margins == TRUE) {
              plot <- plot %>% plotly::layout(margin = list(l = 0, r = 0, t = 0, b = 0))
            }

            # If there is already custom buttons, append the edit one
            # (issue #2210 plotly/plotly.R)
            if (inherits(plot$x$config$modeBarButtons, "list")) {
              for (y in 1:length(plot$x$config$modeBarButtons[[1]])) {
                if (plot$x$config$modeBarButtons[[1]][[y]] == "toImage") {
                  plot$x$config$modeBarButtons[[1]][[y]] <- NULL
                  break
                }
              }
              plot$x$config$modeBarButtons[[1]] <- append(
                plot$x$config$modeBarButtons[[1]],
                list(octocat)
              )
            } else { # Otherwise, apply the button regularly
              plot <- plot %>%
                plotly::config(
                  modeBarButtonsToAdd = list(octocat),
                  modeBarButtonsToRemove = c("zoomIn2d", "toImage")
                )
            }
            plot
          })
        } else {
          render <- renderFunc(func())
        }
      }

      if (is.null(render2) && !is.null(func2)) {
        if (plotlib2 == "plotly") {
          render2 <- renderFunc2({
            # By default remove plotly logo from all plots
            plot <- func2() %>%
              plotly::config(
                displaylogo = FALSE,
                scrollZoom = TRUE
              ) %>%
              plotly::plotly_build()
            # If there is already custom buttons, append the edit one
            # (issue #2210 plotly/plotly.R)
            if (inherits(plot$x$config$modeBarButtons, "list")) {
              for (y in 1:length(plot$x$config$modeBarButtons[[1]])) {
                if (plot$x$config$modeBarButtons[[1]][[y]] == "toImage") {
                  plot$x$config$modeBarButtons[[1]][[y]] <- NULL
                  break
                }
              }
              plot$x$config$modeBarButtons[[1]] <- append(
                plot$x$config$modeBarButtons[[1]],
                list(octocat)
              )
            } else { # Otherwise, apply the button regularly
              plot <- plot %>%
                plotly::config(
                  modeBarButtonsToAdd = list(octocat),
                  modeBarButtonsToRemove = c("zoomIn2d", "toImage")
                )
            }
            plot
          })
        } else {
          render2 <- renderFunc2(func2())
        }
      }

      if (is.null(card)) {
        output$renderfigure <- render
        output$renderpopup <- render2
      } else {
        output[[paste0(
          "renderfigure",
          card
        )]] <- render
        output[[paste0(
          "renderpopup",
          card
        )]] <- render2
      }

      shiny::observeEvent(input$copy_info, {
        shinyjs::runjs(
          paste0(
            "addTick('",
            ns("copy_info"),
            "')"
          )
        )
      })


      ## --------------------------------------------------------------------------------
      ## ---------------------------- RETURN VALUE --------------------------------------
      ## --------------------------------------------------------------------------------

      list(
        plotfun = func,
        plotfun2 = func2,
        .tmpfiles = c(pdf = PDFFILE, html = HTMLFILE),
        render = render,
        render2 = render2,
        download.pdf = download.pdf,
        download.png = download.png,
        download.html = download.html,
        download.csv = download.csv,
        download.excel = download.excel,
        saveHTML = saveHTML,
        renderFunc = renderFunc
      )
    }
  )
}


colBL <- "#00448855"
colRD <- "#88004455"

plotlyExport <- function(p, file = paste0(filename, ".pdf"), format = tools::file_ext(file),
                         width = NULL, height = NULL, scale = 1, server = NULL) {
  is.docker <- file.exists("/.dockerenv")
  is.docker
  export.ok <- FALSE

  if (class(p)[1] != "plotly") {
    message("[plotlyExport] ERROR : not a plotly object")
    return(NULL)
  }
  ## remove old
  unlink(file, force = TRUE)

  ## See if Kaleido is available
  if (1 && !export.ok) {
    ## https://github.com/plotly/plotly.R/issues/2179
    reticulate::py_run_string("import sys")
    err <- try(suppressMessages(plotly::save_image(p, file = file, width = width, height = height, scale = scale)))
    export.ok <- class(err) != "try-error"
    if (export.ok) message("[plotlyExport] --> exported with plotly::save_image() (kaleido)")
    export.ok <- TRUE
  }
  if (1 && !export.ok) {
    ## works only for non-GL plots
    err <- try(plotly::export(p, file, width = width, height = height))
    export.ok <- class(err) != "try-error"
    if (export.ok) message("[plotlyExport] --> exported with plotly::export() (deprecated)")
  }
  if (0 && !export.ok) {
    tmp <- paste0(tempfile(), ".html")
    htmlwidgets::saveWidget(p, tmp)
    err <- try(webshot2::webshot(url = tmp, file = file, vwidth = width * 100, vheight = height * 100))
    export.ok <- class(err) != "try-error"
    if (export.ok) message("[plotlyExport] --> exported with webshot2::webshot()")
  }
  if (!export.ok) {
    message("[plotlyExport] WARNING: export failed!")
    if (format == "png") png(file)
    if (format == "pdf") pdf(file)
    par(mfrow = c(1, 1))
    frame()
    text(0.5, 0.5, "Plotly export error", cex = 2)
    dev.off()
  }

  message("[plotlyExport] file.exists(file)=", file.exists(file))
  export.ok <- export.ok && file.exists(file)
  return(export.ok)
}



## ================================================================================
## ====================== END OF FILE +============================================
## ================================================================================
