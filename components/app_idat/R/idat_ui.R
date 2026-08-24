##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

idat_ui <- function(id) {
  ns <- shiny::NS(id) ## namespace

  title <- div("IDAT Converter", style = "font-size: 18px;")

  ui <- bslib::page_fillable(
    padding = 0,
    div(class = "navbar navbar-static-top", div(title, class = "container-fluid"),
      style = "margin-top: 24px;"),
    bslib::layout_columns(
      col_widths = c(2, 10),
      class = "p-3",
      gap = "3.5rem",
      height = "calc(100vh - 72px)",
      bslib::layout_columns(
        col_widths = 12,
        fill = FALSE,
        shiny::fileInput(ns("idat"),
          shiny::tagList("IDAT files (or a ZIP):", idat_info("upload")),
          multiple = TRUE,
          accept = c(".idat", ".gz", ".zip", ".csv")
        ),
        shiny::div(
          shiny::HTML(paste0(
            "Both <b>_Grn.idat</b> and <b>_Red.idat</b> per sample, plus a ",
            "sample sheet CSV if you have one. Uploads are capped at ",
            IDAT_MAX_UPLOAD_MB, " MB and a 450k pair is ~16 MB, so split ",
            "large cohorts."
          )),
          style = "font-size: 11px; color: #888; margin-top: -12px;"
        ),
        shiny::uiOutput(ns("files_summary")),
        shiny::selectInput(ns("method"),
          shiny::tagList("Preprocessing:", idat_info("preprocessing")),
          choices = idat_normalization_choices(),
          selected = "noob"
        ),
        shiny::div(
          paste(
            "Beyond this preprocessing step, no further normalization is",
            "applied to the beta values."
          ),
          style = "font-size: 11px; color: #888; margin-top: -12px; margin-bottom: 12px;"
        ),
        shiny::selectInput(ns("platform"),
          shiny::tagList("Array:", idat_info("array")),
          choices = idat_platform_choices(),
          selected = "auto"
        ),
        bslib::accordion(
          open = FALSE,
          class = "mb-3",
          bslib::accordion_panel(
            "Advanced",
            shiny::numericInput(ns("detect_p"),
              shiny::tagList("Detection p-value cut:", idat_info("detect_p")),
              value = 0.01, min = 0, max = 1, step = 0.005
            ),
            ## .form-text is Bootstrap's help-text class - it already sits
            ## under its input. The negative margin this replaced pulled the
            ## text up over the field.
            shiny::div("Probes above this become NA in that sample.",
              class = "form-text mb-3"
            ),
            shiny::numericInput(ns("max_fail"),
              shiny::tagList("Drop a probe failing in >", idat_info("max_fail")),
              value = 0.05, min = 0, max = 1, step = 0.01
            ),
            shiny::div(
              paste(
                "Fraction of samples. Raising it drops more probes, which",
                "costs epigenetic-clock coverage. Sex-chromosome probes are",
                "never dropped."
              ),
              class = "form-text"
            )
          )
        ),
        div(
          style = "display: flex; flex-direction: column; gap: 0;",
          shiny::uiOutput(ns("progress")),
          bslib::input_task_button(ns("convert"),
            shiny::tagList("Convert", idat_info("convert")),
            icon = icon("arrows-rotate"),
            label_busy = "Converting...",
            class = "btn-primary mb-2",
            width = "100%"
          ),
          shiny::uiOutput(ns("download_ui"))
        )
      ),
      bslib::layout_columns(
        col_widths = 12,
        div(
          style = "padding-left: 45px; height: 100%;",
          shiny::uiOutput(ns("table_area"), style = "height: 100%;")
        )
      )
    )
  )

  return(ui)
}
