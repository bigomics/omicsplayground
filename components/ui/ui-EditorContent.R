##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

getEditorContent <- function(plot_type = "default", ns, ns_parent, title, cards = FALSE, outputFunc = NULL, width.2 = NULL, height.2 = NULL) {
  # Default editor content
  default_content <- shiny::div(
    class = "popup-modal",
    modalUI(
      id = ns("plotPopup2"),
      title = title,
      size = "fullscreen",
      footer = NULL,
      bslib::layout_column_wrap(
        style = bslib::css(grid_template_columns = "1fr 5fr"),
        bslib::accordion(
          id = ns("plot_options_accordion"),
          # Plot Type & Color Options
          bslib::accordion_panel(
            "Basic Options",
            shiny::textInput(
              ns_parent("title"),
              "Title",
              value = NULL
            )
          ),
          bslib::accordion_panel(
            "Color Scheme",
            bslib::layout_column_wrap(
              width = 1/2,
              colourpicker::colourInput(
                ns_parent("color_up"), "Up",
                "#f23451"
              ),
              colourpicker::colourInput(
                ns_parent("color_down"), "Down",
                "#3181de"
              )
            )
          ),

          # Axis Options
          bslib::accordion_panel(
            "Text sizes",
            bslib::layout_column_wrap(
              width = 1/2,
              numericInput(ns_parent("label_size"), "Labels", value = 4),
              numericInput(ns_parent("marker_size"), "Points", value = 1),
              numericInput(ns_parent("axis_text_size"), "Axis text", value = 14)
            )
          ),

          # Additional Settings
          bslib::accordion_panel(
            "Labels",
            checkboxInput(ns_parent("color_selection"), "Color just selection", value = FALSE),
            checkboxInput(ns_parent("custom_labels"), "Custom labels", value = FALSE),
            textAreaInput(ns_parent("label_features"), "Label features", value = "")
          )
        ),
        shiny::div(
          class = "popup-plot",
          if (cards) {
            outputFunc[[2]](ns("renderfigure_2"), width = width.2, height = height.2, click = ns("plot_click")) %>%
              bigLoaders::useSpinner()
          } else {
            # outputFunc(ns("renderfigure_2")) %>%
            #   bigLoaders::useSpinner()
          }
        )
      )
    )
  )

  # Return content based on plot type
  switch(plot_type,
    "default" = default_content
  )
} 