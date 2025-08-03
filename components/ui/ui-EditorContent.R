##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

getEditorContent <- function(plot_type = "volcano", ns, ns_parent, title, cards = FALSE, outputFunc = NULL, width.2 = NULL, height.2 = NULL) {

  # Default editor content
  volcano_content <- shiny::div(
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

          bslib::accordion_panel(
            "Margins",
            checkboxInput(ns_parent("margin_checkbox"), "Custom margins", value = FALSE),
            conditionalPanel(
              condition = "input.margin_checkbox",
              ns = ns_parent,
              numericInput(ns_parent("margin_left"), "Left", value = 10),
              numericInput(ns_parent("margin_right"), "Right", value = 10),
              numericInput(ns_parent("margin_top"), "Top", value = 10),
              numericInput(ns_parent("margin_bottom"), "Bottom", value = 10)
            )
          ),

          bslib::accordion_panel(
            "Aspect Ratio",
            checkboxInput(ns_parent("aspect_ratio_checkbox"), "Custom aspect ratio", value = FALSE),
            conditionalPanel(
              condition = "input.aspect_ratio_checkbox",
              numericInput(ns_parent("aspect_ratio"), NULL, value = 0.5, min = 0.1, max = 10),
              ns = ns_parent
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
          }
        )
      )
    )
  )

  # Heatmap specific content
  heatmap_content <- shiny::div(
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
          # Text Sizes
          bslib::accordion_panel(
            "Labels",
            bslib::layout_column_wrap(
              width = 1/2,
              numericInput(ns_parent("label_size"), "Label size:", value = 10),
              numericInput(ns_parent("annot_cex"), "Annotation size:", value = 12)
            ),
            shiny::numericInput(
              ns_parent("column_names_rot"),
              "Column names rotation",
              value = 45,
              min = 0,
              max = 90
            ),
            shiny::numericInput(
              ns_parent("rownames_width"),
              "Row names width",
              value = 40,
              min = 10,
              max = 200
            )            
          ),
          # Color Scheme
          bslib::accordion_panel(
            "Color Scheme",
            bslib::layout_column_wrap(
              width = 1/2,
              colourpicker::colourInput(
                ns_parent("color_high"), "High",
                "#f23451"
              ),
              colourpicker::colourInput(
                ns_parent("color_mid"), "Mid",
                "#eeeeee"
              ),
              colourpicker::colourInput(
                ns_parent("color_low"), "Low",
                "#3181de"
              )
            )
          ),
          # Margins
          bslib::accordion_panel(
            "Margins",
            checkboxInput(ns_parent("margin_checkbox"), "Custom margins", value = FALSE),
            conditionalPanel(
              condition = "input.margin_checkbox==true",
              ns = ns_parent,
              numericInput(ns_parent("margin_left"), "Left", value = 10),
              numericInput(ns_parent("margin_right"), "Right", value = 10),
              numericInput(ns_parent("margin_top"), "Top", value = 10),
              numericInput(ns_parent("margin_bottom"), "Bottom", value = 10)
            )
          )
        ),
        shiny::div(
          class = "popup-plot",
          if (cards) {
            outputFunc[[2]](ns("renderfigure_2"), width = width.2, height = height.2) %>%
              bigLoaders::useSpinner()
          }
        )
      )
    )
  )

  # Barplot specific content
  barplot_content <- shiny::div(
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
          # Color Scheme
          bslib::accordion_panel(
            "Color Scheme",
            bslib::layout_column_wrap(
              width = 1,
              colourpicker::colourInput(
                ns_parent("bar_color"), "Bar Color",
                "#3181de"
              )
            )
          ),

          # Bars Order
          bslib::accordion_panel(
            "Bars Order",
            shiny::selectInput(
              ns_parent("bars_order"),
              "Sort bars by:",
              choices = c("Alphabetical" = "alphabetical",
                          "Value (ascending)" = "ascending",
                          "Value (descending)" = "descending",
                          "Custom (drag & drop)" = "custom"),
              selected = "alphabetical"
            ),
            shiny::conditionalPanel(
              condition = paste0("input['", ns_parent("bars_order"), "'] == 'custom'"),
              shiny::div(
                shiny::uiOutput(ns_parent("rank_list"))
              )
            )
          )
        ),
        shiny::div(
          class = "popup-plot",
          if (cards) {
            outputFunc[[2]](ns("renderfigure_2"), width = width.2, height = height.2) %>%
              bigLoaders::useSpinner()
          } else {
            outputFunc(ns("renderfigure_2")) %>%
              bigLoaders::useSpinner()
          }
        )
      )
    )
  )

  # Return content based on plot type
  switch(plot_type,
    "volcano" = volcano_content,
    "heatmap" = heatmap_content,
    "barplot" = barplot_content
  )
} 
