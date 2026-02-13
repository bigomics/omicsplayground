##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

getEditorContent <- function(plot_type = "volcano", ns, ns_parent, title, cards = FALSE, outputFunc = NULL, width.2 = NULL, height.2 = NULL, bar_color_default = "#3181de", palette_default = "muted_light") {
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
              width = 1 / 2,
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
              width = 1 / 2,
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
            textAreaInput(ns_parent("label_features"), "Label features", value = ""),
            bslib::layout_column_wrap(
              width = 1 / 2,
              numericInput(ns_parent("box_padding"), "Box padding", value = 0.1, min = 0, step = 0.05),
              numericInput(ns_parent("min_segment_length"), "Min segment", value = 0, min = 0, step = 0.1)
            ),
            shiny::helpText("Box padding: distance from labels to points. Min segment: 0 forces lines to always appear."),
            bslib::layout_column_wrap(
              width = 1 / 2,
              checkboxInput(ns_parent("label_box"), "Box around labels", value = TRUE),
              selectInput(ns_parent("segment_linetype"), "Line type", choices = 1:6, selected = 1)
            ),
            shiny::helpText("Line types: 1=solid, 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash")
          ),
          # ggprism Theme
          bslib::accordion_panel(
            "Prism Theme",
            checkboxInput(ns_parent("use_ggprism"), "Use ggprism theme", value = FALSE),
            conditionalPanel(
              condition = "input.use_ggprism",
              ns = ns_parent,
              selectInput(
                ns_parent("ggprism_palette"),
                "Theme palette",
                choices = c(
                  "Black & White" = "black_and_white",
                  "Colorblind Safe" = "colorblind_safe",
                  "Office" = "office",
                  "Floral" = "floral",
                  "Earth Tones" = "earth_tones",
                  "Pearl" = "pearl",
                  "Muted Rainbow" = "muted_rainbow",
                  "Candy Bright" = "candy_bright",
                  "Prism Dark" = "prism_dark",
                  "Prism Light" = "prism_light",
                  "Winter Soft" = "winter_soft",
                  "Starry" = "starry",
                  "Viridis" = "viridis",
                  "Plasma" = "plasma",
                  "Inferno" = "inferno",
                  "Magma" = "magma"
                ),
                selected = "black_and_white"
              ),
              checkboxInput(ns_parent("ggprism_border"), "Add border", value = FALSE),
              checkboxInput(ns_parent("ggprism_colors"), "Use prism colors", value = FALSE),
              shiny::hr(),
              shiny::tags$label("Axis guides"),
              selectInput(
                ns_parent("ggprism_axis_guide"),
                NULL,
                choices = c(
                  "Default" = "default",
                  "Minor ticks" = "prism_minor",
                  "Offset axis" = "prism_offset",
                  "Offset + minor ticks" = "prism_offset_minor"
                ),
                selected = "default"
              ),
              shiny::hr(),
              checkboxInput(ns_parent("ggprism_show_legend"), "Show legend", value = FALSE),
              conditionalPanel(
                condition = "input.ggprism_show_legend",
                ns = ns_parent,
                bslib::layout_column_wrap(
                  width = 1 / 2,
                  numericInput(ns_parent("ggprism_legend_x"), "X position", value = 0.95, min = 0, max = 1, step = 0.05),
                  numericInput(ns_parent("ggprism_legend_y"), "Y position", value = 0.95, min = 0, max = 1, step = 0.05)
                ),
                shiny::helpText("Position: 0 = left/bottom, 1 = right/top"),
                checkboxInput(ns_parent("ggprism_legend_border"), "Legend border", value = FALSE)
              )
            )
          ),

          # Hyperbolic Cutoff
          bslib::accordion_panel(
            "Significance Cutoff",
            shiny::radioButtons(
              ns_parent("cutoff_type"),
              "Cutoff type:",
              choices = c(
                "Rectangular (traditional)" = "rectangular",
                "Hyperbolic" = "hyperbolic"
              ),
              selected = "rectangular"
            ),
            shiny::conditionalPanel(
              condition = paste0("input['", ns_parent("cutoff_type"), "'] == 'hyperbolic'"),
              shiny::numericInput(
                ns_parent("hyperbola_k"),
                "Curvature (k):",
                value = 1,
                min = 0.1,
                max = 10,
                step = 0.1
              ),
              shiny::helpText("Smaller k = tighter curve (more stringent)")
            )
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
              width = 1 / 2,
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
              width = 1 / 2,
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
          # Group Order (for split heatmaps)
          bslib::accordion_panel(
            "Group Order",
            shiny::uiOutput(ns_parent("hm_group_order_ui"))
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
                bar_color_default
              )
            )
          ),

          # Bars Order
          bslib::accordion_panel(
            "Bars Order",
            shiny::selectInput(
              ns_parent("bars_order"),
              "Sort bars by:",
              choices = c(
                "Alphabetical" = "alphabetical",
                "Value (ascending)" = "ascending",
                "Value (descending)" = "descending",
                "Custom (drag & drop)" = "custom"
              ),
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

  # Scatterplot specific content
  scatterplot_content <- shiny::div(
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
                ns_parent("scatter_color"), "Point Color",
                "#3181de"
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

  # Feature map specific content
  featuremap_content <- shiny::div(
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
              width = 1 / 2,
              colourpicker::colourInput(
                ns_parent("color_low"), "Low color",
                "#3181de"
              ),
              colourpicker::colourInput(
                ns_parent("color_high"), "High color",
                "#f23451"
              )
            )
          ),
          bslib::accordion_panel(
            "Labels",
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
            outputFunc(ns("renderfigure_2")) %>%
              bigLoaders::useSpinner()
          }
        )
      )
    )
  )

  # Enrichment plot specific content
  enrichment_content <- shiny::div(
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
              width = 1 / 2,
              colourpicker::colourInput(
                ns_parent("color_up"), "Up (positive)",
                "#f23451"
              ),
              colourpicker::colourInput(
                ns_parent("color_down"), "Down (negative)",
                "#3181de"
              )
            ),
            colourpicker::colourInput(
              ns_parent("color_line"), "Enrichment line",
              "#00EE00"
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

  # Clustering scatterplot (categorical palette) specific content
  clustering_content <- shiny::div(
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
            shiny::selectInput(
              ns_parent("palette"), "Color palette",
              choices = c(
                "original",
                "muted_light", "default", "light", "dark",
                "super_light", "super_dark", "muted", "expanded",
                "highlight_blue", "highlight_red", "highlight_orange",
                "custom"
              ),
              selected = palette_default
            ),
            shiny::uiOutput(ns_parent("custom_palette_ui"))
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

  # Grouped barplot: palette selector + bar ordering
  grouped_barplot_content <- shiny::div(
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
            shiny::selectInput(
              ns_parent("palette"), "Color palette",
              choices = c(
                "original",
                "muted_light", "default", "light", "dark",
                "super_light", "super_dark", "muted", "expanded",
                "highlight_blue", "highlight_red", "highlight_orange",
                "custom"
              ),
              selected = palette_default
            ),
            shiny::uiOutput(ns_parent("custom_palette_ui"))
          ),
          bslib::accordion_panel(
            "Bars Order",
            shiny::selectInput(
              ns_parent("bars_order"),
              "Sort bars by:",
              choices = c(
                "Alphabetical" = "alphabetical",
                "Value (ascending)" = "ascending",
                "Value (descending)" = "descending",
                "Custom (drag & drop)" = "custom"
              ),
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

  # Gradient scatterplot: low/high color pickers only (no labels)
  gradient_content <- shiny::div(
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
              width = 1 / 2,
              colourpicker::colourInput(
                ns_parent("color_low"), "Low color",
                "#3181de"
              ),
              colourpicker::colourInput(
                ns_parent("color_high"), "High color",
                "#f23451"
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

  # Significance scatter: colors for significance categories
  significance_content <- shiny::div(
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
            colourpicker::colourInput(
              ns_parent("color_both"), "Significant in both",
              "#5B9B5B"
            ),
            colourpicker::colourInput(
              ns_parent("color_one"), "Significant in one",
              "#e3a45a"
            ),
            colourpicker::colourInput(
              ns_parent("color_ns"), "Not significant",
              "#eeeeee"
            )
          ),
          bslib::accordion_panel(
            "Labels",
            checkboxInput(ns_parent("custom_labels"), "Custom labels", value = FALSE),
            textAreaInput(ns_parent("label_features"), "Label features", value = "")
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

  # Scatter with highlight: point color, highlight color, and labels
  scatter_highlight_content <- shiny::div(
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
            colourpicker::colourInput(
              ns_parent("color_point"), "Point color",
              "#222222"
            ),
            colourpicker::colourInput(
              ns_parent("color_highlight"), "Highlight color",
              "#f23451"
            )
          ),
          bslib::accordion_panel(
            "Labels",
            checkboxInput(ns_parent("custom_labels"), "Custom labels", value = FALSE),
            textAreaInput(ns_parent("label_features"), "Label features", value = "")
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

  # Rank/density plot: fill, outline, and highlight colors
  rank_plot_content <- shiny::div(
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
            colourpicker::colourInput(
              ns_parent("color_fill"), "Fill color",
              "#b8d4f0"
            ),
            colourpicker::colourInput(
              ns_parent("color_line"), "Line color",
              "#3181de"
            ),
            colourpicker::colourInput(
              ns_parent("color_highlight"), "Highlight color",
              "#e3a45a"
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
    "barplot" = barplot_content,
    "correlation" = barplot_content,
    "scatterplot" = scatterplot_content,
    "featuremap" = featuremap_content,
    "enrichment" = enrichment_content,
    "clustering" = clustering_content,
    "grouped_barplot" = grouped_barplot_content,
    "gradient" = gradient_content,
    "significance" = significance_content,
    "scatter_highlight" = scatter_highlight_content,
    "rank_plot" = rank_plot_content
  )
}
