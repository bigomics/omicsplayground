---
name: plot-editor
description: Guide for adding or modifying plot editors in Omics Playground. Use when working on PlotModuleUI, getEditorContent, editor accordion panels, or any plot customization UI involving ns_parent scoping.
---

# Plot Editor System — Omics Playground

## Architecture Overview

The plot editor lets users customize static/ggplot plots (colors, labels, margins, themes) via a fullscreen accordion-panel modal. It uses a **4-layer system**:

```
PlotModuleUI (integration layer)
  → getEditorContent() (editor UI per plot type)
    → render functions in board modules (read editor inputs)
      → playbase plotting functions (accept editor params)
```

### How it works

1. A board module calls `PlotModuleUI()` with `editor = TRUE` and a `plot_type` (one of: `"volcano"`, `"heatmap"`, `"barplot"`, `"correlation"`, `"scatterplot"`, `"featuremap"`, `"enrichment"`, `"clustering"`, `"expression_barplot"`, `"rank_plot"`).
2. `PlotModuleUI` renders a pencil icon button that opens a fullscreen modal.
3. Inside the modal, `getEditorContent()` returns a two-column layout: accordion panels (left) + live plot preview (right).
4. All editor input IDs use `ns_parent()` (the **parent** module's namespace), so the parent's server function reads them directly via `input$<name>`.
5. Editor inputs only affect the **static/ggplot card** (card 2), not the interactive plotly card (card 1).

### Key insight: `ns_parent` scoping

Editor inputs are created with `ns_parent("color_up")` instead of `ns("color_up")`. This places them in the **calling module's** namespace, not inside the `PlotModule` sub-module. The parent's `moduleServer` then reads them as plain `input$color_up`.

## Key Files

| File | Purpose |
|------|---------|
| `components/ui/ui-PlotModule.R` | `PlotModuleUI()` / `PlotModuleServer()` — integration layer |
| `components/ui/ui-EditorContent.R` | `getEditorContent()` — editor UI definitions per plot type |
| `components/ui/ui-ColorDefaults.R` | Global color theme singleton — `COLOR_THEME_DEFAULTS`, `get_color_theme()` |
| `components/board.expression/R/expression_plot_volcano.R` | Volcano editor implementation (full-featured) |
| `components/board.clustering/R/clustering_plot_splitmap.R` | Heatmap/splitmap editor implementation |
| `components/board.dataview/R/dataview_plot_expression.R` | Barplot editor implementation (simplest) |

The plotting functions themselves live in the `playbase` package (e.g. `playbase::ggVolcano()`, `playbase::gx.splitmap()`).

## Global Color Theme System

A singleton reactive store (`ui-ColorDefaults.R`) lets users set **7 semantic theme colors** that automatically propagate to all editor colour pickers across every plot module.

### Singleton API

```r
# Initialise once in server.R
init_color_theme()

# Read theme values anywhere (reactive context)
ct <- get_color_theme()
ct$primary    # "#f23451" — Up / High
ct$secondary  # "#3181de" — Down / Low / scatter
ct$neutral    # "#eeeeee" — Mid / NS
ct$bar_color  # "#A6CEE3" — Bar color
ct$accent     # "#e3a45a" — Significant in one
ct$success    # "#5B9B5B" — Significant in both
ct$line       # "#00EE00" — Enrichment line
ct$palette    # "default" — default palette name
```

### Mapping: theme keys → editor input IDs

`COLOR_THEME_MAPPING` in `ui-ColorDefaults.R` declares which editor inputs each theme key controls:

```r
COLOR_THEME_MAPPING <- list(
  primary   = c("color_up", "color_high"),
  secondary = c("color_down", "color_low", "scatter_color", "rank_color_line"),
  neutral   = c("color_mid", "color_ns"),
  bar_color = c("bar_color"),
  accent    = c("color_one"),
  success   = c("color_both"),
  line      = c("color_line")
)
```

`PlotModuleServer` observes each theme key and calls `updateColourInput(parent_session, input_id, value)` on change. Inputs not present in a given module are silently ignored.

### Theme-aware editor color defaults

`getEditorContent()` snapshots the current theme at UI construction time so lazily-loaded modules start with whatever the user already set:

```r
# Inside getEditorContent() — done once at the top
ct <- shiny::isolate(shiny::reactiveValuesToList(get_color_theme()))

# Then use ct$primary, ct$secondary, etc. instead of hardcoded hex:
colourpicker::colourInput(ns_parent("color_up"),   "Up",   ct$primary)
colourpicker::colourInput(ns_parent("color_down"), "Down", ct$secondary)
```

**Never hardcode hex in editor colour pickers.** Always reference `ct$<key>` so the picker starts at the user's current theme.

### Null-safe theme defaults in render functions

In board module render functions, always fall back to the theme singleton when the editor input hasn't been set yet:

```r
clr_up   <- if (!is.null(input$color_up))   input$color_up   else get_color_theme()$primary
clr_down <- if (!is.null(input$color_down)) input$color_down else get_color_theme()$secondary
```

### Input name collision warning

**Do not reuse a name already in `COLOR_THEME_MAPPING`.** For example, `color_line` is mapped to the enrichment green (`theme$line`). A rank plot that also needs a "line color" must use a distinct name like `rank_color_line` — and add that name to the appropriate mapping key (`secondary` in this case).

## Step-by-Step: Adding a New Editor

### Step 1: Add a new case in `getEditorContent()`

In `components/ui/ui-EditorContent.R`, add a new content block and a new case to the `switch` at the bottom:

```r
# In getEditorContent():

my_new_content <- shiny::div(
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
        # Add accordion panels here (see templates below)
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

# Add to the switch statement:
switch(plot_type,
  "volcano" = volcano_content,
  "heatmap" = heatmap_content,
  "barplot" = barplot_content,
  "my_new_type" = my_new_content    # <-- add here
)
```

### Step 2: Define accordion panels using `ns_parent()`

All input IDs **must** use `ns_parent()`, not `ns()`. Use `ct$<key>` (from the theme snapshot) for colour picker defaults:

```r
bslib::accordion_panel(
  "Color Scheme",
  colourpicker::colourInput(ns_parent("my_color"), "Color", ct$secondary)
),
bslib::accordion_panel(
  "Text Sizes",
  numericInput(ns_parent("my_label_size"), "Label size", value = 12)
)
```

### Step 3: Enable the editor in the board module's UI function

In your board module's UI, call `PlotModuleUI` with:

```r
PlotModuleUI(
  ns("pltmod"),
  title = title,
  plotlib = c("plotly", "ggplot"),   # card 1 = plotly, card 2 = ggplot
  cards = TRUE,
  card_names = c("dynamic", "static"),
  editor = TRUE,                     # <-- enables pencil button
  ns_parent = ns,                    # <-- passes parent namespace
  plot_type = "my_new_type",         # <-- selects your editor
  # ... other args
)
```

### Step 4: Read editor inputs in the render function

In the board module's `moduleServer`, define the render function for card 2 (static/ggplot). Read editor values from `input$` directly (they are in the parent's namespace thanks to `ns_parent`):

```r
base.RENDER <- function() {
  pd <- plot_data()
  shiny::req(pd)

  # Null-safe — fall back to theme singleton
  my_color <- if (is.null(input$my_color)) get_color_theme()$secondary else input$my_color
  my_label_size <- if (is.null(input$my_label_size) || is.na(input$my_label_size)) 12 else input$my_label_size

  # Pass to plotting function
  playbase::myPlotFunction(
    data = pd,
    color = my_color,
    label.cex = my_label_size
  )
}
```

### Step 5: Wire up PlotModuleServer with the dual-card pattern

```r
plot_grid <- list(
  list(plotlib = "plotly", func = plotly.RENDER, func2 = modal_plotly.RENDER, card = 1),
  list(plotlib = "ggplot", func = base.RENDER, func2 = base.RENDER.modal, card = 2)
)

lapply(plot_grid, function(x) {
  PlotModuleServer(
    "pltmod",
    plotlib = x$plotlib,
    func = x$func,
    func2 = x$func2,
    csvFunc = plot_data_csv,
    card = x$card,
    parent_session = session   # needed for theme observers + click-to-label
  )
})
```

### Step 6: Add parameters to the playbase plotting function

In the `playbase` package, add corresponding parameters to your plotting function. Follow the pattern of `playbase::ggVolcano()` which accepts `colors`, `label.cex`, `marker.size`, `use_ggprism`, etc.

## Code Templates

### Editor UI template (accordion panels)

```r
# In getEditorContent() — components/ui/ui-EditorContent.R
# ct is already available from the snapshot at the top of getEditorContent()
my_content <- shiny::div(
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
            colourpicker::colourInput(ns_parent("color_a"), "Color A", ct$primary),
            colourpicker::colourInput(ns_parent("color_b"), "Color B", ct$secondary)
          )
        ),

        bslib::accordion_panel(
          "Text Sizes",
          bslib::layout_column_wrap(
            width = 1 / 2,
            numericInput(ns_parent("label_size"), "Labels", value = 10),
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
```

### Render function template (null-safe input reading)

```r
base.RENDER <- function() {
  pd <- plot_data()
  shiny::req(pd)

  # Null-safe pattern: fall back to theme singleton for colors, check NA for numerics
  color_a <- if (is.null(input$color_a)) get_color_theme()$primary   else input$color_a
  color_b <- if (is.null(input$color_b)) get_color_theme()$secondary else input$color_b
  label_size      <- if (is.null(input$label_size)      || is.na(input$label_size))      10 else input$label_size
  axis_text_size  <- if (is.null(input$axis_text_size)  || is.na(input$axis_text_size))  14 else input$axis_text_size

  p <- playbase::myPlotFunction(
    data = pd,
    color = color_a,
    label.cex = label_size,
    axis.text.size = axis_text_size
  )

  # Apply margins if enabled
  if (isTRUE(input$margin_checkbox)) {
    margin_top    <- ifelse(is.na(input$margin_top),    10, input$margin_top)
    margin_right  <- ifelse(is.na(input$margin_right),  10, input$margin_right)
    margin_bottom <- ifelse(is.na(input$margin_bottom), 10, input$margin_bottom)
    margin_left   <- ifelse(is.na(input$margin_left),   10, input$margin_left)
    p <- p + ggplot2::theme(
      plot.margin = ggplot2::margin(margin_top, margin_right, margin_bottom, margin_left)
    )
  }

  p
}
```

### PlotModuleUI call template

```r
PlotModuleUI(
  ns("pltmod"),
  title = title,
  label = label,
  caption = caption,
  plotlib = c("plotly", "ggplot"),
  info.text = info.text,
  download.fmt = c("png", "pdf", "csv", "svg"),
  width = width,
  height = height,
  cards = TRUE,
  card_names = c("dynamic", "static"),
  editor = TRUE,
  ns_parent = ns,
  plot_type = "my_new_type"
)
```

### PlotModuleServer call template — dual-card

```r
plot_grid <- list(
  list(plotlib = "plotly", func = plotly.RENDER, func2 = modal_plotly.RENDER, card = 1),
  list(plotlib = "ggplot", func = base.RENDER, func2 = base.RENDER.modal, card = 2)
)

lapply(plot_grid, function(x) {
  PlotModuleServer(
    "pltmod",
    plotlib = x$plotlib,
    func = x$func,
    func2 = x$func2,
    csvFunc = plot_data_csv,
    res = c(80, 95),
    pdf.width = 10,
    pdf.height = 8,
    add.watermark = watermark,
    card = x$card,
    parent_session = session
  )
})
```

### PlotModuleServer call template — single-card

```r
PlotModuleServer(
  "pltmod",
  plotlib = "plotly",
  func = plotly.RENDER,
  func2 = modal_plotly.RENDER,
  csvFunc = plot_data,
  download.fmt = c("png", "pdf", "csv", "svg"),
  res = c(90, 170),
  pdf.width = 6,
  pdf.height = 6,
  add.watermark = watermark,
  parent_session = session   # required for theme observers
)
```

## Configurable Default Parameters

`getEditorContent()` and `PlotModuleUI()` accept parameters that let each plot customize editor defaults:

| Parameter | Default | Purpose | Notes |
|-----------|---------|---------|-------|
| `palette_default` | `"default"` | Initial palette selection in clustering/barplot dropdowns | Set per-plot; `"default"` maps to `muted_light` in most clustering plots |

### `bar_color_input_id` — semantic bar color routing

Inside `getEditorContent()`, the bar colour picker's input ID is chosen based on `plot_type`:

```r
bar_color_input_id <- if (plot_type %in% c("correlation", "expression_barplot")) "scatter_color" else "bar_color"
bar_color_init     <- if (plot_type %in% c("correlation", "expression_barplot")) ct$secondary    else ct$bar_color
```

- **`scatter_color`** (→ `theme$secondary`): expression barplots, correlation plots — semantically a "directional secondary" color
- **`bar_color`** (→ `theme$bar_color`): standard barplots — neutral bar fill

Using `scatter_color` for expression/correlation plots means the theme's secondary color propagates to them automatically.

## Plotly-Only Editor Pattern

Not all editors follow the dual-card pattern. Several editors work with **plotly-only** single-card plots. In this case:

- The editor modal still shows a live preview via `renderfigure_2`
- Editor inputs are read directly in the plotly render function (not a separate ggplot function)
- `PlotModuleServer` is called once (not via `lapply` over a `plot_grid`)
- The `output$renderfigure_2` output uses `func()` inside `plotly::renderPlotly({...})`, which IS a reactive context — so `input$` reads inside `func()` trigger re-rendering

**Key gotcha with plotly color mapping:**

When overriding colors in plotly categorical scatterplots, you **must** use the proper plotly pattern:

```r
# CORRECT — colors actually change when palette is updated
plotly::add_trace(color = pheno_factor, colors = color_vector, ...)

# WRONG — plotly uses hex strings as group labels, NOT display colors
plotly::add_trace(color = hex_color_vector, ...)
```

The `color` parameter sets the grouping variable (a factor). The `colors` parameter sets the actual display colors (a character vector of hex values, one per factor level).

### Forwarding editor colors to dynamic plotly renders (volcano)

For **dual-card volcano plots**, the dynamic plotly card normally ignores editor inputs. To keep the dynamic and static renders in sync, pass the editor color inputs explicitly to `plotlyVolcano()` / `plotlyVolcano_multi()`:

```r
plotly_render.RENDER <- function() {
  col_up   <- if (!is.null(input$color_up))   input$color_up   else get_color_theme()$primary
  col_down <- if (!is.null(input$color_down)) input$color_down else get_color_theme()$secondary

  playbase::plotlyVolcano(
    ...,
    color_up_down = TRUE,
    colors = c(up = col_up, notsig = "#707070AA", down = col_down, notsel = "#cccccc88")
  )
}
```

The `colors` named vector is passed through `...` to the underlying trace coloring logic.

## Post-Processing Plotly Traces

For plots where the plotly object is built by a library function you can't modify, override colors by post-processing traces after `plotly_build()`:

```r
plt <- plotly::plotly_build(fig)

# Identify traces by attributes (e.g., line width, color)
for (i in seq_along(plt$x$data)) {
  trace <- plt$x$data[[i]]
  if (!is.null(trace$line$width) && trace$line$width >= 15) {
    plt$x$data[[i]]$line$color <- new_color  # override
  }
}
```

Used in: `enrichment_plot_top_enrich_gsets.R` to recolor colorbar segments (identified by `line$width >= 15`) and enrichment score line (identified by `line$color == "#00EE00"`).

### WARNING: iheatmapr → plotly post-processing does NOT work for body-only color

When `pgx.splitHeatmapFromMatrix()` produces an `Iheatmap` object, converting it to plotly via `iheatmapr::to_plotly_list() %>% plotly::as_widget()` produces multiple heatmap traces — both body cells and annotation tracks. **Do not attempt to selectively recolor body traces by inspecting `colorscale` stop positions** — the heuristic is unreliable and colors annotation tracks too.

**Correct approach:** add a `heatmap_colors` parameter to the playbase function itself and pass colors from the module:

```r
# In clustering_plot_splitmap.R — pass colors to pgx.splitHeatmapFromMatrix()
col_low  <- if (!is.null(input$color_low))  input$color_low  else get_color_theme()$secondary
col_mid  <- if (!is.null(input$color_mid))  input$color_mid  else get_color_theme()$neutral
col_high <- if (!is.null(input$color_high)) input$color_high else get_color_theme()$primary
plt <- playbase::pgx.splitHeatmapFromMatrix(
  ...,
  heatmap_colors = c(col_low, col_mid, col_high)
)
```

The playbase function defaults `heatmap_colors = NULL` (→ original brand_blue/grey/red), so existing call sites are unaffected.

## Barplot Title Color Sync

For barplot editors, link the bar color to the title color for visual consistency:

```r
effective_color <- if (!is.null(bar_color)) bar_color else "#A6CEE3"
fig <- plotly::layout(fig, title = list(font = list(color = effective_color)))
```

For multi-subplot plots with annotations as titles, use the annotation `font` property:
```r
font = list(size = 10 * title.cex, color = effective_color)
```

## Clustering Editor — Categorical Palette Selection

The `"clustering"` editor type handles plots with categorical (discrete) coloring. Unlike continuous-color editors that use `colourInput`, it provides:

1. **Palette dropdown** — `selectInput` with all `omics_pal_d` palettes + `"custom"` + `"custom_gradient"`
2. **Dynamic custom color pickers** — Server-side `renderUI` that generates one `colourInput` per group when "custom" is selected

### "default" palette normalization

The `"original"` palette option has been removed from all selectors. `"default"` is now the canonical "use this plot's natural colors" choice. Each plot maps `"default"` to its appropriate fallback in the render function:

```r
# Clustering plots (clusterpca, phenoplot, clustannot) — remap to muted_light
palette <- if (!is.null(input$palette)) input$palette else "muted_light"
if (palette %in% c("original", "default")) palette <- "muted_light"

# Stacked barplots (freq_top_gsets, correlation_barplot) — skip the override entirely
palette <- input$palette
if (!is.null(palette) && !palette %in% c("original", "default")) {
  # apply named palette
}

# Tissue / overlap plots — add "default" to the pass-through exclusion set
} else if (!is.null(palette) && !palette %in% c("original", "default", "")) {
  # apply named palette
} else {
  # original behavior (tissue group colors, Set2, etc.)
}
```

**Rule of thumb:** if the plot has a "natural" coloring from its data (tissue groups, Set2 categories, plotly default colorway), `"default"` should preserve that. If the plot uses `muted_light` as a hardcoded palette, `"default"` should remap to `muted_light`.

### "custom_gradient" palette

A special palette option added to `omics_pal_d()`. When selected, it reads 3 colors from the theme singleton (`ct$palette_c1`, `ct$palette_c2`, `ct$palette_c3`) and returns a `colorRampPalette` interpolation. All plot modules that call `omics_pal_d(palette)(n)` get gradient interpolated colors automatically — no per-module changes needed.

### Dynamic custom color pickers

Use `renderUI` to generate color pickers based on current groups. Two variants:

**Single phenotype** (e.g., clustering PCA, correlation scatter) — show group names as labels:
```r
output$custom_palette_ui <- shiny::renderUI({
  shiny::req(input$palette == "custom")
  groups <- sort(unique(as.character(pgx$samples[samples, colvar])))
  default_clrs <- rep(omics_pal_d(palette = "muted_light")(8), ceiling(length(groups) / 8))
  pickers <- lapply(seq_along(groups), function(i) {
    colourpicker::colourInput(ns(paste0("custom_color_", i)),
      label = groups[i], value = default_clrs[i])
  })
  shiny::tagList(pickers)
})
```

**Multiple phenotypes** (e.g., phenoplot) — collect group names per color position across all phenotypes:
```r
output$custom_palette_ui <- shiny::renderUI({
  shiny::req(input$palette == "custom")
  # Pre-allocate with fixed size to avoid subscript out of bounds
  level_names <- vector("list", 8)
  for (ph in phenotypes) {
    lvls <- sort(unique(as.character(Y[, ph])))
    for (j in seq_along(lvls)) {
      if (j <= 8) level_names[[j]] <- c(level_names[[j]], lvls[j])
    }
  }
  # Labels show combined group names: "Control, Female" for Color 1
  pickers <- lapply(seq_len(max_levels), function(i) {
    nms <- level_names[[i]]
    label <- if (length(nms) > 0) paste(unique(nms), collapse = ", ") else paste("Color", i)
    colourpicker::colourInput(ns(paste0("custom_color_", i)), label = label, value = default_clrs[i])
  })
  shiny::tagList(pickers)
})
```

### Collecting custom colors in the render function

```r
palette <- if (!is.null(input$palette)) input$palette else "muted_light"
if (palette %in% c("original", "default")) palette <- "muted_light"
custom_colors <- NULL
if (palette == "custom") {
  groups <- sort(unique(as.character(pgx$samples[samples, colvar])))
  custom_colors <- sapply(seq_along(groups), function(j) {
    val <- input[[paste0("custom_color_", j)]]
    if (is.null(val)) omics_pal_d(palette = "muted_light")(8)[(j - 1) %% 8 + 1] else val
  })
}
```

## Existing Implementations Reference

### Volcano (full-featured)

**File:** `components/board.expression/R/expression_plot_volcano.R`

- **Pattern:** dual-card (`cards = TRUE`), card 1 = plotly, card 2 = ggplot
- **Plotting function:** `playbase::ggVolcano()`
- **Editor panels:** Color Scheme, Text Sizes, Margins, Aspect Ratio, Labels, Prism Theme, Significance Cutoff
- **Special features:**
  - Click-to-label: clicking a point in the plotly card adds the gene to `input$label_features` in the editor (uses `parent_session` + `updateTextAreaInput`)
  - Editor `color_up`/`color_down` forwarded to **both** the static ggplot AND the dynamic plotly card via `colors = c(up=, notsig=, down=, notsel=)` argument
  - ggprism theme integration with palette, border, axis guide, legend controls
  - Hyperbolic significance cutoff option
  - Custom label controls (ggrepel box padding, segment line type, min segment length)

**Editor inputs read:**
`color_up`, `color_down`, `label_size`, `marker_size`, `axis_text_size`, `margin_checkbox`, `margin_top/right/bottom/left`, `aspect_ratio_checkbox`, `aspect_ratio`, `color_selection`, `custom_labels`, `label_features`, `box_padding`, `min_segment_length`, `label_box`, `segment_linetype`, `cutoff_type`, `hyperbola_k`, `use_ggprism`, `ggprism_palette`, `ggprism_colors`, `ggprism_border`, `ggprism_axis_guide`, `ggprism_show_legend`, `ggprism_legend_x`, `ggprism_legend_y`, `ggprism_legend_border`

**Also applied to:** `expression_plot_volcanoAll.R`, `expression_plot_volcanoMethods.R`, `enrichment_plot_volcano.R`, `enrichment_plot_volcanoall.R`, `enrichment_plot_volcanomethods.R`

### Heatmap / Splitmap (medium complexity)

**File:** `components/board.clustering/R/clustering_plot_splitmap.R`

- **Pattern:** dual-card (`cards = TRUE`), card 1 = plotly, card 2 = base
- **Plotting function:** `playbase::gx.splitmap()` (base) / `playbase::pgx.splitHeatmapFromMatrix()` (plotly)
- **Editor panels:** Labels, Color Scheme, Margins
- **Special features:**
  - Also has inline `options` (dropdown menu) for plot_type, hm_scale, show_legend, show_rownames, show_colnames — these use `ns()` not `ns_parent()`
  - Column name rotation control
  - Row names width control
  - Both plotly and base render functions read editor color inputs
  - Colors are passed via `heatmap_colors = c(col_low, col_mid, col_high)` directly to `pgx.splitHeatmapFromMatrix()` — **not** post-processed after widget conversion

**Editor inputs read:**
`label_size`, `annot_cex`, `column_names_rot`, `rownames_width`, `color_high`, `color_mid`, `color_low`, `margin_checkbox`, `margin_top/right/bottom/left`

### Barplot (simple)

**File:** `components/board.dataview/R/dataview_plot_expression.R`

- **Pattern:** single card (no `cards = TRUE`), plotly only, `plot_type = "expression_barplot"`
- **Plotting function:** built inline with `plotly::plot_ly()`
- **Editor panels:** Color Scheme, Bars Order
- **Special features:**
  - No dual-card pattern — this is a single plotly plot
  - Drag-and-drop bar reordering via `sortable::rank_list` rendered in `output$rank_list`
  - Uses `scatter_color` input (mapped to `theme$secondary`) instead of `bar_color`
  - Title color sync: links bar color to plot title color

**Editor inputs read:**
`scatter_color`, `bars_order`, `rank_list_basic`

**Also applied to:** `expression_plot_barplot.R`, `expression_plot_topfoldchange.R`, `expression_plot_topgenes.R` (multi-subplot), `enrichment_plot_barplot.R`, `enrichment_plot_geneplot.R`

### Enrichment (color override)

**File:** `components/board.enrichment/R/enrichment_plot_top_enrich_gsets.R`

- **Pattern:** single card, plotly only
- **Editor panels:** Color Scheme (up/down/line colors)
- **Special features:**
  - Post-processes plotly traces after `plotly_build()` to override colorbar segment colors and enrichment score line color
  - Identifies colorbar segments by `line$width >= 15`, enrichment line by `line$color == "#00EE00"`
  - Custom gradient built with `gplots::colorpanel(n, color_down, "#CCCCCC", color_up)`

**Editor inputs read:**
`color_up`, `color_down`, `color_line`

### Clustering / Categorical Scatterplot

**File:** `components/board.clustering/R/clustering_plot_clusterpca.R`

- **Pattern:** single card, plotly only
- **Editor panels:** Color Scheme (palette dropdown + dynamic custom color pickers)
- **Special features:**
  - `selectInput` with all `omics_pal_d` palettes + `"custom"` + `"custom_gradient"` (no `"original"`)
  - Server-side `renderUI` for dynamic `colourInput` per group when "custom" selected
  - Group labels come from the selected `colvar` phenotype
  - `"default"` palette remapped to `"muted_light"` in render

**Editor inputs read:**
`palette`, `custom_color_1..N` (dynamic)

**Also applied to:** `clustering_plot_phenoplot.R` (multi-phenotype, combined group labels), `correlation_plot_scattercorr.R` (uses `palette_default = "default"`, remaps to `"muted_light"`)

### Feature Map (continuous scatterplot)

**File:** `components/board.featuremap/R/featuremap_plot_table_gene_map.R`

- **Pattern:** dual-card, card 1 = plotly, card 2 = ggplot
- **Editor panels:** Color Scheme (low/high continuous colors), Labels (custom labels)
- **Special features:**
  - Continuous color gradient (not categorical)
  - Custom label override via `hilight2_override` parameter on `plotUMAP`

**Editor inputs read:**
`color_low`, `color_high`, `custom_labels`, `label_features`

**Also applied to:** `featuremap_plot_table_geneset_map.R`, `featuremap_plot_gene_sig.R`, `featuremap_plot_gset_sig.R`

## Conventions & Patterns

### Namespace rules

- **Editor inputs:** always use `ns_parent()` for IDs — they must be readable by the parent module's `input$`
- **Inline options** (dropdown menu): use `ns()` — they live in the plot module's own namespace
- **`conditionalPanel` with editor inputs:** pass `ns = ns_parent` so the JS condition resolves correctly:
  ```r
  conditionalPanel(
    condition = "input.margin_checkbox",
    ns = ns_parent,                           # <-- important
    numericInput(ns_parent("margin_left"), "Left", value = 10)
  )
  ```
- **For complex conditions** referencing `ns_parent` inputs by full ID:
  ```r
  condition = paste0("input['", ns_parent("cutoff_type"), "'] == 'hyperbolic'")
  ```

### Null-safe input reading

Always guard against `NULL` and `NA` — inputs may not exist yet on first render:

```r
# For color inputs — fall back to theme singleton
my_color <- if (is.null(input$my_color)) get_color_theme()$secondary else input$my_color

# For numeric inputs (can be NA if user clears the field)
my_num <- if (is.null(input$my_num) || is.na(input$my_num)) 10 else input$my_num

# For checkbox/logical inputs
my_flag <- isTRUE(input$my_flag)
```

### Dual-card pattern

Most editors use a two-card tabbed layout:

- **Card 1 ("dynamic"):** interactive plotly — editor colors are forwarded where possible (see volcano pattern)
- **Card 2 ("static"):** ggplot/base R — publication-ready, reads all editor inputs

The `PlotModuleServer` is called **twice** via `lapply` over a `plot_grid` list, once per card. The `func2` argument provides an alternate render for the fullscreen zoom modal (typically with larger fonts/sizes).

### `func` vs `func2`

- `func`: renders in the inline card (smaller)
- `func2`: renders in the fullscreen zoom modal (larger — often just the same function with bigger font sizes or different dimensions)

### Editor affects static card only (with exceptions)

For **dual-card** modules, editor inputs customize the static/ggplot card (card 2). The plotly card (card 1) remains interactive with its own built-in controls.

**Exceptions (plotly-only editors):** Many plots are single-card plotly-only. In these cases the editor affects the plotly render function directly. Examples: barplot editors, clustering scatterplot editors, enrichment editors, correlation scatter editor. The editor inputs are read inside the plotly render function, and `renderfigure_2` (the editor preview) uses the same `func()` render.

**Dual-card with forwarded colors (volcano):** Even in dual-card modules, it is desirable to forward color inputs to the plotly card so both renders stay in sync visually. See the Volcano section above for the pattern.

### Live preview in editor modal

The editor modal shows a live preview on the right side. For dual-card modules this renders `outputFunc[[2]](ns("renderfigure_2"))` — the same output as card 2 but rendered inside the editor layout. Changes to editor inputs update this preview reactively.
