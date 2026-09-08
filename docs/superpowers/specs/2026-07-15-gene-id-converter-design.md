# Gene ID Converter applet (`app_convert`) — Design

## Purpose

A standalone Shiny module, structured the same way as `app_prism`, that lets a
user pick an organism, paste a list of gene/feature IDs, and get back the
full annotation table (symbol, gene title, biotype, chromosome, etc.) from
`playbase::getGeneAnnotation()`.

Supported organisms for now: Human, Mouse, Rat, Fall armyworm.

## Structure

New folder `components/app_convert/`, mirroring `app_prism`:

```
app_convert/
  R/
    convert_ui.R      # convert_ui(id)
    convert_server.R  # convert_server(id)
```

No other files (no static assets, no extra R files) — same minimal footprint
as `app_prism`.

## UI (`convert_ui`)

- Title bar: "Gene ID Converter".
- `bslib::layout_columns` with a sidebar + main results area (same visual
  weight/pattern as Prism's control column + output column).
- Sidebar controls:
  - `shiny::selectInput(ns("organism"), "Organism:", choices = c("Human",
    "Mouse", "Rat", "Fall armyworm"), selected = "Human")`
  - `shiny::textAreaInput(ns("features"), "Gene/feature IDs (one per
    line):", rows = 12, placeholder = "e.g.\nTP53\nENSG00000141510\n...")`
  - `shiny::actionButton(ns("convert"), "Convert", icon = icon("arrows-rotate"))`
  - `shiny::downloadButton(ns("download"), "Download CSV")`
- Main panel: `DT::dataTableOutput(ns("table"))` showing the full annotation
  result.

## Server (`convert_server`)

- Organism label → scientific name passed to `getGeneAnnotation()`:
  - Human → `Homo sapiens`
  - Mouse → `Mus musculus`
  - Rat → `Rattus norvegicus`
  - Fall armyworm → `Spodoptera frugiperda`

  (Human/Mouse/Rat resolve via `playbase::SPECIES_TABLE`/AnnotationHub;
  Fall armyworm is not in that table, so it only resolves through the
  orthogene/g:profiler fallback inside `getGeneAnnotation()` — this is why
  `methods` is left at its default `c("annothub", "gprofiler")` rather than
  restricted to AnnotationHub only.)

- `eventReactive(input$convert, ...)`:
  - Parse `input$features`: split on newlines/commas/whitespace, trim
    whitespace, drop empty strings, no de-duplication (preserve user's
    input order/count so output rows line up 1:1 with what they pasted).
  - `shiny::validate(shiny::need(length(probes) > 0, "Please paste at least
    one gene/feature ID."))`
  - Call `playbase::getGeneAnnotation(organism = <mapped organism>, probes =
    probes, verbose = FALSE)`.
  - `shiny::validate(shiny::need(!is.null(result), "Could not annotate
    these IDs for the selected organism."))`
  - Return the resulting data frame.

- `output$table`: renders the full result data frame via `DT::renderDataTable`.
- `output$download`: `shiny::downloadHandler()` writing the same data frame
  to CSV (`write.csv(..., row.names = FALSE)`).

## Wiring into the main app

Same pattern as Prism, all DEVMODE-gated:

1. `00SourceAll.R` — regenerate via `dev/create_source_all.R` (it's a
   generated file, not hand-edited) so it picks up
   `app_convert/R/convert_server.R` and `app_convert/R/convert_ui.R`.
2. `app/R/ui.R` — add a hidden nav panel, same shape as Prism's:
   ```r
   bslib::nav_panel_hidden("IDconvert",
     div(convert_ui("idconvert"), class = 'px-4 py-0')
   )
   ```
3. `app/R/server.R` — inside the existing `if (opt$DEVMODE) { ... }` block
   (where `prism_server("prism")` is called), add:
   ```r
   convert_server("idconvert")
   ```

No changes needed to `app_tools`: `tools_ui.R` already has an "ID
SmartConverter" card and `tools_server.R` already wires its button to
`bslib::nav_select("app-sidebar", "IDconvert", session = parent)`. The tool
becomes reachable as soon as the "IDconvert" panel is registered.

## Out of scope

- File upload for the ID list (paste-only for now).
- Editing/expanding the `app_tools` card (already exists and already wired).
- Any organism beyond the four listed.
