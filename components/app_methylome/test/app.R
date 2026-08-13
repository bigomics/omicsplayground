##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Standalone harness for the Methylome Profiler screens, following the
## app_launcher/test convention. Renders the same four screen modules the
## bigdash shell uses, in plain bslib chrome so they can be run and reviewed
## outside the full Playground.
##
## Run with ./run.sh, or:
##   Rscript -e 'shiny::runApp("components/app_methylome/test", port=8794)'

library(shiny)
library(bslib)

source("../R/methylome_screens.R", encoding = "UTF-8")

PGX_PATH <- Sys.getenv(
  "METHYLOME_PGX",
  "/tmp/claude-1000/-home-massagno-bigomics-GitHub-omicsplayground-methyl-analysis/4d9e655f-36af-48da-8177-14e37d553518/scratchpad/meth/GSE43976-methyl-mini.pgx"
)
message("[methylome/test] loading ", PGX_PATH)
PGX <- readRDS(PGX_PATH)
message("[methylome/test] ", PGX$name, " | ", nrow(PGX$X), " x ", ncol(PGX$X))

CSS <- "
body { background:#f0f8ff; font-family:Lato,'Helvetica Neue',Arial,sans-serif; }
.mp-head { background:#fff; border-bottom:2px solid #004ca7; padding:10px 18px;
  display:flex; align-items:center; gap:14px; }
.mp-head b { font-size:16px; }
.mp-head .meta { color:#697586; font-size:12.5px; font-family:ui-monospace,Menlo,monospace; }
.mp-wrap { padding:16px 18px 40px; }
.mp-card { background:#fff; border:1.5px solid #c9ccc7; border-radius:8px;
  overflow:hidden; margin-bottom:14px; }
.mp-card-head { display:flex; align-items:center; gap:9px; padding:10px 13px;
  border-bottom:1px solid #e2e4e0; font-weight:700; font-size:13px; }
.mp-card-body { padding:12px 13px; }
.mp-sp { flex:1; }
.mp-row2 { display:grid; grid-template-columns:1fr 1fr; gap:14px; }
.mp-row4 { display:grid; grid-template-columns:repeat(4,1fr); gap:12px; margin-bottom:14px; }
.mp-tile { border:1px solid #e2e4e0; border-radius:7px; background:#fff; padding:10px 12px; }
.mp-k { color:#697586; font:700 9.5px/1.4 ui-monospace,Menlo,monospace;
  letter-spacing:.08em; text-transform:uppercase; }
.mp-v { margin-top:4px; font-size:21px; font-weight:700; font-variant-numeric:tabular-nums; }
.mp-d { margin-top:2px; color:#697586; font-size:11px; }
.nav-tabs .nav-link.active { font-weight:700; }
table.dataTable { font-size:12px; }
"

ui <- bslib::page_fluid(
  padding = 0,
  tags$head(tags$style(HTML(CSS))),
  div(class = "mp-head",
    tags$b("Methylome Profiler"),
    span(class = "meta", sprintf("%s  ·  %s probes × %s samples  ·  %s",
      PGX$name, format(nrow(PGX$X), big.mark = ","), ncol(PGX$X), PGX$datatype))
  ),
  div(class = "mp-wrap",
    tabsetPanel(
      id = "tabs",
      tabPanel("Sample ledger", value = "ledger",
               div(style = "padding-top:14px", methylome_ledger_ui("ledger"))),
      tabPanel("Epigenetic age", value = "age",
               div(style = "padding-top:14px", methylome_age_ui("age"))),
      tabPanel("Methylome character", value = "character",
               div(style = "padding-top:14px", methylome_character_ui("character"))),
      tabPanel("EWAS", value = "ewas",
               div(style = "padding-top:14px", methylome_ewas_ui("ewas")))
    )
  )
)

server <- function(input, output, session) {
  P <- shiny::reactive(PGX)

  ## Deep-linkable screens: ?tab=ledger|age|character|ewas
  shiny::observe({
    q <- shiny::parseQueryString(session$clientData$url_search)
    if (!is.null(q$tab)) shiny::updateTabsetPanel(session, "tabs", selected = q$tab)
  })

  methylome_ledger_server("ledger", P)
  methylome_age_server("age", P)
  methylome_character_server("character", P)
  methylome_ewas_server("ewas", P)
}

shinyApp(ui, server)
