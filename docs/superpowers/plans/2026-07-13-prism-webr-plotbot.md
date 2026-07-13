# PRISM webR Chartbot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the PRISM board's server-side `eval(parse(text=plotcode))` of LLM-generated ggplot2 code with client-side, sandboxed execution in webR (R-in-WASM), porting the execution model from `~/Playground/testx/plotbot-claude`.

**Architecture:** The LLM call and prompt-building (`omicsai::omicsai_call_llm`) stay server-side and unchanged in spirit. The generated R code, plus a CSV snapshot of the current reactive dataframe, is pushed to the browser via `session$sendCustomMessage`. A new JS module (`prism-webr.js`) runs webR, evaluates the code in-browser, captures an SVG, and injects it into the plot panel. A chat-bubble message log is added above the existing input row. Two triggers feed the same render pipeline: chat-send (calls the LLM) and sidebar changes (mechanically patches the last code, no LLM call) — both end at the same `session$sendCustomMessage("prism-executeCode", ...)` call.

**Tech Stack:** R / Shiny (module `prism_server`/`prism_ui`), `omicsai` (existing LLM client), webR (`https://webr.r-wasm.org/latest/webr.mjs`, loaded as an ES module), testthat (existing repo test framework).

## Global Constraints

- No change to the LLM backend: keep using `omicsai::omicsai_call_llm()` (per spec, "Non-goals").
- Remove `"xkcd"` from the theme dropdown — package not installable under webR/WASM (per spec decision).
- Preloaded webR packages, fixed set: `ggplot2, ggrepel, ggprism, ggpubr, ggsignif` (per spec).
- Custom message types must be prefixed `prism-` to avoid collisions with other custom message handlers already registered elsewhere in the app (verified: no existing handler named `executeCode` or `addMessage`, but `prism-` prefix matches this codebase's existing convention, e.g. `copilot-chat-batch-append`, `bigdash-*`).
- The `prism` Shiny module is mounted exactly once, with id `"prism"` (`components/app/R/ui.R:189`, `components/app/R/server.R:1053`). The new JS file assumes `ns()`-generated ids resolve to the `prism-` prefix; this coupling is documented with a one-line comment in the JS file, not solved generically (per spec: "cheap to do right," not "build for reentrancy").
- New static asset goes in `components/app/R/www/` (the app's existing single shared static folder, served at `static/` via `addResourcePath("static", ...)` in `components/app/R/global.R`), following the flat-file convention already used there (`bigomics-extra.js`, `temp.js`, etc.).
- Follow existing test convention: pure helper functions are defined at file scope (outside `moduleServer`) in `prism_server.R` so `tests/testthat/test-prism-server.R` can `source()` the file directly and call them (same pattern as `tests/testthat/test-copilot-reports-server.R` sourcing `CopilotReportsServer.R`).

---

### Task 1: Add the webR client-side execution script

**Files:**
- Create: `components/app/R/www/prism-webr.js`

**Interfaces:**
- Produces: two Shiny custom message handlers the server will call in Task 3 —
  `prism-executeCode` (payload `{code: string, csv: string}`) and `prism-addMessage`
  (payload `{role: "user"|"assistant"|"error", content: string}`).
- Produces: DOM element ids it expects to find (added to `prism_ui.R` in Task 4):
  `prism-webr-status`, `prism-plot-placeholder`, `prism-plot-container`,
  `prism-plot-error`, `prism-chat-messages`, `prism-chartbot_send`.

- [ ] **Step 1: Write the file**

```javascript
// components/app/R/www/prism-webr.js
//
// Client-side, sandboxed R execution for the PRISM board's chartbot.
// Runs LLM-generated ggplot2 code inside webR (R compiled to WASM) instead of
// evaluating it on the server, so arbitrary generated code never runs server-side.
//
// Assumes the "prism" Shiny module is mounted with module id "prism" (see
// components/app/R/ui.R), so every ns()-generated element id used here is the
// literal "prism-" prefix below.
import { WebR } from 'https://webr.r-wasm.org/latest/webr.mjs';

const ID_PREFIX = 'prism-';
const PACKAGES = ['ggplot2', 'ggrepel', 'ggprism', 'ggpubr', 'ggsignif'];

let webR = null;
let isReady = false;

function el(name) {
  return document.getElementById(ID_PREFIX + name);
}

function status(msg, color = '#64748b') {
  const node = el('webr-status');
  if (node) { node.textContent = msg; node.style.color = color; }
}

function setSendEnabled(enabled) {
  const btn = el('chartbot_send');
  if (btn) btn.disabled = !enabled;
}

function showPlot(svgContent) {
  const placeholder = el('plot-placeholder');
  const container = el('plot-container');
  const errorEl = el('plot-error');

  if (placeholder) placeholder.style.display = 'none';
  if (errorEl) errorEl.style.display = 'none';

  container.innerHTML = svgContent;

  const svg = container.querySelector('svg');
  if (svg) {
    svg.style.maxWidth = '100%';
    svg.style.height = 'auto';
    svg.removeAttribute('width');
    svg.removeAttribute('height');
  }
}

function showError(msg) {
  const node = el('plot-error');
  if (node) {
    node.textContent = 'Plot error: ' + msg;
    node.style.display = 'block';
  }
}

// Escape a string for safe embedding inside an R double-quoted string literal.
function escapeRString(s) {
  return s
    .replace(/\\/g, '\\\\')
    .replace(/"/g, '\\"')
    .replace(/\r\n/g, '\\n')
    .replace(/\n/g, '\\n')
    .replace(/\r/g, '\\n');
}

async function initWebR() {
  try {
    setSendEnabled(false);
    status('Initialising webR runtime…');

    webR = new WebR();
    await webR.init();

    status('Installing webR packages…');
    await webR.installPackages(PACKAGES, { quiet: true });
    await webR.evalR(
      'suppressPackageStartupMessages({' +
      PACKAGES.map((p) => `library(${p})`).join(';') +
      '})'
    );

    isReady = true;
    setSendEnabled(true);
    status('webR ready — ask for a plot!', '#059669');
  } catch (err) {
    status('webR error: ' + err.message, '#dc2626');
    console.error('[prism webR] init error:', err);
  }
}

async function executePlotCode(code, csv) {
  if (!isReady) {
    console.warn('[prism webR] not ready yet, ignoring executeCode');
    return;
  }

  status('Rendering plot…', '#64748b');

  const safeCode = code.replace(/\\/g, '\\\\');
  const safeCsv = escapeRString(csv || '');

  const wrapped = `
local({
  data <- read.csv(text = "${safeCsv}", row.names = 1)
  tmp <- tempfile(fileext = ".svg")
  grDevices::svg(tmp, width = 8, height = 5.5, onefile = FALSE)
  on.exit({ if (dev.cur() > 1) dev.off(); unlink(tmp) }, add = TRUE)
  tryCatch(
    {
      .plotbot_result <- withVisible({ ${safeCode} })
      if (.plotbot_result$visible) print(.plotbot_result$value)
    },
    error = function(e) stop(conditionMessage(e))
  )
  dev.off()
  paste(readLines(tmp, warn = FALSE), collapse = "\\n")
})`;

  try {
    const result = await webR.evalR(wrapped);
    const js = await result.toJs();
    const svg = js.values[0];

    showPlot(svg);
    status('Plot rendered!', '#059669');
  } catch (err) {
    showError(err.message);
    addChatMessage({ role: 'error', content: 'Plot error: ' + err.message });
    status('Plot error — see panel below.', '#d97706');
    console.error('[prism webR] plot error:', err);
  }
}

function addChatMessage({ role, content }) {
  const container = el('chat-messages');
  if (!container) return;

  const div = document.createElement('div');
  div.className = `prism-msg prism-msg-${role}`;
  div.textContent = content;

  container.appendChild(div);
  container.scrollTop = container.scrollHeight;
}

function bootstrap() {
  if (!window.Shiny || !window.Shiny.addCustomMessageHandler) {
    setTimeout(bootstrap, 150);
    return;
  }

  Shiny.addCustomMessageHandler('prism-executeCode', ({ code, csv }) => executePlotCode(code, csv));
  Shiny.addCustomMessageHandler('prism-addMessage', (data) => addChatMessage(data));

  initWebR();
}

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', bootstrap);
} else {
  bootstrap();
}
```

- [ ] **Step 2: Validate JS syntax**

Run: `node --check components/app/R/www/prism-webr.js`
Expected: no output, exit code 0.

- [ ] **Step 3: Commit**

```bash
git add components/app/R/www/prism-webr.js
git commit -m "feat(prism): add webR client-side plot execution script"
```

---

### Task 2: Extract and test pure plot-code helper functions

**Files:**
- Modify: `components/app_prism/R/prism_server.R` (add four functions at file scope, above the existing `prism_server <- function(id) {...}` definition; do not touch the `moduleServer` body yet — that's Task 3)
- Test: `tests/testthat/test-prism-server.R`

**Interfaces:**
- Produces: `PRISM_WEBR_PACKAGES` (character scalar) — `"ggplot2, ggrepel, ggprism, ggpubr, ggsignif"`.
- Produces: `df_to_csv(data)` — takes a data.frame, returns a character scalar of CSV text (header + rownames as first column, via `write.csv(..., row.names = TRUE)`).
- Produces: `strip_code_fences(text)` — takes a character scalar, returns it with any ```` ```r ````, ```` ```R ```` or ```` ``` ```` fence markers removed.
- Produces: `build_prism_prompt(msg, vars, rows, last_plotcode, dataset, pointsize, fontsize, theme, packages)` — returns a single character scalar prompt string.
- Produces: `update_plotcode(last_plotcode, pointsize, fontsize, theme, dataset)` — same behavior as today's closure-based version, but `last_plotcode` is now an explicit first argument instead of a captured variable, so it is a pure function. Returns `NULL` when `last_plotcode == ""`, otherwise a character scalar.
- Consumes: nothing from other tasks.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-prism-server.R`:

```r
## test-prism-server.R

.board_dir <- if (dir.exists("components/app_prism/R")) {
  "components/app_prism/R"
} else {
  "../../components/app_prism/R"
}

source(file.path(.board_dir, "prism_server.R"), local = TRUE)

test_that("df_to_csv round-trips a data.frame including row names", {
  df <- data.frame(x = 1:2, y = c("a", "b"), row.names = c("r1", "r2"))
  csv <- df_to_csv(df)

  expect_type(csv, "character")
  expect_length(csv, 1)

  back <- read.csv(text = csv, row.names = 1)
  expect_equal(back$x, 1:2)
  expect_equal(as.character(back$y), c("a", "b"))
  expect_equal(rownames(back), c("r1", "r2"))
})

test_that("strip_code_fences removes r/R/plain fence markers", {
  expect_equal(strip_code_fences("```r\nggplot(mtcars)\n```"), "\nggplot(mtcars)\n")
  expect_equal(strip_code_fences("```R\nggplot(mtcars)\n```"), "\nggplot(mtcars)\n")
  expect_equal(strip_code_fences("```\nggplot(mtcars)\n```"), "\nggplot(mtcars)\n")
  expect_equal(strip_code_fences("no fences here"), "no fences here")
})

test_that("build_prism_prompt embeds the request, data context and preloaded packages", {
  p <- build_prism_prompt(
    msg = "plot mpg vs wt",
    vars = "mpg, wt",
    rows = "Mazda RX4, Datsun 710",
    last_plotcode = "",
    dataset = "mtcars",
    pointsize = 3,
    fontsize = 12,
    theme = "classic",
    packages = PRISM_WEBR_PACKAGES
  )

  expect_type(p, "character")
  expect_match(p, "plot mpg vs wt", fixed = TRUE)
  expect_match(p, "mpg, wt", fixed = TRUE)
  expect_match(p, "Mazda RX4, Datsun 710", fixed = TRUE)
  expect_match(p, PRISM_WEBR_PACKAGES, fixed = TRUE)
  expect_match(p, "never call library\\(\\)")
})

test_that("update_plotcode returns NULL when there is no last plotcode", {
  expect_null(update_plotcode("", 3, 12, "classic", "mtcars"))
})

test_that("update_plotcode patches theme, point size, font size and title onto the last code", {
  last <- "ggplot(data, aes(x=wt,y=mpg)) + geom_point()"
  out <- update_plotcode(last, 5, 20, "dark", "mtcars")

  expect_match(out, "theme_dark()", fixed = TRUE)
  expect_match(out, "geom_point(size=5)", fixed = TRUE)
  expect_match(out, "element_text(size=20)", fixed = TRUE)
  expect_match(out, "labs(title='mtcars')", fixed = TRUE)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd tests/testthat && Rscript -e 'library(testthat); source("test-prism-server.R")'`
Expected: FAIL — errors such as `could not find function "df_to_csv"`, `"strip_code_fences"`, `"build_prism_prompt"`, or `unused argument` on `update_plotcode("", ...)` (since the current definition takes 4 args, not 5).

- [ ] **Step 3: Add the four pure functions to `prism_server.R`**

In `components/app_prism/R/prism_server.R`, insert the following **above** the `prism_server <- function(id) {` line (after the existing `if(0) {...}` font-install block):

```r
PRISM_WEBR_PACKAGES <- "ggplot2, ggrepel, ggprism, ggpubr, ggsignif"

df_to_csv <- function(data) {
  paste(utils::capture.output(utils::write.csv(data, row.names = TRUE)), collapse = "\n")
}

strip_code_fences <- function(text) {
  gsub("```[rR]|```", "", text)
}

build_prism_prompt <- function(msg, vars, rows, last_plotcode, dataset,
                                pointsize, fontsize, theme, packages) {
  paste(
    "You are asked the following request about modifying a ggplot figure.",
    "Just give the raw plotting code. No explanations.",
    "The code runs in a sandboxed R environment where only these packages",
    "are pre-loaded: ", packages, "-- never call library() and never use any other package.",
    "\nThis is the request of the user: ", msg,
    "\nThese are the variables in the dataframe called 'data': ", vars,
    "\nThese are the rownames of 'data': ", rows,
    "\nThis is the last plotting code of the graph: ", last_plotcode,
    "\nSet default title as dataset name: ", dataset,
    "\nDefault point size unless asked by user: ", pointsize,
    "\nDefault font size unless asked by user: ", fontsize,
    "\nDefault theme unless asked by user: ", theme
  )
}

update_plotcode <- function(last_plotcode, pointsize, fontsize, theme, dataset) {
  if (last_plotcode == "") return(NULL)
  plotcode <- trimws(last_plotcode)
  plotcode <- gsub("\\)$", ") +\n", plotcode)
  plotcode <- paste0(plotcode, "theme_", theme, "()")
  plotcode <- paste0(plotcode, " +\n geom_point(size=", pointsize, ")")
  plotcode <- paste0(plotcode, " +\n theme(text=element_text(size=", fontsize, "))")
  plotcode <- paste0(plotcode, " +\n labs(title='", dataset, "')")
  return(plotcode)
}
```

Do **not** delete the old `update_plotcode` definition or the `msg2 <- paste(...)` block that still live inside `prism_server()`'s `moduleServer` body yet — leave them in place for now. (There will briefly be two `update_plotcode` definitions in the file; the file-scope one defined here will be the one exported and used by the tests. Task 3 removes the old in-module one and rewires everything.)

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd tests/testthat && Rscript -e 'library(testthat); source("test-prism-server.R")'`
Expected: `Test passed` for all 5 `test_that` blocks, no failures.

- [ ] **Step 5: Commit**

```bash
git add components/app_prism/R/prism_server.R tests/testthat/test-prism-server.R
git commit -m "feat(prism): add pure, tested helpers for plot-code prompt/CSV/patch logic"
```

---

### Task 3: Wire `prism_server.R` to drive webR instead of `renderPlot`

**Files:**
- Modify: `components/app_prism/R/prism_server.R`

**Interfaces:**
- Consumes: `PRISM_WEBR_PACKAGES`, `df_to_csv()`, `strip_code_fences()`, `build_prism_prompt()`, `update_plotcode(last_plotcode, pointsize, fontsize, theme, dataset)` from Task 2.
- Consumes: custom message types `prism-executeCode` / `prism-addMessage` and their payload shapes from Task 1's `prism-webr.js`.
- Produces: no more `output$plot1`; the module now drives the plot via `session$sendCustomMessage("prism-executeCode", ...)`.

- [ ] **Step 1: Replace the `prism_server` function body**

Replace the entire `prism_server <- function(id) { moduleServer(id, function(input, output, session) { ... }) }` block (everything from `prism_server <- function(id) {` to the final closing of that function) with:

```r
prism_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    get_dataframe <- reactive({
      if(input$dataset == "mtcars") {
        df <- within(mtcars, {
          vs <- factor(vs)
          am <- factor(am)
          cyl  <- factor(cyl)
          gear <- factor(gear)
        })
      }
      if(input$dataset == "iris") {
        df <- iris
      }
      if(input$dataset == "geiger") {
        pgx <- playdata::GEIGER_PGX
        mm <- playbase::pgx.getMetaMatrix(pgx)
        df <- data.frame(logFC = mm$fc[,1], pv = mm$pv[,1])
      }
      df 
    })
        
    last_plotcode <- ""

    llm_plotcode <- eventReactive( input$chartbot_send, {        

      dbg("chartbot_user_input reacted!") 

      msg <- input$chartbot_user_input
      if(msg == "") return(NULL)

      session$sendCustomMessage("prism-addMessage", list(role = "user", content = msg))

      if(msg == "reset") {
        last_plotcode <<- ""
        shiny::updateTextInput(session, "chartbot_user_input", value = "",
          placeholder = "What do you want to plot?")
        session$sendCustomMessage("prism-addMessage",
          list(role = "assistant", content = "Cleared the plot."))
        return("ggplot() + theme_void()")
      }
      
      data <- get_dataframe()
      vars <- paste(colnames(data),collapse=", ")
      rows <- paste(rownames(data),collapse=", ")      

      pointsize <- input$pointsize
      fontsize <- input$fontsize
      theme <- input$theme
      dataset <- input$dataset
      
      msg2 <- build_prism_prompt(msg, vars, rows, last_plotcode, dataset,
        pointsize, fontsize, theme, PRISM_WEBR_PACKAGES)
      
      ai_model <- getUserOption(session, "llm_model")

      plotcode <- tryCatch({
        raw <- omicsai::omicsai_call_llm(
          msg2,
          omicsai::omicsai_config(model = ai_model)
        )$text
        raw <- paste0(raw, "\n")
        strip_code_fences(raw)
      }, error = function(e) {
        session$sendCustomMessage("prism-addMessage",
          list(role = "error", content = paste("LLM error:", conditionMessage(e))))
        NULL
      })

      shiny::req(plotcode)

      last_plotcode <<- plotcode      
      shiny::updateTextInput(session, "chartbot_user_input", value = "",
        placeholder = "Any edits?")
      session$sendCustomMessage("prism-addMessage",
        list(role = "assistant", content = "Updated the plot."))

      return(plotcode)      
    })

    get_plotcode <- reactive({
      pointsize <- input$pointsize
      fontsize <- input$fontsize
      theme <- input$theme
      dataset <- input$dataset

      dbg("input$chartbot_send = ", input$chartbot_send)
      dbg("input$chartbot_user_input = ", isolate(input$chartbot_user_input))
      
      if(isolate(input$chartbot_user_input) == "") {
        dbg("updating plotcode")
        code <- update_plotcode(last_plotcode, pointsize, fontsize, theme, dataset)
      } else {
        dbg("retrieving plotcode")
        code <- llm_plotcode()
      }
      code
    })

    observe({
      plotcode <- get_plotcode()
      shiny::req(plotcode)
      data <- get_dataframe()
      csv <- df_to_csv(data)
      session$sendCustomMessage("prism-executeCode", list(code = plotcode, csv = csv))
    })

    output$plotcode <- renderUI({
      plotcode <- get_plotcode()
      shiny::req(plotcode)
      plotcode <- sub("ggplot\\(","ggplot(<br>&nbsp;&nbsp;",plotcode)
      plotcode <- gsub("\n","<br>",plotcode)
      HTML(plotcode)      
    })

    output$data1 <- renderDataTable({
      data <- get_dataframe()
      return(data)
    })
    
  })
}
```

This removes the old in-module `update_plotcode()` definition (now only the file-scope one from Task 2 exists), removes `output$plot1 <- renderPlot({...})` and its `require(...)` conditional loads, and removes the old inline `msg2 <- paste(...)` block in favor of `build_prism_prompt()`.

- [ ] **Step 2: Re-run the Task 2 tests to confirm nothing broke**

Run: `cd tests/testthat && Rscript -e 'library(testthat); source("test-prism-server.R")'`
Expected: `Test passed` for all 5 `test_that` blocks (these test the file-scope helpers, which are unaffected by the `moduleServer` body rewrite, but sourcing must still succeed with no duplicate-definition surprises).

- [ ] **Step 3: Parse-check the whole file**

Run: `Rscript -e 'parse("components/app_prism/R/prism_server.R")' `
Expected: no error (confirms valid R syntax after the rewrite).

- [ ] **Step 4: Commit**

```bash
git add components/app_prism/R/prism_server.R
git commit -m "feat(prism): drive plot rendering via webR custom messages instead of server eval()"
```

---

### Task 4: Update `prism_ui.R` — webR script, chat log, drop xkcd theme

**Files:**
- Modify: `components/app_prism/R/prism_ui.R`

**Interfaces:**
- Consumes: `prism-webr.js` from Task 1 (loaded via `<script type="module">`), and the DOM ids it expects (`prism-webr-status`, `prism-plot-placeholder`, `prism-plot-container`, `prism-plot-error`, `prism-chat-messages`, `prism-chartbot_send`) — all produced here via `ns(...)`, which resolves to those exact strings because the module is mounted with id `"prism"`.

- [ ] **Step 1: Remove `"xkcd"` from the theme choices**

In `components/app_prism/R/prism_ui.R`, change:

```r
        shiny::selectInput(ns("theme"), "Theme:", sort(c("gray","bw","light","dark",
          "minimal","classic","xkcd","prism")), selected="gray"),
```

to:

```r
        shiny::selectInput(ns("theme"), "Theme:", sort(c("gray","bw","light","dark",
          "minimal","classic","prism")), selected="gray"),
```

- [ ] **Step 2: Add the webR script tag and CSS to the page head**

Change:

```r
  ui <- bslib::page_fluid(
    div(class = "navbar navbar-static-top", div(title, class = "container-fluid"),
      style="margin-top: 24px;"),
```

to:

```r
  ui <- bslib::page_fluid(
    tags$head(
      tags$script(type = "module", src = "static/prism-webr.js"),
      tags$style(HTML("
        .prism-chat-messages {
          max-height: 160px; overflow-y: auto; padding: 8px;
          display: flex; flex-direction: column; gap: 6px;
          border: 1px solid #e2e8f0; border-radius: 6px; margin-bottom: 8px;
        }
        .prism-msg {
          padding: 8px 12px; border-radius: 10px; max-width: 90%;
          font-size: 0.85rem; line-height: 1.4;
          white-space: pre-wrap; word-break: break-word;
        }
        .prism-msg-user { background:#6366f1; color:#fff; align-self:flex-end; }
        .prism-msg-assistant { background:#f1f5f9; color:#334155; align-self:flex-start; border:1px solid #e2e8f0; }
        .prism-msg-error { background:#fee2e2; color:#991b1b; align-self:flex-start; border:1px solid #fca5a5; }
        .prism-plot-error {
          color:#dc2626; background:#fee2e2; padding:1rem; border-radius:8px;
          font-family:monospace; font-size:0.82rem; display:none; margin-top:8px;
        }
        .prism-webr-status {
          font-size:0.8rem; padding:5px 10px; border-radius:5px;
          background:#f8fafc; border:1px solid #e2e8f0; color:#64748b; margin-bottom:8px;
        }
      "))
    ),
    div(class = "navbar navbar-static-top", div(title, class = "container-fluid"),
      style="margin-top: 24px;"),
```

- [ ] **Step 3: Replace the plot panel markup**

Change:

```r
          bslib::nav_panel(
            title = "plot",
            shiny::plotOutput(ns("plot1"), height='600px'),
            ## shinychat::chat_ui(
            ##   ns("chartbot"),
            ##   style = "max-height: 180px; width: min(800px, 100%);",
            ##   fill = FALSE
            ## )
            div(
              shiny::textInput(ns("chartbot_user_input"),"",
                placeholder = "What do you want to plot?", width=600),
              shiny::actionButton(ns("chartbot_send"),"send",
                icon = icon("arrow-right-from-bracket"))
            )
          ),
```

to:

```r
          bslib::nav_panel(
            title = "plot",
            div(id = ns("webr-status"), class = "prism-webr-status",
              "Initializing webR runtime…"),
            div(id = ns("plot-placeholder"),
              style = "color:#94a3b8; text-align:center; padding:2rem;",
              tags$div(style = "font-size:2.5rem;", "\U0001f4c8"),
              tags$div("Your plot will appear here")
            ),
            div(id = ns("plot-container")),
            div(id = ns("plot-error"), class = "prism-plot-error"),
            div(id = ns("chat-messages"), class = "prism-chat-messages"),
            div(
              shiny::textInput(ns("chartbot_user_input"),"",
                placeholder = "What do you want to plot?", width=600),
              shiny::actionButton(ns("chartbot_send"),"send",
                icon = icon("arrow-right-from-bracket"))
            )
          ),
```

- [ ] **Step 4: Parse-check the file**

Run: `Rscript -e 'parse("components/app_prism/R/prism_ui.R")'`
Expected: no error.

- [ ] **Step 5: Commit**

```bash
git add components/app_prism/R/prism_ui.R
git commit -m "feat(prism): render plot/status/chat via webR-driven DOM instead of plotOutput"
```

---

### Task 5: End-to-end manual verification

**Files:** none (verification only).

- [ ] **Step 1: Launch the app**

Use the `run` skill (or the project's normal dev launcher, e.g. `dev/board.launch`) to start the full Omics Playground app, then open it in a browser and navigate to the PRISM tab.

- [ ] **Step 2: Confirm webR boot sequence**

Expected: the status line shows "Initialising webR runtime…" then "Installing webR packages…" then "webR ready — ask for a plot!" (~30s), and the send button is disabled until ready, matching `components/app/R/www/prism-webr.js`'s `setSendEnabled`/`status` calls.

- [ ] **Step 3: Chat-driven render, all three datasets**

For each of `mtcars`, `iris`, `geiger` in the Dataset selector: type a plotting request (e.g. "Plot mpg vs wt, color by cylinder" for mtcars) and click send.
Expected: a user chat bubble appears, then an assistant "Updated the plot." bubble, the plot panel shows an SVG, and the code wellPanel in the sidebar shows the generated code.

- [ ] **Step 4: Sidebar-driven mechanical update**

After a successful chat-driven plot, change Point size, Font size, or Theme in the sidebar.
Expected: the plot re-renders immediately (patched via `update_plotcode()`, sent to webR) with no new chat bubble and no LLM call.

- [ ] **Step 5: Error path**

Send a prompt likely to produce bad/unrunnable code (e.g. "call a function that does not exist, name it foobarbaz()").
Expected: the red plot-error panel shows the webR error message, and/or an error-role chat bubble appears.

- [ ] **Step 6: Reset**

Type `reset` and send.
Expected: the plot clears to an empty `theme_void()` plot, and an assistant "Cleared the plot." bubble appears.

- [ ] **Step 7: Confirm xkcd theme is gone**

Expected: the Theme dropdown no longer offers "xkcd".

No commit for this task (verification only, no code changes).
