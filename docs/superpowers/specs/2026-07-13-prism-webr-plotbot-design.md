# PRISM board: client-side webR execution for the chartbot

## Goal

Port the client-side, sandboxed R execution model from `~/Playground/testx/plotbot-claude`
(a standalone "AI PlotBot" Shiny app) into the `app_prism` module of Omics Playground.

The PRISM board already has an experimental "chartbot": the user types a plotting
request, `omicsai::omicsai_call_llm()` generates ggplot2 code from a prompt built with
the current dataframe's structure and the last plot code, and that code is currently
run with `eval(parse(text=plotcode))` **on the server** inside `renderPlot`.

`plotbot-claude` solves the same problem differently: it ships the LLM-generated code
to the browser, where **webR** (R compiled to WASM) executes it in a sandboxed runtime
and renders an SVG. No generated code ever runs on the server.

This is the piece worth porting: it removes a real server-side arbitrary-code-execution
risk (and server compute cost) from the PRISM board, while keeping the existing
LLM-calling logic (`omicsai`) and prompt-construction as-is.

## Non-goals

- No change to the LLM backend (`omicsai::omicsai_call_llm` stays).
- No change to prompt-construction logic beyond declaring the fixed set of R packages
  available in the browser sandbox (see below).
- Not porting `plotbot-claude`'s `ellmer`/OpenAI-direct calling — PRISM already uses
  `omicsai` and that stays.

## Architecture / data flow

Today:

```
chat text -> build prompt -> omicsai_call_llm() -> strip fences -> eval(parse(...)) on server -> renderPlot
```

Target:

```
chat text -> build prompt (unchanged) -> omicsai_call_llm() -> strip fences
  -> ship {code, csv-of-dataframe} to browser
  -> webR (WASM R) evals code in sandbox, renders SVG
  -> SVG injected into plot panel; code panel + chat log updated
```

Two triggers feed the same client-side render pipeline, preserving current behavior:

- **Chat send** -> `llm_plotcode()` calls the LLM, returns fresh code, becomes the new
  "last plotcode."
- **Sidebar change** (dataset/theme/pointsize/fontsize) -> `update_plotcode()`
  mechanically patches the *existing* last plotcode (no LLM call), exactly as today.

Either path ends by pushing `{code, csv}` to the browser via
`session$sendCustomMessage`, which runs webR and updates: (a) the plot panel (SVG),
(b) the existing code wellPanel (still server-rendered from the reactive, unchanged),
and (c) a chat bubble log (new).

### Data transfer into the sandbox

webR's WASM R has no access to server-side R objects. The current reactive `data`
data.frame (mtcars, iris, or the computed geiger fold-change table) is serialized to
CSV via `write.csv` and sent alongside the code in the same custom message, every time
(no caching). JS embeds it as `data <- read.csv(text = "...")` ahead of the plot code
before evaluating. These datasets are small (built-in R datasets or a two-column
fold-change table), so re-sending on every render is cheap and avoids state-tracking
complexity.

## Component changes

### `components/app/R/www/prism-webr.js` (new file)

Follows the existing flat-file convention of the shared static `www/` folder (served
at `static/`, per `addResourcePath("static", ...)` in `components/app/R/global.R`).
Adapted from `plotbot-claude/www/webr-plot.js`:

- `initWebR()`: installs `ggplot2, ggrepel, ggprism, ggpubr, ggsignif`, attaches
  libraries via `library()` in the webR session, flips a ready flag, and disables the
  send button / shows an "initializing…" status until ready (~30s cold start on first
  load — unavoidable, same cost `plotbot-claude` pays). This fixes a small UX gap in
  the source app, which silently no-ops on early sends.
- `executePlotCode(code, csv)`: gains a second argument versus the source. Wraps
  `data <- read.csv(text = "...")` ahead of the plot code, evals in webR inside the
  existing SVG-capture wrapper, injects the resulting SVG, and updates the code panel.
- `addChatMessage({role, content})`: unchanged bubble-append logic (user/assistant/error
  roles), reused as-is.
- DOM ids used by this script are namespaced (matching the module's `ns()` prefix)
  even though `app_prism` is mounted once per session today — cheap to do right,
  avoids collisions if that ever changes.

### `components/app_prism/R/prism_ui.R`

- Add `tags$script(type="module", src="static/prism-webr.js")`.
- Remove `"xkcd"` from the theme `selectInput` choices — the `xkcd`/`extrafont`
  packages aren't installable for WASM (no font loading in-browser), so this theme
  can never render client-side. Not special-cased at runtime; simply not offered.
- Replace `plotOutput(ns("plot1"))` with plain `div`s for plot container / placeholder
  / error, matching `plotbot-claude`'s DOM structure.
- Replace the bare text input with a chat-bubble message list
  (`div(id=ns("chat-messages"))`) placed above the existing input/send-button row.
  The existing `textInput`/`actionButton` stay; only a message log is added above them.
- Sidebar controls and the code `wellPanel` are unchanged.

### `components/app_prism/R/prism_server.R`

- `get_dataframe()`, `llm_plotcode()`, `update_plotcode()`, `get_plotcode()` remain the
  single source of truth for "current R code," largely unchanged.
- Remove `output$plot1 <- renderPlot({...})` and the `require(...)` conditional package
  loads (nothing evals server-side anymore).
- Add an observer on `get_plotcode()` (fires for both chat-send and sidebar-change)
  that serializes `get_dataframe()` to CSV via `write.csv` and sends
  `session$sendCustomMessage("executeCode", list(code = plotcode, csv = csv))`.
- In the chat-send branch specifically, also emit `addMessage` custom messages: the
  user's text immediately on send, and either a short assistant acknowledgement (e.g.
  "Updated the plot.") or an error bubble if the LLM call or code-extraction fails.
- Update the `msg2` prompt text: replace "Load package libraries if needed" with an
  explicit statement of the fixed package set preloaded in the browser sandbox
  (`ggplot2, ggrepel, ggprism, ggpubr, ggsignif`), mirroring `plotbot-claude`'s system
  prompt constraint, so the LLM doesn't generate `library()` calls for packages that
  don't exist in webR.
- The `"reset"` special case (`msg == "reset"` -> empty plot) is preserved, just routed
  through the same custom-message pipeline instead of `renderPlot`.

## Error handling

- **webR eval failure** (bad LLM code, missing package, syntax error): shown in the
  existing red `#plot-error`-style panel, plus an error-role chat bubble — same pattern
  as `plotbot-claude`'s `showError()`.
- **LLM/network failure**: caught by the existing `tryCatch` around
  `omicsai_call_llm`; surfaced as an error bubble (in addition to the existing
  `showNotification`, where applicable).
- **webR not ready yet**: send button disabled and a status message shown until
  `initWebR()` completes, instead of silently dropping the send.

## Known limitations (accepted)

- No `xkcd` theme (removed from the dropdown).
- Any package the LLM reaches for outside the five pre-installed ones fails cleanly in
  the error panel — not specially detected or blocked ahead of time, just constrained
  via the prompt.
- ~30s cold-start cost for webR initialization on first load, same as the source app.

## Testing plan

Manual verification (no unit tests exist for this module today, and the changed
surface is UI/JS plus one server file):

- Load the PRISM board, confirm the webR init status sequence and that send is
  disabled until ready.
- Run each of the three sample datasets (mtcars, iris, geiger) through a chat prompt;
  confirm the plot panel, code panel, and chat log all update.
- Change sidebar dataset/theme/pointsize/fontsize controls; confirm the existing plot
  is patched and re-rendered without an LLM call.
- Send an intentionally-bad prompt; confirm the error bubble and error panel both
  appear.
- Send `"reset"`; confirm the plot clears.
