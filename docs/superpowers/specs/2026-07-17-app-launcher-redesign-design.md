# Smart Tools app-launcher redesign — Design

## Purpose

Redesign the "Smart Tools" page (`components/app_tools/R/tools_ui.R`) from
its current 3-card `bslib::layout_columns` layout into a mobile-style app
launcher: a grid of colorful rectangular tiles, one per tool, similar to a
phone's home screen.

Apps covered today: ID Converter, Qsee/Bsee, SmartPrism. No new apps are
being added as part of this work, but the grid should read reasonably if
more are added later.

## Layout

- Header block unchanged: "Smart Tools" title, subtitle, and the
  `tools_search` search box, all in their current positions/styling.
- Below it, a responsive CSS grid replaces the current
  `bslib::layout_columns(col_widths = 4, ...)`:
  `grid-template-columns: repeat(auto-fill, minmax(160px, 1fr))`, with a
  consistent gap. This yields 4+ tiles per row on wide screens and reflows
  down on narrower ones, and grows gracefully as more apps are added later
  (no code change needed to accommodate a 4th/5th tile).
- Implemented as a new `.app-launcher-grid` container class in
  `scss/components/_toolbox.scss`, applied via `tags$div(class =
  "app-launcher-grid", ...)` in `tools_ui.R` (not `bslib::layout_columns`,
  since the auto-fill/wrap behavior isn't something `layout_columns`'
  fixed `col_widths` gives us).

## Tile design

Each app is one rounded-rectangle tile ("Style A" from the interactive
mockup, approved by user):

- Rounded corners (~18px radius), gradient background, one distinct hue per
  app (not derived from the existing brand palette — chosen for visual
  variety/scannability, phone-icon style). Suggested starting hues: blue
  for ID Converter, teal/green for Qsee/Bsee, purple for SmartPrism.
- `shiny::icon(...)` (FontAwesome) centered near the top of the tile —
  reuses the icon mechanism already used elsewhere in the app (e.g.
  `idconvert_server.R`), no new asset/dependency introduced. The existing
  PNG applet images (`converter.png`, `qsee-bsee.png`, `smartprism.png`)
  are **not** reused in the new tiles.
- Bold app label below the icon.
- One-line description below the label (reuse/trim the existing
  descriptions to fit the tile width).
- Hover state: slight lift + shadow increase, for affordance.
- No separate "Run" button/pill — the entire tile is clickable.

## Interaction

- Implemented as a `shiny::actionButton` sized/styled via CSS to fill the
  whole tile (stripped of default Bootstrap button chrome: border,
  background, padding overridden by the `.app-tile` class), rather than a
  plain clickable `<div>` — this keeps standard Shiny click semantics
  (keyboard focus, `input$x` click counting) for free.
- Input IDs are kept **exactly as they are today** —
  `ns("runtool_idconvert")`, `ns("runtool_qc")`, `ns("runtool_prism")` — so
  `components/app_tools/R/tools_server.R` requires **no changes**.

## Styling

New rules added to `scss/components/_toolbox.scss`:

- `.app-launcher-grid` — the responsive grid container described above.
- `.app-tile` — base rectangular tile: rounded corners, padding, box-shadow,
  hover transition, flex layout for icon/label/description, white text.
- Per-app color variants (e.g. `.app-tile-blue`, `.app-tile-teal`,
  `.app-tile-purple`) applied via an extra class per tile, or via a small
  inline `background` style driven from an R-side color value — whichever
  reads cleaner in `tools_ui.R` at implementation time.
- The existing `.action-pill` rule is dropped from `tools_ui.R`'s markup;
  it will be removed from `_toolbox.scss` too, after confirming (via grep)
  it isn't referenced anywhere else in the codebase.

## Files touched

- `components/app_tools/R/tools_ui.R` — rewritten to build the launcher
  grid instead of the 3-card layout.
- `scss/components/_toolbox.scss` — add `.app-launcher-grid`, `.app-tile`,
  and color-variant rules; remove `.action-pill` if confirmed unused
  elsewhere.
- `components/app_tools/R/tools_ui_classic.R` — **new** file, an untouched
  copy of the current `tools_ui.R` contents with the function renamed to
  `tools_ui_classic`. Not wired into `00SourceAll.R` or called anywhere —
  it exists purely as an easy-to-find, easy-to-diff revert reference, in
  addition to git history (the current `tools_ui.R` is already committed
  cleanly on this branch, so `git checkout` also works).
- `components/app_tools/R/tools_server.R` — **unchanged**.

## Out of scope

- No changes to `tools_server.R`, `runmonitor.R`, or any other page in the
  app.
- No new apps/tiles added to the grid.
- No changes to how/where the Smart Tools page itself is navigated to.
