# Conventions

## OmicsBoardUI layout

When placing a `bslib::layout_columns` directly inside `OmicsBoardUI` (as the board's content row), the inner `layout_columns` **must** have an explicit `height = "100%"`. Without it, the fill chain breaks and the board collapses to content height instead of filling the browser window.

This applies to all boards, including Across Datasets, Qsee/Bsee, etc.