# board.wgcna Component

## What it is

The `board.wgcna` component implements Weighted Gene Co-expression Network Analysis (WGCNA) for OmicsPlayground. It identifies modules of highly correlated genes, relates them to sample traits, and provides comprehensive visualizations for network-based gene screening. This is a powerful systems biology approach for finding gene clusters and biomarker candidates.

## Main functions/classes

### Server Functions
- **`WgcnaBoard()`** - Main server module orchestrating WGCNA analysis
  - Computes or loads pre-computed WGCNA results
  - Manages module selection and trait correlation
  - Coordinates all WGCNA visualization sub-modules

### Computation Functions
- **`compute_wgcna()`** - Performs WGCNA calculation
  - Calls `playbase::pgx.wgcna()` with user parameters
  - Parameters: ngenes, minmodsize, power, gset.filter
  - Returns module assignments, eigengenes, TOM matrix
- **`wgcna()`** - Reactive value storing WGCNA results
  - Auto-loads pre-computed results if available
  - Triggers recomputation when parameters change

### Visualization Modules (in R/ folder)

#### Network Visualizations
- **`wgcna_plot_gdendogram.R`** - Gene dendrogram with module colors
- **`wgcna_plot_TOMheatmap.R`** - Topological overlap matrix heatmap
- **`wgcna_plot_correlation_network.R`** - Interactive network graph of module genes
- **`wgcna_plot_module_graph.R`** - Module-level network visualization

#### Module Analysis
- **`wgcna_plot_module_barplot.R`** - Module size distribution
- **`wgcna_plot_module_heatmap.R`** - Expression heatmap for module genes
- **`wgcna_plot_module_membership.R`** - Membership vs. gene significance
- **`wgcna_plot_module_significance.R`** - Module significance for traits

#### Eigengene Analysis
- **`wgcna_plot_eigengene_clustering.R`** - Hierarchical clustering of eigengenes
- **`wgcna_plot_eigengene_heatmap.R`** - Eigengene expression heatmap
- **`wgcna_plot_MTrelationships.R`** - Module-Trait correlation heatmap

#### Quality Control
- **`wgcna_plot_s_independence.R`** - Scale-free topology fit and connectivity
- **`wgcna_plot_sampledendrogram.R`** - Sample clustering for QC
- **`wgcna_plot_gclustering.R`** - Gene clustering (UMAP/tSNE)

#### Tables and Enrichment
- **`wgcna_table_genes.R`** - Interactive table of module genes
- **`wgcna_table_enrichment.R`** - Pathway enrichment for modules
- **`wgcna_html_module_summary.R`** - AI-generated module summaries

### Specialized WGCNA Types (subdirectories)
- **`consensusWGCNA/`** - Consensus network across multiple datasets
- **`multiwgcna/`** - Multi-condition WGCNA
- **`preservationWGCNA/`** - Module preservation analysis

## Project relationships

### Dependencies
- **WGCNA package** - Core network analysis algorithms
- **playbase::pgx.wgcna()** - OmicsPlayground wrapper for WGCNA
- **components/ui/ui-PlotModule.R** - All visualizations use PlotModuleServer
- **components/ui/ui-TableModule2.R** - Gene and enrichment tables
- **plotly, visNetwork** - Interactive network visualizations

### Used by
- **app/server.R** - Main application integrates WgcnaBoard
- Requires **PGX object** with expression data and sample information

### Data flow
1. PGX object provides expression matrix and sample traits →
2. WGCNA computes correlation network and modules →
3. Results stored in `pgx$wgcna` slot →
4. Visualizations query module genes, eigengenes, TOM →
5. User selects modules/traits for detailed exploration

## Modification guide

### Changing WGCNA parameters
1. Edit default parameters in `wgcna_ui.R`:
   - `ngenes` - Number of genes to include
   - `minmodsize` - Minimum module size
   - `power` - Soft-thresholding power
2. Update `compute_wgcna()` to pass new parameters
3. Adjust UI sliders/inputs for parameter selection

### Adding new visualizations
1. Create new plot module: `wgcna_plot_newviz.R`
2. Define UI: `wgcna_plot_newviz_ui()`
3. Define server: `wgcna_plot_newviz_server(wgcna, pgx, ...)`
4. Call in `wgcna_server.R` with appropriate reactive inputs
5. Add to UI in `wgcna_ui.R` within appropriate tab panel

### Customizing module colors
- WGCNA assigns colors automatically (blue, turquoise, brown, etc.)
- To customize: modify `numericlabels` parameter or post-process
- Colors are used consistently across all visualizations

### Adding new enrichment methods
1. Edit `wgcna_table_enrichment_server()`
2. Add new genesets to enrichment analysis
3. Update table columns to display new results
4. Consider caching for performance

### Optimizing performance
- Pre-compute WGCNA during PGX creation: set in `playbase::pgx.createPGX()`
- Cache TOM matrix (can be large) - stored as `TOM`, `svTOM`, or `wTOM`
- Use `summary=TRUE` for reduced memory footprint
- Consider downsampling genes for very large datasets

### AI Report Architecture (`R/ai.report/`)

Self-contained submodule in `R/ai.report/`. Files sourced alphabetically by `wgcna_ui.R`.

| File | Purpose |
|------|---------|
| ai_report_controls.R | Controls UI/server (mode toggle, module selector, style) |
| ai_report_data_extract.R | Private helpers: per-module extraction core, symbol resolution, enrichment selection, artifact detection |
| ai_report_data.R | Pure data functions (no Shiny, no LLM): `wgcna_build_report_tables`, `wgcna_build_summary_params`, `wgcna_rank_modules`, `wgcna_build_methods` |
| ai_report_server.R | Thin coordinator — single entry point for wgcna_server.R |
| ai_report_ui.R | Entry points + 3-card layout (text \| diagram + infographic) |
| ai_text_server.R | Shiny module only: `wgcna_ai_text_server` (summary & report modes) |
| wgcna_diagram_helpers.R | Diagram/infographic prompt builders and style |

**Data flow:**
```
extract_module_data()          ← shared per-module extraction (ai_report_data_extract.R)
    ▲                ▲
    │                │
wgcna_build_report_tables()   wgcna_build_summary_params()   (ai_report_data.R)
    │                │
    └──── both consumed by ────┘
              │
    wgcna_ai_text_server()     ← Shiny module (ai_text_server.R)
              │
    wgcna_ai_report_server()   ← coordinator (ai_report_server.R)
```

**Board integration** is two calls:
- UI: `wgcna_ai_report_ui(ns("ai_report"))` (tab) + `wgcna_ai_report_inputs_ui(ns("ai_report"))` (sidebar)
- Server: `wgcna_ai_report_server("ai_report", wgcna, pgx, session, watermark)`

**To replicate in another board**:
1. Copy `R/ai.report/` to the new board
2. Replace `extract_module_data()` with board-specific data extraction
3. Replace `wgcna_build_summary_params()` / `wgcna_build_report_tables()` with board-specific formatting
4. Rename coordinator entry point (e.g., `enrichment_ai_report_server()`)
5. Keep controls, layout, diagram, and infographic unchanged

**Note**: multiwgcna and consensusWGCNA boards will have their own board-specific extraction and report templates.

### Important considerations
- **Memory**: TOM matrices are large - use sparse representations when possible
- **Computation time**: WGCNA is slow for >10k genes - provide progress feedback
- **Module interpretation**: Grey module = unassigned genes (exclude from some analyses)
- **Soft threshold**: Critical parameter - visualize scale-free topology fit
- **Sample size**: WGCNA needs sufficient samples for reliable correlations
- **Batch effects**: Correct before WGCNA or correlations will be spurious

### Key reactive values
- `wgcna()` - Complete WGCNA results object
- `input$selected_module` - Currently selected module for detailed views
- `input$selected_trait` - Currently selected trait for correlations
- `enrichTableModule()` - Enrichment results for selected module

### WGCNA result structure
```r
wgcna$me.genes       # List of genes per module
wgcna$datME          # Module eigengenes (samples x modules)
wgcna$datTraits      # Sample traits matrix
wgcna$TOM            # Topological overlap matrix
wgcna$geneTree       # Gene dendrogram
wgcna$moduleColors   # Module assignment per gene
wgcna$power          # Soft-thresholding power used
```
