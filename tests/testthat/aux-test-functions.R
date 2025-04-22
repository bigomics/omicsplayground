searchTabs <- function(board) {
    tabs <- list(
        biomarker = c("Feature selection", "Feature-set ranking"),
        dataview = c("Gene overview", "Sample QC", "Data table", "Samples", "Contrasts"),
        clustering = c("Heatmap", "PCA/tSNE", "Parallel"),
        featuremap = c("Gene", "Geneset"),
        expression = c("Overview", "Top features", "Volcano by comparison", "Volcano by method"),
        correlation = c("Correlation", "Graph"),
        enrichment = c("Enrichment", "Geneset expression", "Enrichment by comparison", "Volcano by comparison", "Volcano by method"),
        signature = c("Volcano plots", "Enrichment", "Overlap/similarity", "Markers"),
        pathway = c("WikiPathways", "Reactome", "GO graph"),
        intersection = c("Pairwise scatter", "Signature clustering"),
        compare = c("Compare expression", "Foldchange", "Gene Correlation"),
        connectivity = c("FC correlation", "FC Heatmap", "Meta-network"),
        singlecell = c("Cell type", "Mapping", "Markers"),
        wgcna = c("WGCNA", "Modules", "Eigengenes", "Enrichment"),
        mofa = c("Overview", "Response", "Weights", "Enrichment"),
        lasagna = c("Multi-layer model", "Multi-partite graph"),
        deepnet = c("Model training", "Biomarker heatmap", "Gradient vs. foldchange")
    )
    return(tabs[[board]])
}

generate_js_click_code <- function(data_value) {
  js_code <- sprintf(
    "
    const targetElement = document.querySelector('a[data-value=\"%s\"]');
    if (targetElement) {
    targetElement.click();
    }
    ",
    data_value
  )
  return(js_code)
}
