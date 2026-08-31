##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Methylome Profiler shell. Same construction as qsee_ui(): a bigdash page
## with one sidebarItem per screen and one bigTabItem holding it. The screens
## themselves are the methylome_ui_*() functions below, one per tab, and the
## settings panels live in their own *_inputs() functions declared alongside
## them in the tab. Every panel is a PlotModule or TableModule, so info/methods,
## maximize, and the png/pdf/svg/csv downloads come from the platform rather
## than being reimplemented here.

## Screens size themselves off the flex chain below, not off the viewport.
## Rows stay "auto" for the alert and 1 for the content, otherwise
## layout_columns splits the height evenly and the alert eats a whole share.
##
## This used to be calc(100vh - 96px) per screen plus a 28px top pad, tuned by
## hand against the old shell. Both are gone: the platform now sets
## `.big-tab .bslib-grid.html-fill-item { height: 100% !important }`, which
## overrode the calc anyway (it measured 557px in a 900px window, not the 804px
## the arithmetic promised), and the app container no longer sits above the
## viewport, so the pad was reclaiming nothing.
MP_SUBTAB_HEIGHT <- "100%"

## Every tab is the same shape: one flex wrapper around one layout_columns.
## The wrapper is what mp_fill_css() below hangs the fill chain on.
mp_tab <- function(...) shiny::div(...)

## bigTabs renders as a plain block that sizes to its content, so without this
## the whole app stopped at its natural height and left the bottom of the
## window empty - 661px of a 900px viewport. Same fill chain qsee_ui() carries
## for the same shell, extended one level to the per-screen navset_tab, which
## sits deeper than the platform's `.fullheight-page > div > .tabbable` rules
## reach. Written against the id so it cannot leak into the Dashboard's tabs.
mp_fill_css <- function(id) {
  shiny::tags$style(shiny::HTML(sprintf("
    #%1$s-big-tabs { height: 100%%; display: flex; flex-direction: column; }
    #%1$s-big-tabs > .big-tab:not(.d-none) {
      flex: 1 1 auto; min-height: 0; display: flex; flex-direction: column;
    }
    #%1$s-big-tabs > .big-tab > div {
      flex: 1 1 auto; min-height: 0; display: flex; flex-direction: column;
    }
    #%1$s-big-tabs .tabbable { flex: 1 1 auto; min-height: 0; }
  ", id)))
}

methylome_ui <- function(id = "methylome") {
  ns <- shiny::NS(id)

  shiny::tagList(
  mp_fill_css(id),
  bigdash::bigPage(
    id = id,
    navbar = shiny::div(style = "visibility: hidden; display: none;"),
    sidebar = bigdash::sidebar(
      "Methylome Profiler",
      id = id,
      ## Three screens, ordered by what depends on what. Cell composition used
      ## to sit fifth of six while being the prerequisite for two screens above
      ## it - the EWAS fit and intrinsic age acceleration both refuse without
      ## it and tell the user to go back for it - so it now precedes both.
      bigdash::sidebarItem("Data overview", ns("overview-tab")),
      bigdash::sidebarItem("Sample profiles", ns("profiles-tab")),
      bigdash::sidebarItem("EWAS", ns("ewas-tab"))
    ),
    ## tabSettings need a settings container to render into; without this the
    ## right-hand settings panel never appears. qsee_ui() declares the same.
    settings = bigdash::settings("Settings", id = id),
    bigdash::bigTabs(
      id = id,
      bigdash::bigTabItem(ns("overview-tab"),
        methylome_overview_inputs(id), methylome_ui_overview(id)),
      bigdash::bigTabItem(ns("profiles-tab"),
        methylome_profiles_inputs(id), methylome_ui_profiles(id)),
      bigdash::bigTabItem(ns("ewas-tab"),
        methylome_ewas_inputs(id), methylome_ui_ewas(id))
    )
  ))
}

mp_panel_qc <- function(ns) {
  bslib::nav_panel(
    title = "Sample QC", value = "qc",
    bslib::layout_columns(
    col_widths = 12,
    height = MP_SUBTAB_HEIGHT,
    row_heights = list("auto", 1),
    bs_alert("Per-sample quality and identity for a methylation cohort. Nothing on this screen needs a contrast; it answers what is true of each sample."),
    bslib::layout_columns(
      height = "100%",
      col_widths = c(7, 5),
      methylome_table_ledger_ui(
        ns("ledger_tbl"),
        title = "Per-sample checks",
        info.text = "One row per sample. Bimodality is the fraction of probes outside the intermediate 0.3-0.7 band and sits near 0.9 in a healthy methylome - a low value usually means failed normalisation or degraded DNA. Imprint drift is the mean absolute departure of imprinted DMRs from their expected 0.5. Sex check compares the sex predicted from the X/Y probes against a recorded sex column, and MISMATCH is the cheapest sample-swap signal an array gives you; it reads 'no record' when the sheet does not say. Verdict is relative, not absolute: CHECK marks a sample more than 3 MADs from the cohort median on drift or bimodality, or a sex mismatch - so PASS means 'not an outlier in this cohort', and under five samples nothing is ever flagged. DNAm age is Horvath, read off the stored clock fit; the Epigenetic age screen is where clocks are chosen and withheld below a coverage floor. Computed by playbase.epigenetics::sample_qc(). Sex is predicted from median X and Y beta and is withheld entirely unless the dataset carries at least 100 sex-chromosome probes, since GEO matrices are often deposited with filterXY already applied; the DNAm age column comes from the stored methylclock fit.",
        caption = "Per-sample quality, identity and epigenetic age.",
        height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
      ),
      methylome_plot_betadist_ui(
        ns("ledger_dens"),
        title = "Beta distribution",
        caption = "Density of beta values per sample.",
        info.text = "Beta covers [0,1] and should be strongly bimodal: most probes sit near 0.2 (unmethylated) or near 0.8 (methylated), with little density between. Every sample should trace nearly the same curve. One that flattens out, or whose modes sit noticeably inward of the rest, has usually failed normalisation or come from degraded DNA - and it is the same sample the Bimodality column flags. This is a shape check, not a quantitative one: a real differential-methylation effect is far too small to be visible here.",
        info.methods = "Per-sample kernel density over a fixed random subsample of 6,000 probes, drawn with a common seed so the curves stay comparable between redraws. At most 24 samples are drawn; beyond that the panel is unreadable and the Bimodality column is the better summary. Dashed guides at beta 0.2 and 0.8. Densities are computed on beta rather than M-values, since the two modes this panel checks for are a property of the bounded beta scale. Computed with stats::density() on beta; the per-sample checks beside it come from playbase.epigenetics::sample_qc().",
        height = c("100%", "700px"), width = c("auto", "100%")
      )
    )
  ))
}

## Structure before testing: whatever dominates PC1 is what the EWAS will
## report unless it is adjusted for. The scatter says how the samples fall
## apart; the heatmap beside it says what they fall apart by.
mp_panel_structure <- function(ns) {
  bslib::nav_panel(
    title = "Sample structure", value = "structure",
    bslib::layout_columns(
    col_widths = 12,
    height = MP_SUBTAB_HEIGHT,
    row_heights = list("auto", 1),
    bs_alert("Unsupervised structure, before any contrast is chosen. If a technical variable explains a leading component, adjust for it in the EWAS model - or accept that the result partly measures it."),
    bslib::layout_columns(
      height = "100%",
      col_widths = c(7, 5),
      methylome_plot_clustpca_ui(
        ns("overview_pca"),
        title = "Sample clustering",
        label = "a",
        caption = "PCA, t-SNE and UMAP of the cohort, coloured by phenotype.",
        info.text = "Every sample as a point, positioned by its whole methylome rather than by any single probe. Points that sit together have similar methylomes. Colour by a phenotype and ask whether the grouping you intend to test is the grouping the data already shows: when it is not, the leading structure belongs to something else, and the heatmap beside this panel names it. Switch between PCA, t-SNE and UMAP in the plot options - t-SNE and UMAP separate clusters more readably but their distances between clusters mean nothing, while PCA distances and axis percentages are quantitative. The Arrows control overlays a biplot: <b>features</b> draws the CpGs that define the two components on screen, labelled by probe and gene, and <b>phenotypes</b> draws the direction in which each sample variable increases - so a batch or a sex effect can be read as a direction across the cloud, not just as a colour.",
        info.methods = "PCA is computed in this panel, on M-values: the 1,000 most variable probes, rows centred, an exact singular value decomposition, and the percentage on each axis is that component's share of the variance of those probes. Loadings come from the same decomposition, which is what makes the feature arrows possible; the phenotype arrows are the correlation of each one-hot encoded sample variable with each component, and a continuous variable is split at its median before encoding, so its arrow reads as \"towards the high half\". Arrow lengths are rescaled to fill the cloud - only direction and relative length carry meaning. t-SNE and UMAP are computed here too, from those same components - which is what the platform does as well, embedding its component reduction rather than the raw matrix - so all three layouts answer to the same feature setting. Both are seeded, so a layout does not reshuffle between redraws. The heatmap beside this panel tests these same components, and the <b>Use all features</b> setting drives both. Implemented in mp_pca_components(): base svd() on the centred M-value matrix, with Rtsne::Rtsne() and uwot::umap() run on the resulting component scores for the other two layouts.",
        info.references = list(
          list(
            "Du P, Zhang X, Huang CC, et al. (2010) Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis. BMC Bioinformatics 11:587.",
            "https://doi.org/10.1186/1471-2105-11-587"
          ),
          list(
            "Melville J (2024) uwot: The Uniform Manifold Approximation and Projection (UMAP) Method for Dimensionality Reduction.",
            "https://doi.org/10.32614/CRAN.package.uwot"
          )
        ),
        info.extra_link = NULL,
        height = c("100%", "700px"), width = c("auto", "100%"),
        parent = ns(NULL)
      ),
      methylome_plot_pccovar_ui(
        ns("overview_pccov"),
        title = "Components vs sample variables",
        caption = "What each principal component is made of.",
        info.text = "One test per cell: does this sample variable explain this principal component? A dark cell high in the grid is the finding that matters - it means a leading axis of your data is that variable. When it is the phenotype you intend to test, the experiment worked. When it is a slide, a plate, a scan date or a sex column, that variable is confounded with everything you are about to report, and belongs in the EWAS model as a covariate. Cell composition is the usual culprit in blood, and the Sample profiles screen estimates it. The number in each cell is the share of that component the variable accounts for; a variable can be significant and still explain very little.",
        info.methods = "The same principal components the scatter beside it draws - one decomposition per dataset, so the two panels always report the same PC1. By default those come from the 1,000 most variable probes, which is what the platform stores and what minfi::mdsPlot uses; the <b>Use all features</b> setting decomposes every complete probe instead. That setting matters here: on the demo cohort the default cut leaves cell type and sex untouched but drops the slide (batch) association from p 8e-09 to not significant, because batch is a small shift over very many probes rather than a high-variance effect. Columns tested are those the EWAS accepts as covariates. A numeric column with more than ten distinct values is treated as continuous and tested by Pearson correlation; anything else is a grouping of two to twelve levels tested by Kruskal-Wallis. Shading is -log10 p capped at 16, the printed number is the variance of the component explained (r-squared, or eta-squared for a grouping), and bold marks Benjamini-Hochberg q < 0.05 across the whole grid. This panel is independent of the Layout control: t-SNE and UMAP have no components for a variable to be tested against, so the rows here are always principal components and never the axes on the plot.",
        height = c("100%", "700px"), width = c("auto", "100%")
      )
    )
  ))
}

mp_panel_age <- function(ns) {
  bslib::nav_panel(
    title = "Epigenetic age", value = "age",
    bslib::layout_columns(
    col_widths = 12,
    height = MP_SUBTAB_HEIGHT,
    row_heights = list("auto", 1.15, 1),
    bs_alert("Epigenetic clocks for every sample. Pick the clocks in the settings panel on the right, by family or one at a time. A clock computed on a fraction of its probes returns a confident wrong number, so anything below the coverage floor is withheld rather than shown."),
    bslib::layout_columns(
      height = "100%",
      col_widths = c(7, 5),
      methylome_plot_agecor_ui(
        ns("age_cor"),
        title = "DNAm age vs chronological age",
        caption = "Each clock against the reported age of the sample.",
        info.text = "The headline check on a clock: how closely it tracks the age actually recorded for each sample. Points on the dashed identity line are samples whose methylome matches their years; the solid line is the clock's calibration in this cohort. Read r and MAE together - a clock reading ten years high on every sample still correlates at 0.99, and only the median absolute error shows it. A tight line with a large offset means the clock works but was trained on another tissue; scatter with no line means it has nothing to say here.",
        info.methods = "Up to four of the selected clocks, each against the first sample column that looks like a chronological age - numeric, named age or chron, ignoring anything named for DNAm age, acceleration or a group, and needing at least three values that vary. Dashed line is parity, solid line the least-squares fit, and each title carries the Pearson correlation and the median absolute error in years; median rather than mean, so one failed sample does not set the headline number. The clocks themselves are fitted by methylclock at dataset creation (Pelegi-Siso 2021, Bioinformatics 37:1759); the coverage floor only decides which are shown.",
        height = c("100%", "700px"), width = c("auto", "100%")
      ),
      methylome_plot_clocks_ui(
        ns("age_clocks"),
        title = "Clock agreement",
        caption = "Every selected clock against every other.",
        info.text = "Different clocks were trained on different tissues, different CpG sets and, in some cases, a different target: Horvath and Hannum predict years, PhenoAge was trained on mortality risk and TL on telomere length. Where clocks of the same kind agree, the estimate is trustworthy; where they diverge, the tissue match or the normalisation deserves a look. Disagreement between a chronological clock and a biological one is expected and is not a fault - it is roughly what age acceleration measures.",
        info.methods = "Pairwise scatter of every selected clock that passed the coverage floor, so the panel needs at least two usable ones. Identity line dotted in the lower panels, Pearson correlation printed in the upper panels. All clocks come from the same fit, so the correlations are across samples, never across probe sets. Ages come from methylclock, fitted at dataset creation through playbase.epigenetics::compute_clocks() with coefficients from methylclockData; correlations are stats::cor().",
        height = c("100%", "700px"), width = c("auto", "100%")
      )
    ),
    bslib::layout_columns(
      height = "100%",
      col_widths = c(5, 7),
      methylome_plot_agegroup_ui(
        ns("age_group"),
        title = "Age acceleration by phenotype",
        caption = "Residual of DNAm age on chronological age, split by group.",
        info.text = "Age acceleration is what is left of the epigenetic age once chronological age is accounted for: positive means the methylome looks older than the person is. This is the point where a per-sample number becomes a comparison. It inherits everything the clock inherits, so a group difference in blood composition shows up here as accelerated ageing unless the adjustment is ticked - which is why that box exists. Pick the clock and the phenotype in the settings panel.",
        info.methods = "Residuals of a linear fit of the selected clock on chronological age, boxplotted against the selected phenotype with samples jittered over it, a two-sample t-test for two levels and a Kruskal-Wallis test otherwise. Ticking cell-composition adjustment adds the estimated proportions to that fit, giving intrinsic acceleration in the sense of Chen 2016 (Aging 8:1844); the largest proportion is dropped so the design stays full rank. Without a chronological age column no residual can be taken and the raw DNAm age is shown instead, which is not acceleration. Ages come from methylclock via playbase.epigenetics::compute_clocks(); the fit is stats::lm() and the tests stats::t.test() and stats::kruskal.test().",
        height = c("100%", "700px"), width = c("auto", "100%")
      ),
      methylome_table_coverage_ui(
        ns("age_cov"),
        title = "Clocks and coverage",
        info.text = "What each selected clock is, its median estimate across the cohort, and whether this dataset carries enough of its CpGs to compute it at all. Two columns say what is lost: the fraction of the clock's CpGs absent here, and the fraction of its total absolute model weight those CpGs carried. They diverge - a clock can lose 5% of its probes and 40% of its weight if the missing ones are the heavy ones - which is why the ageing literature reports both (Lussier 2024, Clin Epigenetics 16:166). A clock below the coverage floor reads 'withheld' rather than being estimated from partial data, because a clock fitted on a fraction of its probes returns a confident wrong number. A dash means the coefficients are not reachable for that clock. Both the floor and the clock selection are display filters over a fit that used neither. Coverage is computed against each clock's published coefficient set in methylclockData, the same source methylclock fits from.",
        caption = "Per-clock coverage and whether the estimate is usable.",
        height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
      )
    )
  ))
}

mp_panel_context <- function(ns) {
  bslib::nav_panel(
    title = "Genomic context", value = "context",
    bslib::layout_columns(
    col_widths = 12,
    height = MP_SUBTAB_HEIGHT,
    row_heights = list("auto", 1),
    bs_alert("Where methylation sits in the genome. Both panels read the CpG island and gene-position annotation that playbase produces for every methylation dataset."),
    bslib::layout_columns(
      height = "100%",
      col_widths = c(6, 6),
      methylome_plot_context_ui(
        ns("char_island"),
        title = "Mean beta by relation to CpG island",
        caption = "Methylation across islands, shores, shelves and open sea.",
        info.text = "CpG islands are normally unmethylated while the surrounding shores, shelves and open sea are progressively more methylated, so the bars should fall away from the island in a clear staircase. A flat profile means something has gone wrong upstream - failed normalisation, or an annotation that did not match the array. This is a property of the cohort as a whole, not of any contrast; it is the same shape in a tumour and in healthy tissue, and only the depth of the staircase moves.",
        info.methods = "Per-probe mean beta across all samples, averaged within each Relation_to_Island category from the array annotation. Probes carrying a semicolon-separated list are assigned their first entry. Bars are coloured by the beta value itself, on the same hypo-to-hyper scale as the ideogram, and dotted guides sit at 0.2 and 0.8. Categories are the Relation_to_Island field of the Illumina manifest, read through the IlluminaHumanMethylation450k/EPIC annotation packages.",
        height = c("100%", "700px"), width = c("auto", "100%"), label = "a"
      ),
      methylome_plot_context_ui(
        ns("char_gene"),
        title = "Mean beta by position in gene",
        caption = "Methylation across promoter, exon, body and UTR.",
        info.text = "Promoters (TSS1500, TSS200, 5'UTR, 1st exon) are normally unmethylated while gene bodies and 3'UTRs are methylated, and this promoter-to-body gradient is the first thing a methylation analyst checks. It is also the reason the EWAS screen tests promoter and body separately rather than averaging a whole gene: the two halves sit at opposite ends of this plot, so a whole-gene average cancels. An absent gradient points at the annotation rather than at the biology.",
        info.methods = "Per-probe mean beta across all samples, averaged within each UCSC RefGene group from the array annotation. Probes annotated to several groups are assigned their first, so each probe contributes once. Probes with no gene position are pooled as 'unannotated'; that bar is usually intergenic open sea and is expected to be high. Groups are the UCSC_RefGene_Group field of the Illumina manifest, read through the IlluminaHumanMethylation450k/EPIC annotation packages.",
        height = c("100%", "700px"), width = c("auto", "100%"), label = "b"
      )
    )
  ))
}

mp_panel_chromosomes <- function(ns) {
    bslib::nav_panel(
      title = "Chromosomes", value = "chromosomes",
      bslib::layout_columns(
        col_widths = 12,
        height = MP_SUBTAB_HEIGHT,
        row_heights = list("auto", 1),
        bs_alert("Beta per chromosome across the cohort. Autosomes should sit level with each other; the sex chromosomes are expected to differ, and a lone autosome standing apart points at copy number, a failed array region or an annotation mismatch."),
        bslib::layout_columns(
          height = "100%",
          col_widths = 12,
          epigenomics_plot_boxplot_beta_ui(
            id = ns("boxplotBeta"),
          title = "Beta chromosomal profiles",
          label = "",
          caption = "Boxplots of beta values per chromosome.",
          info.text = "Beta per chromosome across the cohort. A normal autosome sits in the same place as its neighbours and the boxes are tight; the sex chromosomes are expected to differ, and X is the one chromosome where a wide box is informative rather than alarming. A single autosome standing apart is worth chasing - copy number, a failed array region, or an annotation mismatch. With a phenotype selected the boxes split by group, which is where a global methylation shift between groups would show.",
          info.methods = "For each chromosome the mean beta is computed within each sample first, so a box shows the spread of per-sample means rather than the spread of individual probes - which is why the boxes are far narrower than the beta distribution on the Sample ledger. Probes whose annotation names a centromere are dropped. Which chromosomes are drawn, which samples are excluded and which phenotype splits the boxes are all set in the settings panel; the panel menu also removes X and Y.",
          info.references = NULL, info.extra_link = NULL,
            height = c("100%", "600px"), width = c("100%", "100%")
          )
        )
      )
    )
}

mp_panel_ideograms <- function(ns) {
    bslib::nav_panel(
      title = "Ideograms", value = "ideograms",
      bslib::layout_columns(
        col_widths = 12,
        height = MP_SUBTAB_HEIGHT,
        row_heights = list("auto", 1),
        bs_alert("Methylation drawn along the chromosomes themselves. Pick up to six chromosomes in the settings panel; with a two-level phenotype selected the lower track becomes the difference between the groups."),
        epigenomics_plot_methylIdeogram_ui(
          id = ns("methylIdeogram"),
          title = "Methylation ideograms",
          label = "",
          caption = "Chromosome ideogram of average beta per genomic region.",
          info.text = "Beta laid out along the chromosomes themselves, so large-scale features - hypomethylated blocks, whole-arm shifts, a differentially methylated region large enough to see - are visible in a way a Manhattan or a density plot cannot show. The trend line is what to read: the cloud of individual CpGs is bimodal everywhere and says little on its own. With a two-level phenotype selected, the lower track becomes the difference between the groups, and its axis is autoscaled, so check the numbers before reading a band as large.",
          info.methods = "Drawn with karyoploteR (Gel & Serra 2017, Bioinformatics 33:3088) against hg19, or mm10 for mouse. Beta and annotation are restricted to their shared CpGs, and probes with no coordinate or a position outside the reference chromosome length are dropped. Panel 1 is the per-CpG mean beta across the retained samples, subsampled to 200,000 points, coloured blue below 0.2, yellow between and red above 0.8, with dashed guides at those values. The trend over it bins CpGs in fixed 1 Mb windows and smooths them with a LOESS weighted by the CpGs per bin; it is pooled across samples unless the per-phenotype average is ticked in the panel menu, which splits it in two. Panel 2 is the per-bin difference between the two groups when a two-level phenotype is selected, and a log-scaled probe-density track otherwise.",
          info.references = NULL, info.extra_link = NULL,
          height = c("100%", "600px"), width = c("100%", "100%")
        )
      )
    )
}

mp_panel_composition <- function(ns) {
  bslib::nav_panel(
    title = "Cell composition", value = "composition",
    bslib::layout_columns(
    col_widths = 12,
    height = MP_SUBTAB_HEIGHT,
    row_heights = list("auto", 1, 1),
    bs_alert("Estimated cell proportions, projected onto a reference panel. Cell composition is the dominant confounder in bulk-tissue methylation, so this screen exists as much to supply covariates to the EWAS model as to be read on its own."),
    bslib::layout_columns(
      height = "100%",
      col_widths = c(7, 5),
      methylome_plot_composition_ui(
        ns("comp_bars"),
        title = "Composition per sample",
        caption = "Estimated proportion of each reference cell type.",
        info.text = "One stacked bar per sample. Samples are ordered by the phenotype picked in the settings, so a difference in composition between groups reads as a visible block rather than being scattered along the axis - and that block is the thing to look for, because it is what will otherwise be reported as differential methylation. A single sample with an extreme profile is either real biology or a projection that failed on a poor-quality array; the Sample ledger will usually say which.",
        info.methods = "Constrained projection of each sample's betas onto a reference panel of cell-type-discriminating CpGs (Houseman 2012, BMC Bioinformatics 13:86), run through methylclock's meffilEstimateCellCountsFromBetas. The reference panel is chosen in the settings and has to match the tissue: the blood panels are for whole blood, the cord panels for neonatal samples, with saliva and brain panels alongside. Proportions are not forced to sum to exactly one, so small deviations are expected; they are estimates from an external panel, not measurements of this cohort.",
        height = c("100%", "700px"), width = c("auto", "100%")
      ),
      methylome_table_composition_ui(
        ns("comp_tbl"),
        title = "Estimated proportions",
        info.text = "The per-sample proportions behind the stacked bars, as numbers - the panel to read when one bar looked wrong, and the CSV to take away for a model fitted elsewhere. Rows are samples, columns are the reference panel's cell types, so the panel chosen in the settings decides which columns exist at all. Values are the raw output of the projection and are not rescaled: a row summing to 0.95 or 1.04 is normal, a row far off that is a failed projection rather than an unusual sample. Inside the app these same estimates are what the EWAS 'Adjust for cell composition' box and the intrinsic age-acceleration option consume. Estimated by methylclock::meffilEstimateCellCountsFromBetas(), meffil's implementation of the Houseman constrained-projection method, onto the selected reference panel. Rows sum to approximately one by construction, so the values are compositional and not independent of one another.",
        caption = "Per-sample cell-type proportions.",
        height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
      )
    ),
    methylome_plot_compgroup_ui(
      ns("comp_group"),
      title = "Composition by phenotype",
      caption = "Each cell type compared across the groups of the contrast.",
      info.text = "This is the confounding check. If composition differs between the groups being compared, an unadjusted differential-methylation result is largely a picture of that difference rather than of the phenotype - most published blood CpG-age associations lose significance once composition is adjusted for (Jaffe & Irizarry 2014, Genome Biol 15:R31). Anything that shifts the neutrophil-to-lymphocyte ratio does it: infection, obesity, smoking, age. A clear difference here is an instruction to tick the adjustment in the EWAS model, not a result to report.",
      info.methods = "Boxplot of each estimated proportion against the selected phenotype, samples jittered over it, with a two-sample t-test for two-level phenotypes and a Kruskal-Wallis test otherwise. The p-values are per cell type and are not corrected across the panel, and the proportions are not independent of each other - they are read as a warning flag, not as a result. Proportions are estimated by methylclock::meffilEstimateCellCountsFromBetas(), a Houseman-style constrained projection onto a sorted-cell reference; the tests are stats::t.test() and stats::kruskal.test().",
      height = c("100%", "700px"), width = c("auto", "100%")
    )
  ))
}

## The three screens that answer "what is in this dataset" were three sidebar
## entries: Sample ledger, Methylome character, Methylome landscape. None of
## them needs a contrast, a button or a fitted model - they render on load,
## which is exactly what separated them from every other screen - and two of
## the three had no settings panel at all, so the right-hand drawer rendered
## empty on a third of the app. They are one screen with four sub-tabs.
methylome_ui_overview <- function(id = "methylome") {
  ns <- shiny::NS(id)
  mp_tab(bslib::navset_tab(
    id = ns("overview_subtab"),
    mp_panel_qc(ns),
    mp_panel_structure(ns),
    mp_panel_context(ns),
    mp_panel_chromosomes(ns),
    mp_panel_ideograms(ns)
  ))
}

## Composition first, deliberately: it is the prerequisite for the intrinsic
## option on the age sub-tab, and for the EWAS model's cell adjustment.
methylome_ui_profiles <- function(id = "methylome") {
  ns <- shiny::NS(id)
  mp_tab(bslib::navset_tab(
    id = ns("profiles_subtab"),
    mp_panel_composition(ns),
    mp_panel_age(ns)
  ))
}

## Merged screens need their controls gated to the sub-tab that owns them,
## or the drawer shows every screen's settings at once. Same helper shape as
## methylome_ewas_inputs()'s on_tab().
methylome_overview_inputs <- function(id) {
  ns <- shiny::NS(id)
  bigdash::tabSettings(
    shiny::conditionalPanel(
      condition = "input.overview_subtab == 'structure'",
      ns = ns,
      withTooltip(
        shiny::checkboxInput(ns("structure_allfeatures"), "Use all features", FALSE),
        "Off, the components come from the 1,000 most variable probes - the convention, and what the platform stores. On, every complete probe is decomposed: slower, and the only way to see a batch effect, which is a small shift spread over very many probes rather than a high-variance one.",
        placement = "left"
      )
    ),
    shiny::conditionalPanel(
      condition = "input.overview_subtab == 'chromosomes' || input.overview_subtab == 'ideograms'",
      ns = ns,
      methylome_landscape_inputs(id, wrap = FALSE)
    )
  )
}

methylome_profiles_inputs <- function(id) {
  ns <- shiny::NS(id)
  bigdash::tabSettings(
    shiny::conditionalPanel(
      condition = "input.profiles_subtab == 'composition'", ns = ns,
      methylome_deconv_inputs(id, wrap = FALSE)
    ),
    shiny::conditionalPanel(
      condition = "input.profiles_subtab == 'age'", ns = ns,
      methylome_age_inputs(id, wrap = FALSE)
    )
  )
}

methylome_ui_ewas <- function(id = "methylome") {
  ns <- shiny::NS(id)
  mp_tab(bslib::navset_tab(
    ## id + explicit values so the settings panel can show only the
    ## controls that belong to the visible sub-tab.
    id = ns("ewas_subtab"),

    bslib::nav_panel(
      title = "Manhattan & hits", value = "hits",
      bslib::layout_columns(
        col_widths = 12,
        height = MP_SUBTAB_HEIGHT,
        row_heights = list("auto", 1.05, 1.2),
        bs_alert("Differential methylation for the model fitted in the settings panel, shown genome-wide. Set the threshold on the right; the lines, the coloured points, the heatmap and the table all follow it."),
        bslib::layout_columns(
          height = "100%",
          col_widths = c(7, 5),
          methylome_plot_manhattan_ui(
            ns("ewas_manhattan"),
            title = "Manhattan",
          caption = "Genome-wide significance by chromosome position.",
          info.text = "Every tested feature at its genomic position against significance. Features passing the threshold are coloured by direction, red where methylation increases and blue where it decreases; an F-test across three or more groups has no direction and gets a single accent colour instead. Read the shape as well as the peaks - a real EWAS is a mostly flat field with a few towers, while signal smeared across whole chromosomes is confounding or batch, and the QQ panel is where to confirm it. On a gene-region fit each point is a region, drawn at the median position of its probes.",
            info.methods = "Features are ordered by chromosome and position from the array annotation; the stored chromosome field is a cytoband, so the arm and band are stripped first, and anything without a chromosome or a position is dropped. Alternating grey shades separate chromosomes. A q-value cut-off has no fixed position on a -log10(p) axis, so the dashed line is drawn at the largest p-value that still passes the current threshold - including the minimum |delta-beta| filter, so the line always sits exactly where the coloured points begin. When nothing passes, no line is drawn and the corner says so instead of printing a hit count. Positions come from the Illumina manifest; the model behind the p-values is limma::lmFit() with eBayes(trend = TRUE) on M-values.",
            height = c("100%", "700px"), width = c("auto", "100%")
          ),
          methylome_plot_volcano_ui(
            ns("ewas_volcano"),
            title = "Volcano",
            caption = "Effect size against significance.",
            info.text = "The view a Manhattan cannot give: how large the change is, not just how significant. A CpG can be overwhelmingly significant and move methylation by half a percent - in a population EWAS that is the normal case rather than a disappointment, so read the x axis before the y. Points are coloured by direction once they pass the threshold and the strongest are labelled by gene. There is no volcano for an F-test across three or more groups: its effect is an unsigned range, and drawing it here would imply a direction the test never established.",
            info.methods = "Delta beta on the x axis against -log10(p) on the y. The model is fitted on M-values, which is the appropriate scale for testing, and the effect is reported on the beta scale, which is the one that means anything (Du 2010, BMC Bioinformatics 11:587). For a two-group contrast delta beta is a plain difference of group means, unadjusted by the covariates, so it can sit beside a much weaker p-value once they are added. For a continuous outcome it is the slope per unit of the exposure, taken from the same design refitted on the beta scale, and the axis says so. The horizontal line follows the same largest-passing-p rule as the Manhattan; vertical dashed lines appear only when a minimum |delta-beta| is set. Up to eight of the most significant hits are labelled. Fitted with limma::lmFit() and eBayes(trend = TRUE) on M-values.",
            height = c("100%", "700px"), width = c("auto", "100%")
          )
        ),
        bslib::layout_columns(
          height = "100%",
          col_widths = c(7, 5),
          methylome_table_hits_ui(
            ns("ewas_hits"),
            title = "CpGs passing the threshold",
            info.text = "Every feature called at the current cut-off, most significant first. Delta beta is the change in methylation itself, not the M-value log fold change the model was fitted on. For a categorical outcome with three or more levels the model is an F-test, which has no direction: the column becomes an unsigned beta range, Direction reads '-', and the F statistic is shown. Studies and Reported carry the EWAS Catalog evidence; 'not reported' means absent from that catalog, which is not the same as novel biology. When Test at is a gene region, rows are GENE:promoter or GENE:body rather than CpGs, delta beta is a difference of medians across the region and is not comparable to a probe-level delta beta, and Probes shows how many probes the median was taken over - it ranges from 3 to over a thousand. See the Test at tooltip for how the regions are defined. Fitted with limma::lmFit() and eBayes(trend = TRUE) on M-values; delta beta is a raw difference of group means, or on a continuous outcome the slope from that same design refitted on the beta scale. Annotation columns come from the Illumina manifest, trait columns from the bundled EWAS Catalog release.",
            caption = "Significant CpGs, effect and significance first. Studies and Reported at the end carry the EWAS Catalog evidence: sort Studies ascending for the hits nobody has reported yet; each Reported trait links out to its own catalog page.",
            height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
          ),
          methylome_plot_hitmap_ui(
            ns("ewas_hitmap"),
            title = "Top hits, sample by sample",
            caption = "Beta of the most significant CpGs across the fitted samples.",
            info.text = "The hit list as a picture. Rows are the features the model called, columns are the samples it was fitted on, ordered by group. A real result reads as a block that changes at the group boundary; a scattered pattern with no boundary means the hits are carried by a few samples, and the stripcharts on the QQ & context sub-tab show which. Some separation is guaranteed - these are the features selected for separating the groups - so the question is whether it is clean, and whether it is the same samples row after row.",
            info.methods = "Up to fifty of the most significant features passing the current threshold. Colour is beta on a fixed 0-1 scale, blue unmethylated and red methylated, rather than a row z-score: the absolute level is the interpretable quantity, and a z-score makes a 2% difference look like a 60% one. Nothing is clustered - rows stay in p-value order and columns are ordered by the model's own groups, or by tertile of the exposure for a continuous outcome, with the bar above the matrix carrying that grouping. On a gene-region fit the rows are rebuilt from their member probes with the same median the fit used. Features and their statistics come from the limma fit (limma::lmFit, eBayes(trend = TRUE)); the heatmap itself is drawn with base graphics.",
            height = c("100%", "700px"), width = c("auto", "100%")
          )
        )
      )
    ),

    bslib::nav_panel(
      title = "Regions & pathways", value = "regions",
      bslib::layout_columns(
        col_widths = 12,
        height = MP_SUBTAB_HEIGHT,
        row_heights = list("auto", 1),
        bs_alert("Region-level views of the current model. At probe level: DMR calling and gene-set testing, both run on demand from the settings panel. At gene-region level: how each gene's promoter compares with its body, which is the comparison only that fit can make."),
        ## Which half of this sub-tab is live follows the Test at picker, not the
        ## fitted model, so switching resolution swaps the panels immediately and
        ## the panel's own validate() covers the gap until the refit lands.
        ##
        ## Both modes sit inside one plain div rather than being siblings of the
        ## alert. As direct children of the grid, the hidden one still reserved
        ## its row and left a tall blank band above whichever mode was showing.
        shiny::div(
          style = "height: 100%; min-height: 0;",
        shiny::conditionalPanel(
          condition = "input.ewas_collapse != 'probe'", ns = ns,
          style = "height: 100%;",
          bslib::layout_columns(
            height = "100%",
            col_widths = c(6, 6),
            methylome_plot_divergence_ui(
              ns("ewas_divplot"),
              title = "Promoter against body",
              caption = "Each gene's promoter change against its gene-body change.",
              info.text = "The comparison only a gene-region fit can make. Promoter methylation tends to silence a gene while gene-body methylation tracks with expression instead - the DNA methylation paradox (Jones 1999, Trends Genet 15:34; Yang 2014, Cancer Cell 26:577) - so a gene whose two halves move in opposite directions is a specific finding rather than a stronger version of the same one. Points off the dotted diagonal are genes where promoter and body disagree; the two off-axis quadrants are where they disagree in sign. Red points are significant in both halves and opposing - those are the ones worth reading. Grey points did not pass the threshold in both halves.",
              info.methods = "One point per gene carrying both a promoter and a gene-body feature - about two thirds of genes on a 450K array once the three-probe floor is applied. Axes are the delta beta of each feature, so a difference of medians across that region on the beta scale. Significance is the same threshold set for the rest of the screen, applied to each half independently; a gene counts as opposing when the two delta betas have opposite signs, which is a statement about direction and not about magnitude. The signed difference between the two halves is MeGDP, defined body minus promoter by Li 2019 (Genome Res 29:270), and reported here in that direction. The layout of this plot is not taken from any publication. Note the promoter half inherits the array manifest's RefGene annotation, and over a quarter of Infinium promoter assignments have been reported as questionable (Bizet 2022, Epigenetics 17:2434). An F-test across three or more groups reports an unsigned beta range and is refused rather than compared for direction. Both features come from one limma fit (limma::lmFit, eBayes(trend = TRUE)) over gene regions built on mCSEA's partition of the manifest.",
              height = c("100%", "700px"), width = c("auto", "100%")
            ),
            methylome_table_divergence_ui(
              ns("ewas_divtbl"),
              title = "Genes ranked by promoter-body divergence",
              info.text = "Every gene with both a promoter and a gene-body feature, findings first: significant in both halves, then opposing in sign, then by the size of the gap. Promoter and Body are each region's delta beta on the beta scale, coloured by direction; MeGDP is the signed difference body minus promoter, in the direction Li 2019 (Genome Res 29:270) defines it. Opposing marks the genes whose halves move in opposite directions - the pattern that makes averaging a whole gene meaningless. Probes columns say how many probes each median was taken over, which is worth checking before believing a row: a three-probe promoter and a four-hundred-probe body are not equally well measured. Both effects come from one limma fit (limma::lmFit, eBayes(trend = TRUE)) over gene regions built on mCSEA's partition of the manifest.",
              caption = "Paired promoter and body results per gene.",
              height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
            )
          )
        ),
        shiny::conditionalPanel(
          condition = "input.ewas_collapse == 'probe'", ns = ns,
          style = "height: 100%;",
          bslib::layout_columns(
            col_widths = 12, height = "100%", row_heights = list(1.15, 1),
        methylome_plot_dmrregion_ui(
          ns("ewas_dmrplot"),
          title = "Region detail",
          caption = "Methylation across the selected region. Click a row in the table below.",
          info.text = "A region in the table is only a set of numbers; this is what it looks like. A real DMR is a run of neighbouring CpGs all moving the same way and separating the groups across the whole region. One strong CpG dragging its neighbours along produces a similar table row and an obviously different picture, which is the reason this panel sits above the table rather than below it. Click any row to draw that region; with nothing selected the most significant one is shown.",
          info.methods = "The region view the field draws for a called DMR, following DMRcate's DMR.plot (Peters 2015, Epigenetics Chromatin 8:6) and coMET (Martin 2015, BMC Bioinformatics 16:131). Beta of every CpG in the region plus 30% flanking context on each side, against genomic position. Faint lines are individual samples - only those the model was actually fitted on - and heavy lines the mean of each group, or of the low, middle and top tertile when the outcome is continuous. The shaded box is the called region, so the flanking probes are visibly outside it, and the track under the axis carries each probe's CpG island context. Beta, not M-values, so the vertical distance between the group means is the methylation difference itself. Regions are called by dmrff::dmrff() on the per-CpG estimates and standard errors from the limma fit; the ideogram track is karyoploteR.",
          height = c("100%", "700px"), width = c("auto", "100%")
        ),
        bslib::layout_columns(
          height = "100%",
          col_widths = c(6, 6),
          methylome_table_dmr_ui(
            ns("ewas_dmr"),
            title = "Differentially methylated regions",
            info.text = "Contiguous runs of CpGs that move together, the unit methylation phenotypes occur in - promoter island hypermethylation, hypomethylated blocks - rather than isolated probes. Single-CpG regions are dropped: those are already in the hits table. Effect (M) is on the M-value scale the model was fitted on, so it carries a direction and a magnitude but is not a percentage of methylation; click a row to see the region on the beta scale above. Adjusted p is corrected over the regions tested, not over the probes. Regions are called by dmrff, which meta-analyses each probe's estimate and standard error across its neighbours, reading the methylation matrix only to account for the correlation between them (Suderman 2018, bioRxiv preprint - never journal-published). It needs a signed effect per probe, so an F-test across three or more groups is refused, as is a fit collapsed to gene regions. Called by dmrff::dmrff() from the per-CpG effect estimates and standard errors of the limma fit - not from the p-values, which is what lets it meta-analyse adjacent probes rather than merely count them.",
            caption = "Regions called from the fitted model by dmrff.",
            height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
          ),
          methylome_table_enrich_ui(
            ns("ewas_enrichgs"),
            title = "Gene sets (probe-bias corrected)",
            info.text = "Gene-set enrichment of the hit list, corrected for the number of probes per gene. Genes carry very different probe counts, so an uncorrected test returns the same probe-dense developmental sets whatever the biology (Geeleher 2013, Bioinformatics 29:1851). gometh, from missMethyl (Phipson 2016, Bioinformatics 32:286), reweights with Wallenius' noncentral hypergeometric distribution and also handles CpGs mapping to several genes. The test is on the hit list, so it moves with the threshold - a set that appears and disappears as you nudge the cut-off is not a finding. Background is the probes actually tested, not the whole array. Needs at least ten hits and a probe-level fit; a gene-region fit has already averaged away the probe ids the correction runs on. Tested with missMethyl::gometh(), whose probes-per-gene correction is the point of using it: genes covered by more CpGs are likelier to contain a hit by chance, and an uncorrected test returns the same probe-dense developmental terms on almost any hit list.",
            caption = "GO or KEGG terms enriched among the significant CpGs.",
            height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
          )
        )
          )
        )
        )
      )
    ),

    bslib::nav_panel(
      title = "QQ & context", value = "qq",
      bslib::layout_columns(
        col_widths = 12,
        height = MP_SUBTAB_HEIGHT,
        row_heights = list("auto", 1, 1.35),
        bs_alert("Whether the result holds up: how far the p-values depart from the null, where in the genome the hits concentrate, and what the top hits look like sample by sample."),
        bslib::layout_columns(
          height = "100%",
          col_widths = c(5, 7),
          methylome_plot_qq_ui(
            ns("ewas_qq"),
            title = "QQ and inflation",
            caption = "Observed against expected p-values, with lambda.",
            info.text = "How far the whole p-value distribution departs from the null. Points on the diagonal mean nothing beyond chance. A curve that lifts off early and stays up is the signature of a systematic difference between the groups - unadjusted cell composition, batch, population structure - while a tail departing only at the very top is what genuine signal looks like. Lambda puts one number on that departure, but it rises with real widespread signal as well as with confounding, so it is a prompt to look rather than a pass or fail. Where it is high, the answer is usually another covariate, not a stricter threshold.",
            info.methods = "Observed -log10 p-values against the uniform expectation. Lambda is the median chi-square statistic implied by the p-values divided by its null expectation (Devlin & Roeder 1999, Biometrics 55:997). Where the bacon package is installed, the bias and inflation of an empirical null fitted to the moderated t statistics are printed beside it (van Iterson 2017, Genome Biol 18:19). Both are diagnostics only and nothing is rescaled by them: applying bacon to an unadjusted model would launder confounding into well-calibrated-looking p-values, which is worse than an honest lambda. They are withheld on the F-test branch, where the statistic is an F rather than a signed z and bacon's empirical null does not apply.",
            height = c("100%", "700px"), width = c("auto", "100%")
          ),
          methylome_plot_enrichment_ui(
            ns("ewas_enrich"),
            title = "Where the hits sit",
            caption = "Genomic context of the significant CpGs.",
            info.text = "Whether the significant CpGs concentrate in islands, shores, shelves or open sea. Shore enrichment with island depletion is the classic pattern in differential methylation (Irizarry 2009, Nat Genet 41:178): islands at active promoters are held unmethylated and have little room to move, while the shores flanking them are where tissue and disease differences sit. A hit list piled into open sea instead is usually reporting something other than regulation - SNP-driven probes, batch, or a copy-number difference between the groups.",
            info.methods = "Odds ratio of each Relation_to_Island category among the features passing the threshold against the probes actually tested in this contrast, not against the whole array: masked probes and probes dropped for missing values are not in the background, and using the array instead would report the masking as enrichment. A half-count is added to each cell so an empty category stays finite. Bars are the log2 odds ratio; no p-value is given, because neighbouring probes are not independent. The panel refuses on a gene-region fit - an aggregated feature's context is a majority vote over its probes, and a gene body routinely spans island, shore, shelf and open sea, so the odds ratio would be measuring the vote. Categories are the manifest's Relation_to_Island field; the odds ratio is computed directly with a half-count added to each cell, no package involved.",
            height = c("100%", "700px"), width = c("auto", "100%")
          )
        ),
        methylome_plot_stripcharts_ui(
          ns("ewas_strips"),
          title = "Top CpGs, sample by sample",
          caption = "Per-sample beta for the most significant CpGs, split by group.",
          info.text = "Beta of every sample at each of the top features. This is the panel that shows whether a hit is real or driven by a couple of outliers, and at what absolute level of methylation it sits - a 0.05 difference at beta 0.50 and at beta 0.02 are very different claims, and only one of them survives a second cohort. Look for group clouds that are separated but overlapping; a hit whose groups are perfectly separated by a handful of points is usually a probe problem rather than biology.",
          info.methods = "The most significant features passing the current threshold, one panel each, points jittered within group with a bar at the group mean. For a continuous outcome the panel becomes a scatter against the exposure with its least-squares line, which is the honest picture of whether the slope rides on a few extreme samples. The y axis is fixed to the full [0,1] beta range rather than zooming to the data, so panels stay comparable and the absolute level stays visible; dotted guides at 0.2 and 0.8. On a gene-region fit each panel is the region's median over its probes. How many are drawn is set in the settings panel. Betas are plotted as stored; the ranking and the p-values come from the limma fit (limma::lmFit, eBayes(trend = TRUE)).",
          height = c("100%", "700px"), width = c("auto", "100%")
        )
      )
    ),

    bslib::nav_panel(
      title = "EWAS Catalog", value = "catalog",
      bslib::layout_columns(
        col_widths = 12,
        height = MP_SUBTAB_HEIGHT,
        row_heights = list("auto", 1),
        bs_alert("What your hit list is already known for, against the EWAS Catalog (release 2026-04-09, bundled - no lookup happens over the network). Both panels follow the outcome and threshold set on the other sub-tabs. For a single CpG's own record, follow the link on its row in the hits table."),
        bslib::layout_columns(
          height = "100%",
          col_widths = c(7, 5),
          methylome_plot_traitfreq_ui(
            ns("ewas_traitfreq"),
            title = "Already reported for",
            caption = "Traits your hit list is collectively already associated with.",
            info.text = "What the whole hit list, taken together, is already known for. A smoking EWAS whose top bar is Smoking is a result behaving as it should. A BMI EWAS whose top bar is Age is a confounding warning no other panel here will give you. The subtitle carries the denominator - how many of your hits appear in the catalog at all - and it matters: most CpGs have never been reported for anything, so a short list of bars over a small denominator says little. Bars count CpGs, so one heavily studied CpG cannot dominate the chart on its own.",
            info.methods = "Every CpG passing the threshold is looked up in the EWAS Catalog (Battram 2022, Wellcome Open Res 7:41), release 2026-04-09, bundled with the app - no request leaves the container. Traits are counted by how many of your hits carry them. This is a count and deliberately not an enrichment test: whether a CpG is in the catalog at all reflects how many EWASs have been run and on what, so a p-value against that background would be measuring publication history rather than biology. Looked up with playbase.epigenetics::ewas_catalog_traits() against the release bundled in playdata; no package fetches anything at runtime.",
            height = c("100%", "700px"), width = c("auto", "100%")
          ),
          methylome_table_traits_ui(
            ns("ewas_traits"),
            title = "Every trait, not just the top ones",
            info.text = "The full trait breakdown of the hit list, which the chart can only show the top of. CpGs is how many of your hits are reported for the trait; Studies is how many catalogued study analyses reported those pairs. The two disagree in a useful way: a trait can be heavily studied on one CpG, or lightly studied across hundreds of them, and only the first is evidence that your hit list as a whole overlaps it. Sort or search for a specific trait, and follow the arrow to its page on ewascatalog.org. Counts come from the bundled EWAS Catalog release 2026-04-09 (Battram 2022, Wellcome Open Res 7:41); absence from it means nobody has published on that CpG, not that the association is new. Read through playbase.epigenetics::ewas_catalog_traits() against the release bundled in playdata; no request leaves the container.",
            caption = "Traits ranked by how much of your hit list they already account for.",
            height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
          )
        )
      )
    )
  ))
}
