##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Methylome Profiler shell. Same construction as qsee_ui(): a bigdash page
## with one sidebarItem per screen and one bigTabItem holding it. Every panel
## is a PlotModule or TableModule, so info/methods, maximize, and the
## png/pdf/svg/csv downloads come from the platform rather than being
## reimplemented here.

## bigdash pulls its app container 15px above the viewport - it assumes a
## navbar occupies that strip, and this board hides its navbar - so without a
## correction the first element on every tab is clipped by the app top bar.
## Push the content back down, then size the tab to the space that is really
## left. Rows are "auto" for the alert and 1 for the content, otherwise
## layout_columns splits the height evenly and the alert eats a whole share.
MP_TAB_PAD <- "padding-top: 28px;"
MP_TAB_HEIGHT <- "calc(100vh - 40px)"
## Sub-tabbed screens lose another strip of height to the nav_tab chrome.
MP_SUBTAB_HEIGHT <- "calc(100vh - 96px)"

## Every tab is the same shape: a nudged wrapper around one layout_columns.
mp_tab <- function(...) shiny::div(style = MP_TAB_PAD, ...)

methylome_ui <- function(id = "methylome") {
  ns <- shiny::NS(id)

  bigdash::bigPage(
    id = id,
    navbar = shiny::div(style = "visibility: hidden; display: none;"),
    sidebar = bigdash::sidebar(
      "Methylome Profiler",
      id = id,
      bigdash::sidebarItem("Sample ledger", ns("ledger-tab")),
      bigdash::sidebarItem("Epigenetic age", ns("age-tab")),
      bigdash::sidebarItem("Methylome character", ns("character-tab")),
      bigdash::sidebarItem("Cell composition", ns("deconv-tab")),
      bigdash::sidebarItem("EWAS", ns("ewas-tab"))
    ),
    ## tabSettings need a settings container to render into; without this the
    ## right-hand settings panel never appears. qsee_ui() declares the same.
    settings = bigdash::settings("Settings", id = id),
    bigdash::bigTabs(
      id = id,

      ## ------------------------------------------------------------ ledger --
      bigdash::bigTabItem(
        ns("ledger-tab"),
        mp_tab(bslib::layout_columns(
          col_widths = 12,
          height = MP_TAB_HEIGHT,
          row_heights = list("auto", 1),
          bs_alert("Per-sample quality and identity for a methylation cohort. Nothing on this screen needs a contrast; it answers what is true of each sample."),
          bslib::layout_columns(
            height = "100%",
            col_widths = c(7, 5),
            methylome_table_ledger_ui(
              ns("ledger_tbl"),
              title = "Per-sample checks",
              info.text = "One row per sample. Bimodality is the fraction of probes outside the intermediate 0.3-0.7 band; a healthy methylome sits near 0.9. Imprint drift is the mean absolute deviation of imprinted DMRs from their expected 0.5. DNAm age is the first clock with complete probe coverage.",
              caption = "Per-sample quality, identity and epigenetic age.",
              height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
            ),
            methylome_plot_betadist_ui(
              ns("ledger_dens"),
              title = "Beta distribution",
              caption = "Density of beta values per sample.",
              info.text = "Beta values cover [0,1] and should be bimodal, enriched near 0.2 (hypomethylated) and 0.8 (hypermethylated). A sample that flattens out has usually failed normalisation.",
              info.methods = "Per-sample kernel density over a fixed random subsample of 6,000 probes, drawn with a common seed so samples remain comparable. Dashed guides at beta 0.2 and 0.8.",
              height = c("100%", "700px"), width = c("auto", "100%")
            )
          )
        ))
      ),

      ## --------------------------------------------------------------- age --
      bigdash::bigTabItem(
        ns("age-tab"),
        ## Clock selection lives in the board namespace: it drives all four
        ## panels on this screen.
        methylome_age_inputs(id),
        mp_tab(bslib::layout_columns(
          col_widths = 12,
          height = MP_TAB_HEIGHT,
          row_heights = list("auto", 1.15, 1),
          bs_alert("Epigenetic clocks for every sample. Pick the clocks in the settings panel on the right, by family or one at a time. A clock computed on a fraction of its probes returns a confident wrong number, so anything below the coverage floor is withheld rather than shown."),
          bslib::layout_columns(
            height = "100%",
            col_widths = c(7, 5),
            methylome_plot_agecor_ui(
              ns("age_cor"),
              title = "DNAm age vs chronological age",
              caption = "Each clock against the reported age of the sample.",
              info.text = "The headline check on any epigenetic clock: how closely it tracks the age actually recorded for the sample. Points on the dashed identity line are samples whose methylome matches their years; the fitted line shows the clock's overall calibration in this cohort.",
              info.methods = "Up to four selected clocks, each plotted against the first chronological age column found in the sample sheet. Dashed line is parity, solid line is the least-squares fit, and the panel title carries the Pearson correlation.",
              height = c("100%", "700px"), width = c("auto", "100%")
            ),
            methylome_plot_clocks_ui(
              ns("age_clocks"),
              title = "Clock agreement",
              caption = "Every selected clock against every other.",
              info.text = "Different clocks were trained on different tissues and CpG sets. Where they agree the estimate is trustworthy; where they diverge, the sample or its normalisation deserves a look.",
              info.methods = "Pairwise scatter of the selected clocks, identity line dotted in the lower panels, Pearson correlation printed in the upper panels.",
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
              info.text = "Age acceleration is what is left of the epigenetic age once the chronological age is accounted for: positive means the methylome looks older than the person is. This is the point where a per-sample metric becomes a comparison the rest of the platform can use. Pick the clock and the phenotype in the settings panel.",
              info.methods = "Residuals of a linear fit of the selected clock on chronological age, boxplotted against the selected phenotype, with samples jittered over it, a two-sample t-test for two levels and a Kruskal-Wallis test otherwise. Ticking cell-composition adjustment adds the estimated proportions to that fit, giving intrinsic acceleration - without it, a cohort whose groups differ in blood composition reports that difference as an accelerated methylome. Without a chronological age column the raw DNAm age is shown instead.",
              height = c("100%", "700px"), width = c("auto", "100%")
            ),
            methylome_table_coverage_ui(
              ns("age_cov"),
              title = "Clocks and coverage",
              info.text = "What each selected clock is, its median estimate, and whether this dataset carries enough of its CpGs to compute it at all.",
              caption = "Per-clock coverage and whether the estimate is usable.",
              height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
            )
          )
        ))
      ),
      ## --------------------------------------------------------- character --
      bigdash::bigTabItem(
        ns("character-tab"),
        mp_tab(bslib::layout_columns(
          col_widths = 12,
          height = MP_TAB_HEIGHT,
          row_heights = list("auto", 1),
          bs_alert("Where methylation sits in the genome. Both panels read the CpG island and gene-position annotation that playbase produces for every methylation dataset."),
          bslib::layout_columns(
            height = "100%",
            col_widths = c(6, 6),
            methylome_plot_context_ui(
              ns("char_island"),
              title = "Mean beta by relation to CpG island",
              caption = "Methylation across islands, shores, shelves and open sea.",
              info.text = "CpG islands are normally unmethylated while the surrounding shores, shelves and open sea are progressively more methylated. A flat profile here means something has gone wrong upstream.",
              info.methods = "Per-probe mean beta across all samples, averaged within each Relation_to_Island category from the array annotation. Bars are coloured by the beta value itself.",
              height = c("100%", "700px"), width = c("auto", "100%"), label = "a"
            ),
            methylome_plot_context_ui(
              ns("char_gene"),
              title = "Mean beta by position in gene",
              caption = "Methylation across promoter, exon, body and UTR.",
              info.text = "Promoters (TSS200, TSS1500, 1st exon) are normally unmethylated while gene bodies and 3'UTRs are methylated. This promoter-to-body gradient is the first thing a methylation analyst checks.",
              info.methods = "Per-probe mean beta across all samples, averaged within each UCSC RefGene group from the array annotation. Probes annotated to several groups are assigned their first.",
              height = c("100%", "700px"), width = c("auto", "100%"), label = "b"
            )
          )
        ))
      ),

      ## ------------------------------------------------------ composition --
      bigdash::bigTabItem(
        ns("deconv-tab"),
        methylome_deconv_inputs(id),
        mp_tab(bslib::layout_columns(
          col_widths = 12,
          height = MP_TAB_HEIGHT,
          row_heights = list("auto", 1, 1),
          bs_alert("Estimated cell proportions, projected onto a reference panel. Cell composition is the dominant confounder in bulk-tissue methylation, so this screen exists as much to supply covariates to the EWAS model as to be read on its own."),
          bslib::layout_columns(
            height = "100%",
            col_widths = c(7, 5),
            methylome_plot_composition_ui(
              ns("comp_bars"),
              title = "Composition per sample",
              caption = "Estimated proportion of each reference cell type.",
              info.text = "One stacked bar per sample. Samples are ordered by the first categorical phenotype, so a difference in composition between groups shows up as a visible block rather than being scattered across the axis.",
              info.methods = "Houseman-style constrained projection of each sample's betas onto a reference panel of cell-type-discriminating CpGs, via meffilEstimateCellCountsFromBetas. Proportions are not forced to sum to exactly one; small deviations are expected.",
              height = c("100%", "700px"), width = c("auto", "100%")
            ),
            methylome_table_composition_ui(
              ns("comp_tbl"),
              title = "Estimated proportions",
              info.text = "The per-sample proportions, downloadable as CSV for use elsewhere.",
              caption = "Per-sample cell-type proportions.",
              height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
            )
          ),
          methylome_plot_compgroup_ui(
            ns("comp_group"),
            title = "Composition by phenotype",
            caption = "Each cell type compared across the groups of the contrast.",
            info.text = "This is the confounding check. If composition differs between the groups being compared, an unadjusted differential-methylation result is largely a picture of that difference rather than of the phenotype - adjust for these proportions in the EWAS model.",
            info.methods = "Boxplot of each estimated proportion against the selected phenotype, samples jittered over it, with a two-sample t-test for two-level phenotypes and a Kruskal-Wallis test otherwise.",
            height = c("100%", "700px"), width = c("auto", "100%")
          )
        ))
      ),

      ## -------------------------------------------------------------- EWAS --
      ## Two sub-tabs: the result (where the signal is, and which CpGs) and
      ## the supporting views (is it inflated, where does it sit, and do the
      ## top hits hold up per sample).
      bigdash::bigTabItem(
        ns("ewas-tab"),
        ## Inputs are declared in the board namespace, not a sub-module: the
        ## threshold drives both sub-tabs.
        methylome_ewas_inputs(id),
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
              bs_alert("The differential methylation already computed for this dataset, shown genome-wide. Set the threshold in the settings panel on the right; the line, the coloured points and the table all follow it."),
              methylome_plot_manhattan_ui(
                ns("ewas_manhattan"),
                title = "Manhattan",
                caption = "Genome-wide significance by chromosome position.",
                info.text = "Every tested CpG plotted at its genomic position against significance. Probes passing the threshold are coloured by direction: red where methylation increases, blue where it decreases.",
                info.methods = "Probes are ordered by chromosome and position using the array annotation; the chromosome field is a cytoband so the arm and band are stripped first. Alternating grey shades separate chromosomes. A q-value cut-off has no fixed position on a -log10(p) axis, so the dashed line is drawn at the largest p-value that still passes, and is omitted entirely when nothing does.",
                height = c("100%", "700px"), width = c("auto", "100%")
              ),
              methylome_table_hits_ui(
                ns("ewas_hits"),
                title = "CpGs passing the threshold",
                info.text = "Every CpG called significant at the current cut-off, most significant first. Delta beta is the difference in mean beta between the two groups of the contrast - the change in methylation itself, not the M-value log fold change the model was fitted on.",
                caption = "Significant CpGs with gene, genomic context and effect size.",
                height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
              )
            )
          ),

          bslib::nav_panel(
            title = "Regions & pathways", value = "regions",
            bslib::layout_columns(
              col_widths = 12,
              height = MP_SUBTAB_HEIGHT,
              row_heights = list("auto", 1.15, 1),
              bs_alert("Region-level and gene-set views of the current model. Both run on demand from the settings panel, against the outcome and threshold set on the other sub-tabs."),
              methylome_plot_dmrregion_ui(
                ns("ewas_dmrplot"),
                title = "Region detail",
                caption = "Methylation across one called region, by genomic position.",
                info.text = "A region in the table is only a set of numbers; this is what it looks like. A real DMR is a run of neighbouring CpGs all moving the same way, separating the groups across the whole region. One CpG pulling its neighbours along produces a similar table row and an obviously different picture.",
                info.methods = "Beta values of every CpG in the region plus 30% flanking context on each side, drawn against genomic position. Faint lines are individual samples, heavy lines the mean of each group - or of the low, middle and top tertile when the outcome is continuous. The shaded box is the called region, so flanking probes are visibly outside it, and the track under the axis carries the CpG island context of each probe.",
                height = c("100%", "700px"), width = c("auto", "100%")
              ),
              bslib::layout_columns(
                height = "100%",
                col_widths = c(6, 6),
                methylome_table_dmr_ui(
                  ns("ewas_dmr"),
                  title = "Differentially methylated regions",
                  info.text = "Contiguous runs of CpGs that move together, which is the unit methylation phenotypes actually occur in - promoter island hypermethylation, hypomethylated blocks - rather than isolated probes. Single-CpG regions are excluded; those are already in the hits table.",
                  caption = "Regions called from the fitted model by dmrff.",
                  height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
                ),
                methylome_table_enrich_ui(
                  ns("ewas_enrichgs"),
                  title = "Gene sets (probe-bias corrected)",
                  info.text = "Gene-set enrichment of the hit list, corrected for the number of probes per gene. Genes carry very different probe counts, so an uncorrected test returns the same probe-dense developmental sets whatever the biology; gometh reweights with Wallenius noncentral hypergeometric and also handles CpGs mapping to several genes.",
                  caption = "GO or KEGG terms enriched among the significant CpGs.",
                  height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
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
                  info.text = "The genomic inflation factor lambda summarises how far the p-value distribution departs from the null. It rises with genuine widespread signal as well as with confounding, so it is a prompt to look rather than a pass or fail.",
                  info.methods = "Observed -log10 p-values against the uniform expectation. Lambda is the median chi-square statistic divided by its null expectation.",
                  height = c("100%", "700px"), width = c("auto", "100%")
                ),
                methylome_plot_enrichment_ui(
                  ns("ewas_enrich"),
                  title = "Where the hits sit",
                  caption = "Genomic context of the significant CpGs.",
                  info.text = "Whether the significant CpGs concentrate in islands, shores, shelves or open sea. Shore enrichment with island depletion is the classic pattern in differential methylation.",
                  info.methods = "Odds ratio of each Relation_to_Island category among probes passing the threshold versus the probes actually tested in this contrast, not the whole array. A half-count is added to each cell so empty categories remain finite.",
                  height = c("100%", "700px"), width = c("auto", "100%")
                )
              ),
              methylome_plot_stripcharts_ui(
                ns("ewas_strips"),
                title = "Top CpGs, sample by sample",
                caption = "Per-sample beta for the most significant CpGs, split by group.",
                info.text = "Beta value of every sample at each of the top CpGs. This is the panel that shows whether a hit is real or driven by a couple of outliers, and at what absolute level of methylation it sits - a 0.05 difference at beta 0.50 and at beta 0.02 are very different claims.",
                info.methods = "The most significant CpGs passing the current threshold, one panel each, points jittered within group with a bar at the group mean. The y axis is fixed to the full [0,1] beta range rather than zooming to the data, so panels stay comparable and the absolute level stays visible. Dotted guides at 0.2 and 0.8. How many CpGs are drawn is set in the settings panel.",
                height = c("100%", "700px"), width = c("auto", "100%")
              )
            )
          )
        ))
      )
    )
  )
}
