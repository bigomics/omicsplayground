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
      bigdash::sidebarItem("EWAS", ns("ewas-tab"))
    ),
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
        mp_tab(bslib::layout_columns(
          col_widths = 12,
          height = MP_TAB_HEIGHT,
          row_heights = list("auto", 1, "auto"),
          bs_alert("Epigenetic age from the wateRmelon clocks. A clock computed on a fraction of its probes returns a confident wrong number, so any clock missing probes is withheld rather than shown with a caveat."),
          bslib::layout_columns(
            height = "100%",
            col_widths = c(6, 6),
            methylome_plot_clocks_ui(
              ns("age_clocks"),
              title = "Clock agreement",
              caption = "Two independent clocks plotted against each other.",
              info.text = "Different clocks were trained on different tissues and probe sets. Where they agree, the age estimate is trustworthy; where they diverge, the sample or its normalisation deserves a look.",
              info.methods = "The first two clocks with complete probe coverage are plotted against each other with the identity line and Pearson correlation.",
              height = c("100%", "700px"), width = c("auto", "100%")
            ),
            methylome_plot_agegroup_ui(
              ns("age_group"),
              title = "DNAm age by phenotype",
              caption = "Epigenetic age split by the first two-level phenotype.",
              info.text = "Epigenetic age compared across the groups defined in the sample sheet. This is the point where a per-sample metric becomes a comparison the rest of the platform can use.",
              info.methods = "Boxplot of the first fully covered clock against the first two-level phenotype column, with the individual samples jittered over it.",
              height = c("100%", "700px"), width = c("auto", "100%")
            )
          ),
          methylome_table_coverage_ui(
            ns("age_cov"),
            title = "Per-clock coverage",
            info.text = "How many of each clock's required probes are present in this dataset. Clocks with any missing probes have their age withheld.",
            caption = "Clock probe coverage and whether the estimate is usable.",
            height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
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

      ## -------------------------------------------------------------- EWAS --
      bigdash::bigTabItem(
        ns("ewas-tab"),
        ## Inputs are declared in the board namespace, not a sub-module: the
        ## threshold drives the Manhattan, the enrichment and the hit table.
        methylome_ewas_inputs(id),
        mp_tab(bslib::layout_columns(
          col_widths = 12,
          height = MP_TAB_HEIGHT,
          row_heights = list("auto", 1, 1.45),
          bs_alert("The differential methylation already computed for this dataset, shown genome-wide. This is the one screen here that uses a contrast."),
          methylome_plot_manhattan_ui(
            ns("ewas_manhattan"),
            title = "Manhattan",
            caption = "Genome-wide significance by chromosome position.",
            info.text = "Every tested CpG plotted at its genomic position against significance. Significant probes are coloured by direction: red where methylation increases, blue where it decreases.",
            info.methods = "Probes are ordered by chromosome and position using the array annotation; the chromosome field is a cytoband so the arm and band are stripped first. Alternating grey shades separate chromosomes. The dashed line is the conventional 1e-7 array-wide threshold rather than a naive Bonferroni over all probes.",
            height = c("100%", "700px"), width = c("auto", "100%")
          ),
          bslib::layout_columns(
            height = "100%",
            col_widths = c(7, 5),
            methylome_table_hits_ui(
              ns("ewas_hits"),
              title = "CpGs passing the threshold",
              info.text = "Every CpG called significant at the current cut-off, most significant first. Delta beta is the difference in mean beta between the two groups of the contrast - the change in methylation itself, not the M-value log fold change the model was fitted on.",
              caption = "Significant CpGs with gene, genomic context and effect size.",
              height = c("100%", TABLE_HEIGHT_MODAL), width = c("auto", "100%")
            ),
            bslib::layout_columns(
              height = "100%",
              col_widths = 12,
              row_heights = list(1, 1),
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
              info.methods = "Odds ratio of each Relation_to_Island category among significant probes versus the probes actually tested in this contrast, not the whole array. A half-count is added to each cell so empty categories remain finite.",
              height = c("100%", "700px"), width = c("auto", "100%")
            )
            )
          )
        ))
      )
    )
  )
}
