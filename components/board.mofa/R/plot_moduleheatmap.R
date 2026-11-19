##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## Maximum characters drawn for the feature labels of the multi-omics heatmap.
## Gene names can be long semicolon-joined symbol lists (e.g.
## "H3C1;H3C4;H3C3;H3C6;H3C11;H3C8;H3C12;..."), which take over the whole plot.
## Longer labels are shortened with an ellipsis. The default also drives the
## shortened plot; the expanded plot uses the larger value.
MODULEHEATMAP_MAXCHAR <- c(30, 80)

mofa_plot_moduleheatmap_ui <- function(
  id,
  title = "",
  info.text = "",
  caption = "",
  label = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("split"),
      label = "Split heatmap by data types",
      value = TRUE
    ),
    hr(),
    withTooltip(
      shiny::sliderInput(
        inputId = ns("maxchar"),
        label = "Max label length (chars):",
        min = 10, max = 100,
        value = MODULEHEATMAP_MAXCHAR[1],
        step = 5
      ),
      "Longer feature names are shortened with an ellipsis. Affects the downloaded plot too.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

## Shorten long feature labels (e.g. semicolon-joined gene symbol lists) with an
## ellipsis. Feature IDs are mapped onto the gene names that get drawn, those
## names are shortened, and the data-type prefix is put back so the heatmap can
## still split rows by data type. Row names that duplicate after shortening are
## extended until unique, as heatmaps reject duplicate row names.
moduleheatmap_shorten_labels <- function(x, maxchar = 30, gene_table = NULL) {
  if (is.null(x) || !is.numeric(maxchar) || length(maxchar) != 1 ||
    is.na(maxchar)) {
    return(x)
  }
  ## Keep the data-type prefix ("protein:") and the 3-character ellipsis, so
  ## labels stay identifiable and split correctly even at silly limits.
  maxchar <- max(10L, as.integer(maxchar))
  ids <- rownames(x)
  if (is.null(ids)) {
    return(x)
  }

  labels <- playbase::probe2symbol(
    ids, gene_table, "gene_name",
    fill_na = TRUE, add_datatype = TRUE
  )
  if (is.null(labels) || length(labels) != length(ids)) {
    ## no usable gene table: shorten the feature IDs themselves
    labels <- ids
  }
  labels <- as.character(labels)
  labels[is.na(labels)] <- ""

  ## prefix carries the data type, everything after it is the gene name. The
  ## limit applies to the whole drawn label, prefix included: counting the gene
  ## name alone would let "mir:hsa-mir-106a" through at maxchar = 10.
  prefix <- ifelse(
    grepl("^[A-Za-z0-9]+:", labels),
    sub("^([A-Za-z0-9]+:).*", "\\1", labels),
    ""
  )
  gname <- sub("^[A-Za-z0-9]+:", "", labels)
  budget <- as.integer(pmax(4L, maxchar - nchar(prefix)))

  ## shorten each label to the budget left after its prefix
  trunc_to <- function(txt, width) {
    out <- txt
    for (k in unique(width)) {
      jj <- which(width == k)
      out[jj] <- stringr::str_trunc(txt[jj], k)
    }
    out
  }
  gname <- trunc_to(gname, budget)

  ## Names that only differ past the limit shorten to the same text. Number
  ## them, taking the tag out of the same budget so labels stay within it.
  lab <- paste0(prefix, gname)
  if (anyDuplicated(lab)) {
    dup <- duplicated(lab)
    tagged <- lab
    tagged[dup] <- paste0(lab[dup], " (", seq_len(sum(dup)), ")")
    lab[dup] <- trunc_to(tagged[dup], budget[dup])
    ## If the tag itself was cut away, shorten the gene part and append the
    ## number, keeping the data-type prefix intact.
    if (anyDuplicated(lab)) {
      num <- seq_len(sum(dup))
      base <- trunc_to(gname[dup], pmax(4L, budget[dup] - nchar(paste0(" ", num))))
      lab[dup] <- paste0(prefix[dup], base, " ", num)
    }
  }
  ## Last resort: identical even after numbering. Never touch the prefix, so
  ## the data-type split keeps working; widen only the gene part.
  guard <- 0L
  while (anyDuplicated(lab) && guard < 10L) {
    guard <- guard + 1L
    dup <- duplicated(lab)
    lab[dup] <- paste0(prefix[dup], trunc_to(gname[dup], budget[dup] + guard), " ", guard)
  }
  rownames(x) <- lab
  x
}

mofa_plot_moduleheatmap_server <- function(id,
                                           mofa,
                                           pgx,
                                           input_factor = reactive(1),
                                           show_types = reactive(NULL),
                                           ntop = c(40, 40),
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    render_plot <- function(legend, maxchar, n, show_colnames) {
      res <- mofa()
      k <- input_factor()
      factors <- colnames(res$F)
      if (!is.null(k)) shiny::req(k %in% factors)

      show_types <- show_types()
      dtypes <- names(res$ww)
      if (is.null(show_types)) show_types <- dtypes
      shiny::validate(need(
        length(show_types) > 0,
        "Please select at least one datatype"
      ))

      ## playbase::mofa.plot_heatmap truncates the short feature IDs *before*
      ## mapping them onto the long, semicolon-joined symbol lists, so long
      ## gene names reached the plot unshortened (issue #1642). Relabel the
      ## features here, where the shortening takes effect: mofa$X supplies the
      ## drawn labels, and mofa$ww supplies the row names that pick the top
      ## features, so both must be relabelled identically.
      old_x <- rownames(res$X)
      res$X <- moduleheatmap_shorten_labels(
        res$X,
        maxchar = maxchar, gene_table = pgx$genes
      )
      new_x <- stats::setNames(rownames(res$X), old_x)
      ww <- lapply(names(res$ww), function(d) {
        w <- res$ww[[d]]
        lab <- unname(new_x[paste0(d, ":", rownames(w))])
        if (!anyNA(lab)) rownames(w) <- sub("^[A-Za-z0-9]+:", "", lab)
        w
      })
      res$ww <- stats::setNames(ww, names(res$ww))

      playbase::mofa.plot_heatmap(
        res,
        k = k, ## main=k,
        gene_table = NULL,
        ntop = n,
        split = input$split,
        type = "splitmap",
        annot = "pheno",
        maxchar = maxchar,
        show_types = show_types,
        show_legend = legend,
        show_colnames = show_colnames,
        mar = c(3, 0, 0, 0),
        annot.ht = 0.9,
        cexRow = 0.9
      )
    }

    user_maxchar <- shiny::reactive({
      if (is.null(input$maxchar)) {
        MODULEHEATMAP_MAXCHAR[1]
      } else {
        ## floor matches moduleheatmap_shorten_labels(), which needs room for
        ## the ellipsis itself
        max(10L, as.integer(input$maxchar))
      }
    })

    plot.RENDER <- function() {
      render_plot(
        legend = FALSE, maxchar = user_maxchar(), n = ntop[1],
        show_colnames = FALSE
      )
    }

    plot.RENDER2 <- function() {
      render_plot(
        legend = TRUE,
        maxchar = max(MODULEHEATMAP_MAXCHAR[2], user_maxchar()), n = ntop[2],
        show_colnames = TRUE
      )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,
      pdf.width = 8,
      pdf.height = 12,
      res = c(80, 100),
      add.watermark = watermark
    )
  })
}
