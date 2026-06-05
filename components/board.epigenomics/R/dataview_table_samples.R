## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.

dataview_table_beta_ui <- function(id,
                                   width,
                                   height,
                                   title,
                                   info.text,
                                   caption) {

  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    width = width,
    height = height,
    title = title,
    info.text = info.text,
    caption = caption,
    label = "c"
  )
}

dataview_table_beta_server <- function(id,
                                       pgx,
                                       r.samples = reactive(""),
                                       r.pheno = reactive(""),
                                       scrollY) {

  moduleServer(id, function(input, output, session) {

    table_data <- shiny::reactive({

      shiny::req(pgx$X, pgx$genes, pgx$samples)
      X <- playbase::mToBeta(pgx$X)
      Y <- pgx$samples
      annot <- pgx$genes
      rownames(X) <- sub("_.*", "", rownames(X))
      rownames(annot) <- sub("_.*", "", rownames(annot))
      kk <- intersect(rownames(X), rownames(annot))
      if (length(kk) == 0) return(NULL)
      samples <- r.samples()
      if (!all(samples %in% colnames(X))) return(NULL)
      X <- X[kk, samples, drop = FALSE]
      annot <- annot[kk, , drop = FALSE]
      Y <- Y[samples, , drop = FALSE]

      kk <- grep("chr|chromosome|chrom", tolower(colnames(annot)))[1]
      if (is.na(kk)) return(NULL)
      chr_vals <- as.character(annot[, kk])
      keep <- !is.na(chr_vals) & chr_vals != "" & !grepl("cen", tolower(chr_vals))
      annot <- annot[keep, , drop = FALSE]
      X <- X[keep, , drop = FALSE]
      annot$chrom_tmp <- paste0("chr", sub("^chr", "", sub("[pq].*", "", as.character(annot[, kk]))))

      chroms <- unique(sub("^chr", "", annot$chrom_tmp))
      sex_chr <- intersect(c("X", "Y"), chroms)
      autosomes <- suppressWarnings(sort(as.numeric(setdiff(chroms, sex_chr))))
      chroms <- paste0("chr", c(autosomes[!is.na(autosomes)], sex_chr))

      chrom_colmeans <- function(M) {
        do.call(rbind, lapply(chroms, function(chr) {
          jj <- which(annot$chrom_tmp == chr)
          if (length(jj) == 0) return(matrix(rep(NA, ncol(M)), nrow = 1))
          matrix(round(colMeans(M[jj, , drop = FALSE], na.rm = TRUE), 3), nrow = 1)
        }))
      }

      pheno <- r.pheno()

      if (is.null(pheno) || pheno %in% c("", "<ungrouped>")) {
        dt <- as.data.frame(chrom_colmeans(X))
        colnames(dt) <- paste0("Ave.", colnames(X))
        dt <- cbind(Ave = round(rowMeans(as.matrix(dt), na.rm = TRUE), 3), dt)
      } else {
        pheno <- intersect(pheno, colnames(Y))[1]
        if (is.na(pheno)) return(NULL)
        pheno <- as.character(Y[, pheno])
        kk <- which(!is.na(pheno) & pheno != "")
        pheno <- pheno[kk]
        X <- X[, kk, drop = FALSE]
        groups <- sort(unique(pheno))
        M <- do.call(cbind, lapply(groups, function(g)
          rowMeans(X[, pheno == g, drop = FALSE], na.rm = TRUE)))
        colnames(M) <- groups
        dt <- as.data.frame(chrom_colmeans(M))
        colnames(dt) <- paste0("Ave.", groups)
        dt <- cbind(Ave = round(rowMeans(as.matrix(dt), na.rm = TRUE), 3), dt)
      }

      rownames(dt) <- chroms

      return(dt)

    })

    table.RENDER <- function() {
      dt <- table_data()
      req(dt)
      DT::datatable(dt,
        class = "compact hover",
        rownames = TRUE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = 1),
        options = list(
          dom = "lfrtip",
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = scrollY,
          scrollResize = TRUE,
          deferRender = TRUE
        )
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    modal_table.RENDER <- function() {
      dt <- table_data()
      req(dt)
      DT::datatable(dt,
        class = "compact hover",
        rownames = TRUE,
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        options = list(
          dom = "lfrtip",
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = SCROLLY_MODAL,
          deferRender = TRUE
        )
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "20px", lineHeight = "70%")
    }

    TableModuleServer(
      "datasets",
      func = table.RENDER,
      func2 = modal_table.RENDER,
      selector = "none"
    )
  })

}
