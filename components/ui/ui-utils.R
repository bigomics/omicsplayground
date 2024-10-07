##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


visPrint <- function(visnet, file, width = 3000, height = 3000, delay = 0, zoom = 1) {
  is.pdf <- grepl("pdf$", file)
  if (is.pdf) {
    width <- width * 600
    height <- height * 600
  }
  vis2 <- htmlwidgets::createWidget(
    name = "visNetwork",
    x = visnet$x,
    width = width, height = height,
    package = "visNetwork"
  )
  tmp.html <- paste0(tempfile(), "-visnet.html")
  tmp.png <- paste0(tempfile(), "-webshot.png")
  visNetwork::visSave(vis2, file = tmp.html)
  webshot2::webshot(
    url = tmp.html,
    file = tmp.png,
    selector = "#htmlwidget_container",
    delay = delay,
    zoom = zoom,
    cliprect = "viewport",
    vwidth = width,
    vheight = height
  )
  if (is.pdf) {
    cmd <- paste("convert", tmp.png, "-density 600", file)
    system(cmd)
  } else {
    file.copy(tmp.png, file, overwrite = TRUE)
  }
  unlink(tmp.html)
}


addWatermark.PDF <- function(file) {
  if (system("which pdftk", ignore.stdout = TRUE) == 1) {
    return
  } ## if no pdftk installed...
  mark <- file.path(FILES, "watermark.pdf")
  tmp <- paste0(gsub("file", "plot", tempfile()), ".pdf")
  cmd <- paste("pdftk", file, "stamp", mark, "output", tmp) ## NEED pdftk installed!!!
  cmd
  system(cmd)
  file.copy(tmp, file, overwrite = TRUE)
  unlink(tmp)
}


addWatermark.PNG2 <- function(file, out = file,
                              mark = file.path(FILES, "watermark-logo.png"),
                              logo.scale = 0.045, position = "topright") {
  if (system("which convert", ignore.stdout = TRUE) == 1) {
    return
  } ## if no pdftk installed...
  if (position %in% c(FALSE, "none")) {
    return
  }
  img <- png::readPNG(file)
  w <- dim(img)[2]
  h <- dim(img)[1]
  tmp <- paste0(gsub("file", "plot", tempfile()), ".png")
  logo.height <- max(logo.scale * h, logo.scale * w / 3)
  logo.width <- logo.height * 4
  logo.x <- 0.2 * logo.width
  logo.y <- 0.3 * logo.height
  if (grepl("right", position)) logo.x <- w - 1.2 * logo.width
  if (grepl("bottom", position)) logo.y <- h - 0.3 * logo.height
  cmd.str <- "convert %s \\( %s -thumbnail %.0fx%.0f \\) -geometry +%.0f+%.0f -composite %s"
  cmd <- sprintf(cmd.str, file, mark, logo.width, logo.height, logo.x, logo.y, tmp)
  system(cmd)
  file.copy(tmp, out, overwrite = TRUE)
  unlink(tmp)
}

addWatermark.PDF2 <- function(file, w, h, out = file,
                              mark = file.path(FILES, "watermark-logo.pdf"),
                              logo.scale = 0.045, position = "topright") {
  if (system("which pdftk", ignore.stdout = TRUE) == 1) {
    return
  } ## if no pdftk installed...
  if (position %in% c(FALSE, "none")) {
    return
  }
  tmp1 <- paste0(gsub("file", "plot", tempfile()), ".pdf")
  tmp2 <- paste0(gsub("file", "plot", tempfile()), ".pdf")
  tmp3 <- paste0(gsub("file", "plot", tempfile()), ".pdf")
  logo.height <- max(logo.scale * h, logo.scale * w / 3) * 720
  logo.width <- logo.height * 4
  logo.x <- 0.2 * logo.width
  logo.y <- h * 720 - 1.66 * logo.height
  if (grepl("right", position)) logo.x <- w * 720 - 1.2 * logo.width
  if (grepl("bottom", position)) logo.y <- 0.66 * logo.height

  scale.cmd <- sprintf("gs -o %s -sDEVICE=pdfwrite  -dDEVICEWIDTH=%0.f -dDEVICEHEIGHT=%0.f -dPDFFitPage -dAutoRotatePages=/None -f %s", tmp1, logo.width, logo.height, mark)
  system(scale.cmd, ignore.stdout = TRUE, ignore.stderr = FALSE)
  ## create empty page with same plot size and logo translated to topleft
  cmd1 <- sprintf(
    "gs -o %s -sDEVICE=pdfwrite -g%.0fx%.0f -c '<</PageOffset [%.0f %.0f]>> setpagedevice' -f %s",
    tmp2, w * 720, h * 720, logo.x / 10, logo.y / 10, tmp1
  )
  system(cmd1, ignore.stdout = TRUE, ignore.stderr = FALSE)

  ## merge: overlay logo and plot
  cmd3 <- paste("pdftk", file, "stamp", tmp2, "output", tmp3) ## NEED pdftk installed!!!
  system(cmd3)
  file.copy(tmp3, out, overwrite = TRUE)
  unlink(tmp1)
  unlink(tmp2)
  unlink(tmp3)
}

gadgetize <- function(moduleUI, moduleSERVER, title = "shiny gadget", ...) {
  ## Creates gadget from a Shiny module. Gadget are browser-based UI
  ## applets for a single task that can be run from the R command
  ## line. It is not to be used inside Shiny programs themselves.
  ##

  id <- sub(".*file", "gadget", tempfile()) ## random ID
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar(title),
    miniUI::miniContentPanel(moduleUI(id))
  )
  server <- function(input, output, session) {
    return_obj <- moduleSERVER(id, ...)
    shiny::observeEvent(input$done, {
      shiny::stopApp(return_obj())
    })
  }
  X <- shiny::runGadget(ui, server)
  cat("[gadgetize] names(X)=", names(X), "\n")
  cat("[gadgetize] *** closing gadget ***\n")
  X
}

gadgetize2 <- function(moduleUI, moduleSERVER, title = "shiny gadget",
                       size = "m", ...) {
  ##
  ## Creates modalDialog from a Shiny module similar as used inside
  ## Shiny programs.
  ##


  id <- sub(".*file", "gadget", tempfile()) ## random ID
  ui <- shiny::fluidPage()
  server <- function(input, output, session) {
    return_obj <- moduleSERVER(id, ...)
    shiny::showModal(shiny::modalDialog(
      moduleUI(id),
      footer = shiny::tagList(
        shiny::actionButton("gdgt_close", "X")
      ),
      size = size,
      easyClose = FALSE,
      fade = FALSE
    ))
    shiny::observeEvent(input$gdgt_close, {
      cat("[gadgetize2] *** closing gadget ***\n")
      shiny::stopApp(return_obj())
    })
  }

  pgx <- shiny::runGadget(ui, server)
  cat(names(pgx))
  pgx
}

this.style <- function(id, css, ns = NULL) {
  if (!is.null(ns)) id <- ns(id)
  shiny::tags$head(shiny::tags$style(paste0("#", id, " ", css)))
}

alertDataLoaded <- function(session, ngs) {
  if (!is.null(ngs)) {
    return()
  }
  message("[alertDataLoaded] WARNING:: no PGX object")
}

pgx.randomCartoon <- function() {
  cartoon_list <- list(
    list(slogan = "Visual analytics. See and understand", img = "data-graph-wisdom.jpg"),
    list(slogan = "Fasten your seat belts. Accelerated discovery", img = "cartoon-speedup.jpg"),
    list(slogan = "Analytics anywhere. Anytime.", img = "cartoon-cloudservice.jpg"),
    list(slogan = "Analyze with confidence. Be a rockstar", img = "bigomics-rockstar3.jpg"),
    list(slogan = "Fast track your Bioinformatics", img = "selfservice-checkout2.png"),
    list(slogan = "Integrate more. Dig deeper", img = "cartoon-integration.jpg"),
    list(slogan = "Your analysis doesn't take coffee breaks", img = "gone-for-coffee.png"),
    list(slogan = "Too much data? Help yourself", img = "cartoon-datahelp2.jpg"),
    list(slogan = "Big Friendly Omics", img = "big-friendly-omics1.jpg"),
    list(slogan = "Big Data meets Biology", img = "bigdata-meets.png")
  )

  cartoon <- sample(cartoon_list, 1)[[1]]
  cartoon$img2 <- file.path("cartoons", cartoon$img)
  cartoon$img <- file.path("www/cartoons", cartoon$img)
  cartoon
}

pgx.showCartoonModal <- function(msg = "Loading data...", img.path = "www/cartoons") {
  cartoon_list <- list(
    list(slogan = "Visual analytics. See and understand", img = "data-graph-wisdom.jpg"),
    list(slogan = "Fasten your seat belts. Accelerated discovery", img = "cartoon-speedup.jpg"),
    list(slogan = "Analytics anywhere. Anytime.", img = "cartoon-cloudservice.jpg"),
    list(slogan = "Analyze with confidence. Be a rockstar", img = "bigomics-rockstar3.jpg"),
    list(slogan = "Fast track your Bioinformatics", img = "selfservice-checkout2.png"),
    list(slogan = "Integrate more. Dig deeper", img = "cartoon-integration.jpg"),
    list(slogan = "Your analysis doesn't take coffee breaks", img = "gone-for-coffee.png"),
    list(slogan = "Too much data? Help yourself", img = "cartoon-datahelp2.jpg"),
    list(slogan = "Big Friendly Omics", img = "big-friendly-omics1.jpg"),
    list(slogan = "Big Data meets Biology", img = "bigdata-meets.png")
  )

  randomCartoon <- function() {
    cartoon <- sample(cartoon_list, 1)[[1]]
    cartoon$img2 <- paste0("static/cartoons/", cartoon$img)
    cartoon$img <- file.path(img.path, cartoon$img)
    cartoon
  }

  toon <- randomCartoon()
  shiny::showModal(shiny::modalDialog(
    title = shiny::div(shiny::h2(toon$slogan), shiny::p("with Omics Playground"), style = "text-align:center;"),
    shiny::img(src = toon$img2, class = "img-fluid"),
    footer = fillRow(flex = c(1, NA, 1), " ", msg, " "),
    size = "l",
    easyClose = FALSE,
    fade = TRUE
  ))
}

HandleNoLinkFound <- function(wrapHyperLinkOutput, NoLinkString, SubstituteString) {
  pattern <- paste0("^", NoLinkString, "$")
  special_cases <- grepl(pattern, wrapHyperLinkOutput, perl = TRUE)
  wrapHyperLinkOutput[special_cases] <- SubstituteString
  return(wrapHyperLinkOutput)
}

getSettings <- function(ns, session) {
  # Get board/plot ns
  board_ns <- sub("-.*", "", ns(""))
  plot_ns <- sub(".*-(.*?)-.*", "\\1", ns(""))
  # Get board inputs
  board_inputs <- names(.subset2(session, "parent")$input)[grepl(board_ns, names(.subset2(session, "parent")$input))]
  board_inputs <- board_inputs[substr(board_inputs, 1, nchar(board_ns) + 1) == paste0(board_ns, "-")]
  # Get board settings
  board_settings <- board_inputs[grep("^[^-]*-[^-]*$", board_inputs)]
  # Remove `data_options`, `tabs` `board_info`
  board_settings <- board_settings[!grepl("data_options|tabs|info|options|compute|pdx_runbutton", board_settings)]
  # Get settings values
  board_settings_values <- lapply(board_settings, function(x) {
    val <- .subset2(session, "parent")$input[[x]]
    if (is.null(val)) val <- ""
    if (any(nchar(val) > 30)) val <- paste0(substr(val, 1, 30), "...")
    val <- paste(val, collapse = ", ")
    return(val)
  }) |> unlist()
  # Merge values and input names (without namespacing)
  settings_table <- data.frame(
    setting = sub("^[^-]*-", "", board_settings),
    value = board_settings_values
  )
  # Correct column names (of board settings)
  df_names <- lapply(settings_table[,1], function(x) {
    tspan(inputLabelDictionary(board_ns, x), js = FALSE)
  }) |> unlist()
  settings_table[,1] <- df_names

  # Get plot inputs
  plot_inputs <- board_inputs[grepl(plot_ns, board_inputs)]
  # Get plot settings
  plot_settings <- plot_inputs[grep("^[^-]*-[^-]*-[^-]*$", plot_inputs)]
  # Get plot values
  plot_settings_values <- lapply(plot_settings, function(x) {
    val <- .subset2(session, "parent")$input[[x]]
    if (is.null(val)) val <- ""
    if (any(nchar(val) > 30)) val <- paste0(substr(val, 1, 30), "...")
    val <- paste(val, collapse = ", ")
    return(val)
  }) |> unlist()
  # Merge values and input names (without namespacing)
  plot_table <- data.frame(
    setting = sub("^[^-]*-", "", sub("^[^-]*-[^-]*-", "", plot_settings)),
    value = plot_settings_values
  )

  # Get loaded metadata
  timestamp <- as.character(format(Sys.time(), "%a %b %d %X %Y"))
  version <- scan(file.path(OPG, "VERSION"), character())[1]
  dataset <- LOADEDPGX
  datatype <- DATATYPEPGX
  metadata <- data.frame(
    setting = c("Dataset", "Timestamp", "Data type", "Version"),
    value = c(dataset, timestamp, datatype, version)
  )

  ## setting as string
  df <- rbind(metadata, plot_table, settings_table)
  settings_str <- paste(paste0(df[,1],"=",df[,2]),collapse=";")
  
  list(
    metadata = metadata,
    plot_table = plot_table,
    settings_table = settings_table,
    settings_str = settings_str)
}


addSettings <- function(ns, session, file) {

  ##board_ns <- sub("-.*", "", ns(""))
  settings <- getSettings(ns, session)

  # Merge plot and settings
  df <- rbind(
    settings$metadata, c("", ""),
    c("Plot option", "Value"), settings$plot_table, c("", ""),
    c("Setting", "Value"), settings$settings_table
  )

  ## # Correct column names
  ## df_names <- lapply(df$setting, function(x) {
  ##   tspan(inputLabelDictionary(board_ns, x), js = FALSE)
  ## }) |> unlist()
  ## df$setting <- df_names
  
  # Setup table theme
  table_theme <- gridExtra::ttheme_minimal(
    colhead = list(
      fg_params = list(
        fontface = "bold", # Bold font for headers
        hjust = 0, # Left-align the text
        x = 0 # Align text to the left within the cell
      )
    ),
    core = list(
      fg_params = list(
        fontface = c(
          rep("plain", nrow(settings$metadata) + 1),
          "bold",
          rep("plain", nrow(settings$plot_table) + 1),
          "bold",
          rep("plain", nrow(settings$settings_table))
        ),
        hjust = 0,
        x = 0
      )
    )
  )

  # Compute PDF height using nrow
  height <- nrow(df) * 0.4

  # Print PDF temp table
  df_pdf <- tempfile(fileext = ".pdf")
  final_pdf <- tempfile(fileext = ".pdf")
  pdf(df_pdf, height = height, width = 10)
  gridExtra::grid.table(df, rows = NULL, col = c("Metadata", "Value"), theme = table_theme)
  dev.off()
  # Construct the pdftk command
  pdftk_command <- sprintf("pdftk %s %s cat output %s", file, df_pdf, final_pdf)
  # Execute the command
  system(pdftk_command)
  ## finally copy to final exported file
  dbg("[downloadHandler.PDF] copy PDFFILE", final_pdf, "to download file", file)
  file.copy(final_pdf, file, overwrite = TRUE)
}

addSettings.SAVE <- function(ns, session, file) {
  # Get board/plot ns
  board_ns <- sub("-.*", "", ns(""))
  plot_ns <- sub(".*-(.*?)-.*", "\\1", ns(""))
  # Get board inputs
  board_inputs <- names(.subset2(session, "parent")$input)[grepl(board_ns, names(.subset2(session, "parent")$input))]
  board_inputs <- board_inputs[substr(board_inputs, 1, nchar(board_ns) + 1) == paste0(board_ns, "-")]
  # Get board settings
  board_settings <- board_inputs[grep("^[^-]*-[^-]*$", board_inputs)]
  # Remove `data_options`, `tabs` `board_info`
  board_settings <- board_settings[!grepl("data_options|tabs|info|options|compute|pdx_runbutton", board_settings)]
  # Get settings values
  board_settings_values <- lapply(board_settings, function(x) {
    val <- .subset2(session, "parent")$input[[x]]
    if (is.null(val)) val <- ""
    if (any(nchar(val) > 30)) val <- paste0(substr(val, 1, 30), "...")
    val <- paste(val, collapse = ", ")
    return(val)
  }) |> unlist()
  # Merge values and input names (without namespacing)
  settings_table <- data.frame(
    setting = sub("^[^-]*-", "", board_settings),
    value = board_settings_values
  )

  # Get plot inputs
  plot_inputs <- board_inputs[grepl(plot_ns, board_inputs)]
  # Get plot settings
  plot_settings <- plot_inputs[grep("^[^-]*-[^-]*-[^-]*$", plot_inputs)]
  # Get plot values
  plot_settings_values <- lapply(plot_settings, function(x) {
    val <- .subset2(session, "parent")$input[[x]]
    if (is.null(val)) val <- ""
    if (any(nchar(val) > 30)) val <- paste0(substr(val, 1, 30), "...")
    val <- paste(val, collapse = ", ")
    return(val)
  }) |> unlist()
  # Merge values and input names (without namespacing)
  plot_table <- data.frame(
    setting = sub("^[^-]*-", "", sub("^[^-]*-[^-]*-", "", plot_settings)),
    value = plot_settings_values
  )

  # Get loaded metadata
  timestamp <- as.character(format(Sys.time(), "%a %b %d %X %Y"))
  version <- scan(file.path(OPG, "VERSION"), character())[1]
  dataset <- LOADEDPGX
  datatype <- DATATYPEPGX
  metadata <- data.frame(
    setting = c("Dataset", "Timestamp", "Data type", "Version"),
    value = c(dataset, timestamp, datatype, version)
  )

  # Merge plot and settings
  df <- rbind(
    metadata, c("", ""),
    c("Plot option", "Value"), plot_table, c("", ""),
    c("Setting", "Value"), settings_table
  )

  # Correct column names
  df_names <- lapply(df$setting, function(x) {
    tspan(inputLabelDictionary(board_ns, x), js = FALSE)
  }) |> unlist()
  df$setting <- df_names

  # Setup table theme
  table_theme <- gridExtra::ttheme_minimal(
    colhead = list(
      fg_params = list(
        fontface = "bold", # Bold font for headers
        hjust = 0, # Left-align the text
        x = 0 # Align text to the left within the cell
      )
    ),
    core = list(
      fg_params = list(
        fontface = c(
          rep("plain", nrow(metadata) + 1),
          "bold",
          rep("plain", nrow(plot_table) + 1),
          "bold",
          rep("plain", nrow(settings_table))
        ),
        hjust = 0,
        x = 0
      )
    )
  )

  # Compute PDF height using nrow
  height <- nrow(df) * 0.4

  # Print PDF temp table
  df_pdf <- tempfile(fileext = ".pdf")
  final_pdf <- tempfile(fileext = ".pdf")
  pdf(df_pdf, height = height, width = 10)
  gridExtra::grid.table(df, rows = NULL, col = c("Metadata", "Value"), theme = table_theme)
  dev.off()
  # Construct the pdftk command
  pdftk_command <- sprintf("pdftk %s %s cat output %s", file, df_pdf, final_pdf)
  # Execute the command
  system(pdftk_command)
  ## finally copy to final exported file
  dbg("[downloadHandler.PDF] copy PDFFILE", final_pdf, "to download file", file)
  file.copy(final_pdf, file, overwrite = TRUE)
}

inputLabelDictionary <- function(board_ns, inputId) {
  dictionary <- list(
    dataview = list(
      search_gene = "Gene",
      data_samplefilter = "Filter samples",
      data_groupby = "Group by",
      data_type = "Data type",
      clustsamples = "cluster samples",
      vars = "show variables"
    ),
    clustersamples = list(
      selected_phenotypes = "Show phenotypes",
      hm_splitby = "Split samples by",
      hm_splitvar = "Split samples by (phenotype/gene)",
      hm_average_group = "Split samples by (phenotype/gene) average by group",
      hm_samplefilter = "Filter samples",
      hm_features = "Gene family",
      hm_customfeatures = "Gene family (<custom>)",
      hm_contrast = "Gene family (<contrast>)",
      hm_topmode = "Top mode",
      hm_ntop = "Top N",
      hm_clustk = "K modules",
      hm_scale = "Scale",
      hm_legend = "Show legend",
      hm_cexRow = "cexRow",
      hm_cexCol = "cexCol",
      xann_level = "Reference level",
      xann_odds_weighting = "Reference level ('geneset') Fisher test weighting",
      xann_refset = "Reference level ('geneset') Reference set",
      pca_label = "Label",
      all_clustmethods = "Show all methods",
      plot3d = "Plot 3D",
      showlabels = "Shoe group labels",
      hm_pcaverage = "Average by gene module",
      hm_pcscale = "Scale values",
      gx_grouped = "Group samples",
      hm_filterXY = "Exclude X/Y genes",
      hm_filterMitoRibo = "Exclude mito/ribo genes",
      hmpca.colvar = "Color/label",
      hmpca.shapevar = "Shape",
      hm_clustmethod = "Layout",
      hm_level = "Level"
    ),
    diffexpr = list(
      gx_contrast = "Contrast",
      gx_features = "Gene family",
      gx_fdr = "FDR",
      gx_lfc = "logFC",
      gx_statmethod = "Statistical methods",
      gx_showall = "Show all genes",
      color_up_down = "Color up/down regulated",
      barplot_grouped = "Grouped",
      barplot_logscale = "Log scale",
      barplot_showothers = "Show others",
      gx_logscale = "Log scale",
      gx_grouped = "Group samples",
      gx_showothers = "Show others",
      scale_per_plot = "Scale plots"
    ),
    drug = list(
      dsea_contrast = "Contrast",
      dsea_method = "Analysis type",
      dseatable_filter = "Only annotated drugs",
      dsea_moatype = "Plot type",
      qweight = "q-weighting",
      dsea_normalize = "normalize activation matrix"
    ),
    comp = list(
      contrast1 = "Dataset1",
      dataset2 = "Dataset2: (name)",
      contrast2 = "Dataset2: (contrast)",
      plottype = "Plot type",
      hilighttype = "Highlight genes",
      ntop = "ntop",
      genelist = "Highlight genes (cusom)",
      colorby = "Color by"
    ),
    bio = list(
      pdx_predicted = "Predicted target",
      pdx_filter = "Feature set",
      pdx_samplefilter = "Filter samples",
      pdx_select = "Feature set: <custom> Custom features",
      clust_featureRank_method = "Method"
    ),
    wordcloud = list(
      wc_contrast = "Contrast",
      wordcloud_exclude = "Exclude words",
      wordcloud_colors = "Colors",
      tsne_algo = "Clustering algorithm"
    )
  )
  val <- dictionary[[board_ns]][[inputId]]
  if (is.null(val)) val <- inputId
  return(val)
}


tspan <- function(text, js = TRUE) {
  if (is.null(text)) {
    return(NULL)
  }
  if (length(text) == 0) {
    return(NULL)
  }
  text <- paste(text, collapse = " ")
  if (nchar(text) == 0) {
    return("")
  }
  if (!grepl("gene|counts|transcriptomics|rna-seq|logcpm",
    text,
    ignore.case = TRUE
  )) {
    return(text)
  }
  keys <- c(
    "gene", "Gene", "GENE", "counts", "Counts", "COUNTS",
    "transcriptomics", "Transcriptomics", "RNA-seq",
    "logCPM", "log2p1", "expression", "Expression"
  )
  if (js) {
    i18n.tr <- function(key) shiny::span(class = "i18n", `data-key` = key, key)
  } else {
    i18n.tr <- function(key) i18n$t(key)
  }
  for (k in keys) {
    tt <- i18n.tr(k)
    if (grepl(k, text)) text <- gsub(k, tt, text, ignore.case = FALSE)
  }
  if (js) {
    text <- paste0("<span>", text, "</span>")
    text <- shiny::HTML(text)
  }
  text
}

## forced JS version
jspan <- function(text) tspan(text, js = TRUE)


tspan.SAVE <- function(label) {
  shiny::span(class = "i18n", `data-key` = label, label)
}
