##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

DeepNetInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),

    ## data set parameters
    shiny::selectInput(ns("selected_pheno"), "Select phenotype:", choices = NULL, multiple = FALSE),
    shiny::selectInput(ns("show_conditions"), "Show conditions:", choices = NULL, multiple = TRUE),
    shiny::selectInput(ns("show_datatypes"), "Show datatypes:", choices = NULL, multiple = TRUE),    

    p("Network learning:"),
    shiny::actionButton(ns("step"), "Step1", size="xs"),
    shiny::actionButton(ns("step20"), "Step20", size="xs"),
    shiny::actionButton(ns("step100"), "Step100", size="xs"),    
    shiny::actionButton(ns("reset"), "Reset", size="xs"),        
    shiny::br(),
    shiny::br(),    
    shiny::actionLink(ns("options"), "Network options", icon = icon("cog", lib = "glyphicon")),
    shiny::br(), 
    shiny::conditionalPanel(
      "input.options % 2 == 1",
      ns = ns,
      shiny::tagList(
        shiny::radioButtons(ns("model"), "Model:", c("AE","SAE","MLP"), selected="SAE", inline=TRUE),
        shiny::selectInput(ns("layers"), "Layers:", choices = c("mini","medium","deep"),
                           selected="mini"),
        shiny::sliderInput(ns("latent_dim"),"Latent dimension:",4,80,16,8),
        shiny::checkboxInput(ns("augment"), "augment data (10x)", FALSE),
        # shiny::checkboxInput(ns("scaleinput"), "scale input", TRUE),
        # shiny::checkboxInput(ns("sdweight"), "gradient SD weight", TRUE),
        # shiny::checkboxInput(ns("useBN"), "use BatchNorm", TRUE),        
        #shiny::checkboxInput(ns("dropout"), "use dropout", FALSE),
        shiny::checkboxInput(ns("addnoise"), "add internal noise", TRUE),
        shiny::checkboxInput(ns("useGLU"), "use GLU", FALSE),
        shiny::checkboxInput(ns("addgsets"), "add genesets", FALSE)
        #shiny::checkboxInput(ns("multitarget"), "multi target", FALSE)        
        #shiny::selectInput(ns("optim"), "Optimizer",
        #  choices = c("adam","adamw","sgd","lbfgs"), selected="adam"),
        #shiny::selectInput(ns("actfun"), "Activation function",
        #  choices = c("relu","gelu","leaky","silu"), selected="leaky"),
        #shiny::sliderInput(ns("l1regularization"),"L1 regularization (log10):",-4,4,0),
        #shiny::sliderInput(ns("l2regularization"),"L2 regularization (log10):",-4,4,0)
      )
    )
  )
}

DeepNetUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height
  
  shiny::div(
    boardHeader(title = "Multi-Omics Supervised Auto-Encoder", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ##----------------------------------------------------------------
      shiny::tabPanel(
        "Model training",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>Multi-Omics integration by deep learning network</b>. Here we integrate multi-omics data using a multi-view supervised auto-encoder. This architecture is a combination of (multiple) auto-encoders and a multilayer perceptron (MLP) classifier that uses a merge of the bottleneck layers for prediction.")),
          bslib::layout_columns(
            col_widths = bslib::breakpoints(
              xxxl = c(6, 6),
              xl = c(12, 12),              
              sm = c(12, 12)
            ),
            bslib::layout_columns(
              col_widths = bslib::breakpoints(
                xxxl = c(12, 6, 6),
                xl = c(6, 3, 3),              
                sm = c(12, 12, 12)
              ),
              plot_deepnet_diagram_ui(
                ns("deepnet_diagram"),
                title = "Network model architecture",
                info.text = "Supervised Auto-encoder",
                caption = ""
              ),
              plot_deepnet_clusters_ui(
                ns("deepnet_clusters"),
                title = "Network clustering",
                info.text = "Neural net clustering",
                caption = ""
              ),
              plot_deepnet_aescatter_ui(
                ns("deepnet_aescatter"),
                title = "Signal reconstruction",
                info.text = "Neural net clustering",
                caption = ""
              )
            ),
            bslib::layout_columns(
              col_widths = bslib::breakpoints(
                xl = c(6, 3, 3),              
                sm = c(12, 12, 12, 12)
              ),
              plot_deepnet_biomarkerheatmap_ui(
                ns("deepnet_biomarkerheatmap"),
                title = "Biomarker heatmap",
                info.text = "",
                caption = ""
              ),              
              plot_deepnet_confusionmatrix_ui(
                ns("deepnet_confusionmatrix"),
                title = "Confusion matrix",
                info.text = "Neural net clustering",
                caption = ""
              ),
              plot_deepnet_lossplot_ui(
                ns("deepnet_lossplot"),
                title = "Loss history",
                info.text = "Plot of loss history",
                caption = ""
              )
            )
          )
        )
      ),  ## end of tabPanel

      ##----------------------------------------------------------------
      shiny::tabPanel(
        "Gradient vs. foldchange",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>Gradient vs foldchange. </b>. This board compares input gradients of the network (i.e. the change of prediction with respect to inputs) with the log-foldchange. Good biomarkers should have both high foldchange and input gradient.")),
          bslib::layout_columns(
            col_widths = c(7,5),
            height = "calc(100vh - 180px)",            
            bslib::layout_columns(
              col_widths = c(12,12),
              plot_deepnet_gradients_ui(
                ns("deepnet_gradients"),
                title = "Network gradients",
                info.text = "SNF affinity matrices",
                caption = ""
              ),
              plot_deepnet_gradients_ui(
                ns("deepnet_fcvsgrad"),
                title = "Gradient vs. foldchange",
                info.text = "Foldchange vs. gradient",
                caption = "",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              )
            ),
            table_deepnet_gradients_ui(
              ns("deepnet_table"),
              title = "Network gradients",
              info.text = "SNF affinity matrices",
              caption = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      )  ## end of tabPanel
      
    )
  )
}
