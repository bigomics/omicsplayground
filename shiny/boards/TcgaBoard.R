TcgaInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

TcgaUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("TCGA survival", shiny::uiOutput(ns("TCGAanalysis_UI")))
        )
    )
}

TcgaBoard <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE
    fullH = 800       # row height of panel 
    tabH = '70vh'
    
    ## reactive functions from shared environment
    inputData <- env[["load"]][["inputData"]]
    selected_gxmethods <- env[["expr"]][["selected_gxmethods"]]
    selected_gsetmethods <- env[["enrich"]][["selected_gsetmethods"]]
    
    description = "<b>TCGA analysis (beta)</b>. Correlate your signature with the survival in cancer patients from the TCGA database. Warning: EXPERIMENTAL."
    output$description <- shiny::renderUI(shiny::HTML(description))

    tcga_infotext =
        "This <strong>TCGA analysis module</strong> computes the survival probability in (more than 10000) cancer patients of 32 TCGA cancer types, for your selected contrast. Each cohort is dichotomized into positively and negatively correlated with your signature. The survival probabilities are computed and tested using the Kaplan-Meier method.

"
    
    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================    

    output$inputsUI <- shiny::renderUI({
        ui <- shiny::tagList(
            tags$head(tags$style("#tcga-genelist.form-control {font-size:11px !important;padding:3px;height:200px;}")),
            shinyBS::tipify( shiny::actionLink(ns("tcga_info"), "Info", icon = shiny::icon("info-circle")),
                   "Show more information about this module"),
            shiny::hr(), shiny::br(),
            shinyBS::tipify( shiny::radioButtons( ns('sigtype'),"Signature type:",
                                 choices=c("contrast","genelist"),
                                 selected="contrast", inline=TRUE),
                   "number of top genes to show",
                    placement="right", options = list(container = "body")),            
            shiny::conditionalPanel(
                "input.sigtype == 'contrast'", ns=ns,
                shinyBS::tipify( shiny::selectInput(ns('contrast'),NULL, choices=NULL, multiple=FALSE),
                   "Select the contrast that you want to correlate with survival.",
                   placement="right", options = list(container = "body")
                   ),
            ),
            shiny::conditionalPanel(
                "input.sigtype == 'genelist'", ns=ns,
                shinyBS::tipify( shiny::textAreaInput(ns("genelist"), NULL, value = NULL,
                                      height = "100px", width = "100%", 
                                      rows=4, placeholder="Paste your custom gene list"),
                       "Paste a custom list of genes to be used as features.",
                       placement="bottom")
            ),
            shiny::br(),
            shinyBS::tipify( shiny::actionLink(ns("tcga_options"), "Options",
                               icon=icon("cog", lib = "glyphicon")),
                   "Toggle advanced options.",
                   placement="top", options = list(container = "body"))
            ## shiny::br(),br(),            
            ## shiny::conditionalPanel(
            ##     "input.tcga_options % 2 == 1", ns=ns,
            ##     shinyBS::tipify( shiny::selectInput(ns('tcga_profiledb'),"Profile DB:",
            ##                         choices=NULL, multiple=TRUE),
            ##            "Select external database for reference profiles.",
            ##            placement="top", options = list(container = "body"))
            ## )
        )
        ui
    })
    # shiny::outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!
        
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    shiny::observeEvent( input$tcga_info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>TCGA Analysis Board</strong>"),
            shiny::HTML(tcga_infotext),
            easyClose = TRUE, size="l" ))
    })

    ## update choices upon change of data set 
    shiny::observe({
        ngs <- inputData()
        ##req(ngs)
        if(is.null(ngs)) return(NULL)
        comparisons <- colnames(ngs$model.parameters$contr.matrix)
        comparisons <- sort(comparisons)
        shiny::updateSelectInput(session, "contrast", choices=comparisons,
                          selected = head(comparisons,1))
        
    })
    
    ##================================================================================
    ## TCGA survival
    ##================================================================================

    search.path <- function(paths, file) {
        dir <- paths[which(file.exists(file.path(paths,file)))]
        if(length(dir)==0) return(NULL)
        file.path(dir[1],file)
    }
    
     tcga_tcgasurv.RENDER %<a-% shiny::reactive({

        ngs <- inputData()
        shiny::req(ngs)

        if(input$sigtype == 'contrast') {
            contrast = 1
            contrast <- input$contrast
            shiny::req(contrast)
            res = pgx.getMetaFoldChangeMatrix(ngs, what="meta")
            names(res)
            sig <- res$fc[,contrast]
        } else if(input$sigtype == 'genelist') {
            shiny::req(input$genelist)
            genes <- as.character(input$genelist)
            genes <- strsplit(genes, split='[\t, \n]')[[1]]
            genes <- gsub("[ ]","",genes)
            sig <- rownames(ngs$X) %in% genes 
            names(sig) <- rownames( ngs$X)
        } else {
            stop("[tcga_tcgasurv.RENDER] invalid sigtype")
        }

        ##matrix_file = search.path(c(FILES,ARCHS4.DIR),"tcga_matrix.h5")
        ##FILESX <- sub("lib$","libx",FILES)
        matrix_file = search.path(c(FILES,FILESX),"tcga_matrix.h5")
        if(!file.exists(matrix_file)) {
            shiny::showNotification("FATAL ERROR: could not find tcga_matrix.h5")
            return(NULL)
        }
        dbg("[tcga_tcgasurv.RENDER] reading TCGA data from",matrix_file)
        shiny::showNotification("computing survival probabilities...")
        sortby="name";ntop=100
        sortby <- input$tcga_tcgasurv_sortby
        ntop <- as.integer(input$tcga_tcgasurv_ntop)
        deceased.only <- input$tcga_surv_deceasedonly
        
        pgx.testTCGAsurvival(
            sig,
            matrix_file,
            lib.dir = FILES,
            ntop = ntop,
            sortby.p = FALSE,
            deceased.only = deceased.only,
            min.cases = 10            
        )

    })
    
    tcga_tcgasurv.opts = shiny::tagList(
        shinyBS::tipify( shiny::checkboxInput(ns("tcga_surv_deceasedonly"), "deceased only", FALSE),
               "Only include deceased cases in survival analysis, i.e. exclude censored cases (patients still alive at evaluation time). This compares strictly the deceased cases, early vs late.",
               placement="left", options = list(container = "body")),
        ## shiny::radioButtons(ns("tcga_tcgasurv_sortby"), "sort by:", c("name","p-value"), inline=TRUE),
        shinyBS::tipify(shiny::radioButtons(ns("tcga_tcgasurv_ntop"), "N cor genes:", c(25,100,250,1000),
                            selected=100, inline=TRUE),
               "Number of top genes for calculating the correlation.", placement="left",
               options = list(container = "body"))        
    )
    tcga_tcgasurv_info = "<strong>TCGA survival analysis.</strong> Survival probability of cancer patients in 32 TCGA cancer types. Each cohort is dichotomized into positively and negatively correlated with your signature. The survival probabilities are computed and tested using the Kaplan-Meier method."
    tcga_tcgasurv_caption = tcga_tcgasurv_info
    
    shiny::callModule(
        plotModule,
        id = "tcga_tcgasurv",
        title = "TCGA survival analysis", label="a",
        func = tcga_tcgasurv.RENDER,
        func2 = tcga_tcgasurv.RENDER,
        info.text = tcga_tcgasurv_info,
        ## caption = tcga_tcgasurv_info,
        options = tcga_tcgasurv.opts,
        download.fmt = c("pdf","png"),
        pdf.width = 15, pdf.height = 10,
        height = c(fullH, 750), width = c("auto",1400),
        res=c(80,85),
        add.watermark = WATERMARK
    )

    
    ##================================================================================
    ##========================= OUTPUT UI ============================================
    ##================================================================================

    output$TCGAanalysis_UI <- shiny::renderUI({
        shiny::fillCol(
            ## id = ns("expr_topgenes"),
            height = fullH,
            flex = c(NA,0.02,1), ##height = 370,
            shiny::div(shiny::HTML(tcga_tcgasurv_caption), class="caption"),
            shiny::br(),
            plotWidget(ns("tcga_tcgasurv"))
        )
    })
    # shiny::outputOptions(output, "TCGAanalysis_UI", suspendWhenHidden=FALSE) 
    
} ## end-of-Board 
