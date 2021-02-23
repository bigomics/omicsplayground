TcgaInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
    )
}

TcgaUI <- function(id) {
    ns <- NS(id)  ## namespace
    fillCol(
        height = 750,
        tabsetPanel(
            id = ns("tabs1"),
            tabPanel("TCGA survival", uiOutput(ns("TCGAanalysis_UI")))
        )
    )
}

TcgaBoard <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE
    fullH = 750       # row height of panel 
    tabH = '70vh'
    
    ## reactive functions from shared environment
    inputData <- env[["load"]][["inputData"]]
    selected_gxmethods <- env[["expr"]][["selected_gxmethods"]]
    selected_gsetmethods <- env[["enrich"]][["selected_gsetmethods"]]
    
    description = "<b>TCGA analysis (beta)</b>. Correlate your signature with the survival in cancer patients from the TCGA database. Warning: EXPERIMENTAL."
    output$description <- renderUI(HTML(description))

    tcga_infotext =
        "This <strong>TCGA analysis module</strong> computes the survival probability in (more than 10000) cancer patients of 32 TCGA cancer types, for your selected contrast. Each cohort is dichotomized into positively and negatively correlated with your signature. The survival probabilities are computed and tested using the Kaplan-Meier method.

"

    
    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================
    

    output$inputsUI <- renderUI({
        ui <- tagList(
            tipify( actionLink(ns("tcga_info"), "Info", icon = icon("info-circle")),
                   "Show more information about this module"),
            hr(), br(),             
            tipify( selectInput(ns('tcga_contrast'),'Contrast:', choices=NULL, multiple=FALSE),
                   "Select the contrast that you want to correlate with survival.",
                   placement="right", options = list(container = "body")
                   ),

            tipify( actionLink(ns("tcga_options"), "Options",
                               icon=icon("cog", lib = "glyphicon")),
                   "Toggle advanced options.",
                   placement="top", options = list(container = "body"))
            ## br(),br(),            
            ## conditionalPanel(
            ##     "input.tcga_options % 2 == 1", ns=ns,
            ##     tipify( selectInput(ns('tcga_profiledb'),"Profile DB:",
            ##                         choices=NULL, multiple=TRUE),
            ##            "Select external database for reference profiles.",
            ##            placement="top", options = list(container = "body"))
            ## )
        )
        ui
    })
    outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!
        
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    observeEvent( input$tcga_info, {
        showModal(modalDialog(
            title = HTML("<strong>TCGA Analysis Board</strong>"),
            HTML(tcga_infotext),
            easyClose = TRUE, size="l" ))
    })

    ## update choices upon change of data set 
    observe({
        ngs <- inputData()
        ##req(ngs)
        if(is.null(ngs)) return(NULL)
        comparisons <- colnames(ngs$model.parameters$contr.matrix)
        comparisons <- sort(comparisons)
        updateSelectInput(session, "tcga_contrast", choices=comparisons,
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
    
    tcga_tcgasurv.RENDER %<a-% reactive({

        ngs <- inputData()
        req(ngs)
        
        contrast = 1
        contrast <- input$tcga_contrast
        res = pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        names(res)
        sig <- res$fc[,contrast]
        ##matrix_file = search.path(c(FILES,ARCHS4.DIR),"tcga_matrix.h5")
        ##FILESX <- sub("lib$","libx",FILES)
        matrix_file = search.path(c(FILES,FILESX),"tcga_matrix.h5")
        if(!file.exists(matrix_file)) {
            showNotification("FATAL ERROR: could not find tcga_matrix.h5")
            return(NULL)
        }
        dbg("[tcga_tcgasurv.RENDER] reading TCGA data from",matrix_file)
        showNotification("computing survival probabilities...")
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
    
    tcga_tcgasurv.opts = tagList(
        tipify( checkboxInput(ns("tcga_surv_deceasedonly"), "deceased only", FALSE),
               "Only include deceased cases in survival analysis, i.e. exclude censored cases (patients still alive at evaluation time). This compares strictly the deceased cases, early vs late.",
               placement="left", options = list(container = "body")),
        ## radioButtons(ns("tcga_tcgasurv_sortby"), "sort by:", c("name","p-value"), inline=TRUE),
        tipify(radioButtons(ns("tcga_tcgasurv_ntop"), "N cor genes:", c(25,100,250,1000),
                            selected=100, inline=TRUE),
               "Number of top genes for calculating the correlation.", placement="left",
               options = list(container = "body"))        
    )
    tcga_tcgasurv_info = "<strong>TCGA survival analysis.</strong> Survival probability of cancer patients in 32 TCGA cancer types. Each cohort is dichotomized into positively and negatively correlated with your signature. The survival probabilities are computed and tested using the Kaplan-Meier method."
    tcga_tcgasurv_caption = tcga_tcgasurv_info
    
    callModule(
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
        height = c(fullH-80, 700), width = c("auto",1350),
        res=c(72,85)
    )

    
    ##================================================================================
    ##========================= OUTPUT UI ============================================
    ##================================================================================

    output$TCGAanalysis_UI <- renderUI({
        fillCol(
            ## id = ns("expr_topgenes"),
            height = fullH,
            flex = c(NA,0.05,1), ##height = 370,
            div(HTML(tcga_tcgasurv_caption), class="caption"),
            br(),
            plotWidget(ns("tcga_tcgasurv"))
        )
    })
    outputOptions(output, "TCGAanalysis_UI", suspendWhenHidden=FALSE) 
    
} ## end-of-Board 
