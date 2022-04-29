##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

UploadBoard <- function(id,
                        pgx_dir,
                        pgx,
                        auth,
                        limits = c("samples"=1000,"comparisons"=20,
                                   "genes"=20000, "genesets"=10000,
                                   "datasets"=10),
                        enable_upload = TRUE,
                        enable_save = TRUE,
                        enable_userdir = TRUE
                        )
{
  moduleServer(id, function(input, output, session) 
  {
    ns <- session$ns ## NAMESPACE
    dbg("[UploadBoard] >>> initializing UploadBoard...")

    loadedDataset <- shiny::reactiveVal(0)  ## counts/trigger dataset upload
    
    message("[UploadBoard] in.shinyproxy = ",in.shinyproxy())    
    message("[UploadBoard] SHINYPROXY_USERNAME = ",Sys.getenv("SHINYPROXY_USERNAME"))
    message("[UploadBoard] SHINYPROXY_USERGROUPS = ",Sys.getenv("SHINYPROXY_USERGROUPS"))
    message("[UploadBoard] pgx_dir = ",pgx_dir)
    
    dbg("[UploadBoard] getwd = ",getwd())      

    
    output$navheader <- shiny::renderUI({
      fillRow(
        flex=c(NA,1,NA),
        ##h2(input$nav),
        shiny::div(
          id="navheader-current-section",
          HTML("Upload data &nbsp;"), 
          shiny::actionLink(
            ns("module_info"), "",
            icon = shiny::icon("info-circle"),
            style = "color: #ccc;"
            )
        ),        
        shiny::br(),
        shiny::div(pgx$name, id="navheader-current-dataset")
      )
    })

    shiny::observeEvent( input$module_info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Upload data</strong>"),
            shiny::HTML(module_infotext),
            easyClose = TRUE, size="l" ))
    })

    module_infotext =paste0(
        'Under the <b>Upload data</b> panel users can upload their transcriptomics and proteomics data to the platform. The platform requires 3 data files as listed below: a data file containing counts/expression (counts.csv), a sample information file (samples.csv) and a file specifying the statistical comparisons as contrasts (contrasts.csv). It is important to name the files exactly as shown. The file format must be comma-separated-values (CSV) text. Be sure the dimensions, row names and column names match for all files. On the left side of the panel, users need to provide a unique name and brief description for the dataset while uploading. N.B. Users can now create contrasts from the platform itself, so the contrasts.csv file is optional.

<br><br>
<ol>
<li>counts.csv: Count/expression file with gene on rows, samples as columns.
<li>samples.csv: Samples file with samples on rows, phenotypes as columns.
<li>contrasts.csv: Contrast file with conditions on rows, contrasts as columns.
</ol>

<br><br><br>
<center><iframe width="560" height="315" src="https://www.youtube.com/embed/elwT6ztt3Fo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe><center>

')
    
    ##================================================================================
    ##====================== NEW DATA UPLOAD =========================================
    ##================================================================================
    ##reload_pgxdir()
    
    getPGXDIR <- shiny::reactive({
        ##reload_pgxdir()  ## force reload

        email="../me@company.com"
        email <- auth$email()
        email <- gsub(".*\\/","",email)
        pdir  <- pgx_dir  ## from module input

        ##USERDIR=FALSE
        if(enable_userdir) {
            pdir <- paste0(pdir,"/",email)
            if(!is.null(email) && !is.na(email) && email!="") pdir <- paste0(pdir,'/')
            if(!dir.exists(pdir)) {
                dbg("[LoadingBoard:getPGXDIR] userdir does not exists. creating pdir = ",pdir)
                dir.create(pdir)
                dbg("[LoadingBoard:getPGXDIR] copy example pgx")                
                file.copy(file.path(pgx_dir,"example-data.pgx"),pdir)
            }
        }
        pdir
    })
    
    if(enable_upload) {

        uploaded_pgx <- UploadModuleServer(
            id = "upload_panel",
            FILES = FILES,
            pgx.dirRT = shiny::reactive(getPGXDIR()),
            height = 720,
            ## limits = c(samples=20, comparisons=20, genes=8000),
            limits = limits
        )
        
        shiny::observeEvent( uploaded_pgx(), {
            
            dbg("[observe::uploaded_pgx] uploaded PGX detected!")
            new_pgx <- uploaded_pgx()
            
            dbg("[observe::uploaded_pgx] initializing PGX object")
            new_pgx <- pgx.initialize(new_pgx)
            
            ## update Session PGX
            dbg("[UploadBoard@load_react] **** copying current pgx to session.pgx  ****")        
            for(i in 1:length(new_pgx)) {
                pgx[[names(new_pgx)[i]]] <- new_pgx[[i]]
            }

            DT::selectRows(proxy = DT::dataTableProxy(ns("pgxtable")), selected=NULL)            

            savedata_button <- NULL
            if(enable_save) {                

                dbg("[UploadBoard] observeEvent:savedata reacted")        
                ## -------------- save PGX file/object ---------------
                pgxname <- sub("[.]pgx$","",new_pgx$name)
                pgxname <- gsub("^[./-]*","",pgxname)  ## prevent going to parent folder
                pgxname <- paste0(gsub("[ \\/]","_",pgxname),".pgx")
                pgxname
                
                pgxdir  <- getPGXDIR()
                fn <- file.path(pgxdir,pgxname)
                fn <- iconv(fn, from = '', to = 'ASCII//TRANSLIT')
                
                ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ## Note: Currently we use 'ngs' as object name but want to go
                ## towards 'pgx' as standard name. Actually saving as RDS
                ## should be better.
                ngs=new_pgx
                save(ngs, file=fn)

                remove(ngs)
                remove(new_pgx)
                

                message("[UploadBoard::@savedata] updating PGXINFO")                    
                pgx.initDatasetFolder(pgxdir, force=FALSE, verbose=TRUE)
                ## reload_pgxdir(reload_pgxdir()+1)
            }
            
            ## shiny::removeModal()
            msg1 <- "<b>Ready!</b>"
            ##beepr::beep(sample(c(3,4,5,6,8),1))  ## music!!
            beepr::beep(10)  ## short beep
            
            if(enable_save) {
                msg1 <- "<b>Ready!</b><br>Your data is ready and has been saved in your library. You can now start exploring your data."
            } else {
                msg1 <- "<b>Ready!</b><br>Your data is ready. You can now start exploring your data."
            }
            loadedDataset(loadedDataset()+1)  ## notify new data uploaded

            showModal(
                modalDialog(
                    HTML(msg1),
                    title = NULL,
                    size = "s",
                    footer = tagList(
                        ## savedata_button,
                        ## shiny::actionButton(ns("sharedata"), "Share with others", icon=icon("share-alt")),
                        modalButton("Start!")
                    )
                ))
            ## updateTabsetPanel(session, "tabs",  selected = "Datasets")
        })

    }
    
    ##------------------------------------------------
    ## Board return object
    ##------------------------------------------------
    res <- list(
        loaded = loadedDataset
    )
    return(res)
  })
}
