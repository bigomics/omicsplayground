##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##


QA <- list(
  "What are my top genes?" = "Your top differentially expressed genes across all comparisons are: {topgenes}.",
  "What is a 'meta-q' value?" = "The meta-q value is a aggregated q-value of multiple statistical methods",
  "How do I test my own gene set?" = "You can paste your list of genes/proteins in the 'Test genesets' panel. If you have many custom gene sets, you can upload custom GMT files at the final step of the upload wizard under the 'computation options'.",
  "Can I use LIMMA, EdgeR and DESeq2 for my proteomics data?" = "While these methods were originally conceived for the for differential expression analysis in RNA-seq data, there is increasingly more acceptance that these methods can also be used for proteomics data.",
  "How are duplicated gene/protein names handled in the counts file?" = "Duplicated row identifiers (genes/proteins with same name) are handled by averaging their intensities/counts. If you do not want this, be sure your feature names are unique.",
  "How are missing values handled?" = "Missing values can be imputed using SVD, or left as missing. If you prefer other imputation methods, please impute missing value before uploading to Omics Playground"
)

qa_items <- list()
for(i in 1:length(QA)) {
  q <- names(QA)[i]
  p <- bslib::accordion_panel(
    HTML(paste0("<b>",q,"</b>")),
    QA[[q]],
    value = paste0("qa_",i) )
  qa_items <- c(qa_items, list(p))
}

SearchWidgetUI <- function(id, label="FAQ", class="") {
  ns <- shiny::NS(id)
  tagList(
    shiny::actionButton(
      ns("show_faq"),
      label = label,
##      icon = icon("search"),
      class = paste("btn btn-outline-secondary search_btn",class),
      width = NULL
    ),
    uiOutput(ns("modalUI")),
    tags$style(".accordion-title {font-weight: 400;}")
  )
}

SearchWidgetServer <- function(
    id, pgx ) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    output$modalUI <- renderUI({
      bsutils::modal(
        id = ns("modal"),
        size = "lg",
        bsutils::modalHeader(
          shiny::div( style = "display: inline-block; font-size: 1.7em; padding: 0 15px 15px 0;",
            icon("search")),
          shiny::div( style = "display: inline-block; width: 100%;", 
            shiny::textInput(ns("question"), label=NULL,
              placeholder="Ask a question...",
              width = "100%"))
        ),
        bsutils::modalBody(
          # First shown by default
          bslib::accordion(id = ns("accordion"), !!!qa_items, multiple=FALSE)
        )
      )
    })
        
    observeEvent(input$show_faq, {
      bsutils::modal_show(ns("modal"))
    })
    
    score.order <- eventReactive( input$question, {
      q <- input$question
      q <- gsub("['?]","",q)
      q <- gsub("[-_;]"," ",q)
      keywords <- tolower(strsplit(q, split=' ')[[1]])
      score_question <- function(q) sum(unlist(sapply(keywords, agrep, q,
        ignore.case = TRUE, fixed=FALSE)))
      score <- sapply(tolower(names(QA)), score_question )
      order(score, decreasing=TRUE)
    })
    
    prepare_answer <- function(a) {
      markers <- playbase::pgx.getMarkerGenes(pgx)
      top.markers <- paste(head(names(sort(table(unlist(markers)),
        decreasing=TRUE))),collapse=" ")
      a <- sub("\\{topgenes\\}", top.markers, a)
      a
    }
    
    observeEvent( score.order(), {
      req(score.order())
      nn <- score.order()
      nn <- head(nn,8)
      for(i in 1:length(nn)) {
        j <- nn[i]
        answer <- prepare_answer(QA[[j]])      
        bslib::accordion_panel_update(
          id = "accordion", target=paste0("qa_",i),
          title = names(QA)[j], answer )
      }
    })
    
    # observe dataset and
  }) ## end of moduleServer
}
