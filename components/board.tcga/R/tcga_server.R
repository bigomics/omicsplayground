## DEAN ATTALI code recommendation:
## example of how a board module should be written (TcgaBoard) 
##
## https://github.com/bigomics/omicsplayground/pull/20/commits/bd943d84d316d76dca9140f2fd3610b3d1dfc950

TcgaBoard <- function(id, inputData)
{
    moduleServer(id, function(input, output, session)
    {
	ns <- session$ns
	fullH <- 800
	tabH <- "70vh"
	
	tcga_tcgasurv_info <- div(
		"This", tags$strong("TCGA analysis module"),
		"computes the survival probability in (more than 10000) cancer patients of 32 TCGA cancer types, for your selected contrast.",
		"Each cohort is dichotomized into positively and negatively correlated with your signature.",
		"The survival probabilities are computed and tested using the Kaplan-Meier method."
	)
	
	observeEvent(input$tcga_info, {
		showModal(
			modalDialog(
				title = tags$strong("TCGA Analysis Board"),
				tcga_tcgasurv_info,
				easyClose = TRUE,
				size = "l"
			)
		)
	})
	
	observe({
		ngs <- inputData()
		if (is.null(ngs)) return(NULL)
		comparisons <- colnames(ngs$model.parameters$contr.matrix)
		comparisons <- sort(comparisons)
		updateSelectInput(session, "contrast", choices = comparisons, selected = head(comparisons, 1))
	})
	
	tcga_tcgasurv_test <- reactive({
		matrix_file <- search_path(c(FILES, FILESX), "tcga_matrix.h5")
		if (is.null(matrix_file)) {
			showNotification("FATAL ERROR: could not find tcga_matrix.h5")
			return(NULL)
		}
		
		ngs <- inputData()
		req(ngs)
		
		if (input$sigtype == "contrast") {
			req(input$contrast)
			
			res <- pgx.getMetaFoldChangeMatrix(ngs, what = "meta")
			sig <- res$fc[, input$contrast]
		} else if (input$sigtype == "genelist") {
			req(input$genelist)
			
			genes <- as.character(input$genelist)
			genes <- strsplit(genes, split = "[\t, \n]")[[1]]
			genes <- gsub("[ ]", "", genes)
			sig <- rownames(ngs$X) %in% genes
			names(sig) <- rownames(ngs$X)
		}
		
		showNotification("computing survival probabilities...")
		pgx.testTCGAsurvival(
			sig,
			matrix_file,
			lib.dir = FILES,
			ntop = as.integer(input$tcga_tcgasurv_ntop),
			sortby.p = FALSE,
			deceased.only = input$tcga_surv_deceasedonly,
			min.cases = 10
		)
	})
	
	tcga_tcgasurv_opts <- tagList(
		withTooltip(
			checkboxInput(ns("tcga_surv_deceasedonly"), "deceased only", FALSE),
			paste(
				"Only include deceased cases in survival analysis,",
				"i.e. exclude censored cases (patients still alive at evaluation time).",
				"This compares strictly the deceased cases, early vs late."
			),
			placement = "left",
			options = list(container = "body")
		),
		
		withTooltip(
			radioButtons(ns("tcga_tcgasurv_ntop"), "N cor genes:", c(25, 100, 250, 1000), selected = 100, inline = TRUE),
			"Number of top genes for calculating the correlation.",
			placement = "left",
			options = list(container = "body")
		)
	)
	
	callModule(
		plotModule,
		id = "tcga_tcgasurv",
		title = "TCGA survival analysis",
		label = "a",
		func = tcga_tcgasurv_test,
		func2 = tcga_tcgasurv_test,
		info.text = as.character(tcga_tcgasurv_info),
		options = tcga_tcgasurv_opts,
		download.fmt = c("pdf", "png"),
		pdf.width = 15,
		pdf.height = 10,
		height = c(fullH, 750),
		width = c("auto", 1400),
		res = c(80, 85),
		add.watermark = WATERMARK
	)


    })
}
