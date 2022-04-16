## DEAN ATTALI code recommendation:
## example of how a board module should be written (TcgaBoard) 
##
## https://github.com/bigomics/omicsplayground/pull/20/commits/bd943d84d316d76dca9140f2fd3610b3d1dfc950



TcgaInputs <- function(id) {
	ns <- NS(id)
	
	bigdash::tabSettings(
		withTooltip(
			actionLink(ns("tcga_info"), "Info", icon = icon("info-circle")),
			"Show more information about this module"
		),
		
		hr(),
		br(),
		
		withTooltip(
			radioButtons(
				ns("sigtype"),
				"Signature type:",
				choices = c("contrast","genelist"),
				selected = "contrast",
				inline = TRUE
			),
			"number of top genes to show",
			placement = "right",
			options = list(container = "body")
		),
		
		conditionalPanel(
			"input.sigtype == 'contrast'",
			ns = ns,
			withTooltip(
				selectInput(ns("contrast"), NULL, choices = NULL, multiple = FALSE),
				"Select the contrast that you want to correlate with survival.",
				placement = "right",
				options = list(container = "body")
			),
		),
		
		conditionalPanel(
			"input.sigtype == 'genelist'",
			ns = ns,
			withTooltip(
				textAreaInput(
					ns("genelist"),
					NULL,
					value = NULL,
					height = "100px",
					width = "100%",
					rows = 4,
					placeholder = "Paste your custom gene list"
				),
				"Paste a custom list of genes to be used as features.",
				placement = "bottom"
			)
		),
		
		br(),
		
		withTooltip(
			actionLink(ns("tcga_options"), "Options", icon = icon("cog", lib = "glyphicon")),
			"Toggle advanced options.",
			placement = "top",
			options = list(container = "body")
		)
	)
}

TcgaUI <- function(id) {
	ns <- NS(id)
	
	fillCol(
		height = 750,
		tabsetPanel(
			id = ns("tabs1"),
			tabPanel(
				"TCGA survival",
				fillCol(
					height = 800,
					flex = c(NA, 0.02, 1),
					div(
						class = "caption",
						tags$strong("TCGA survival analysis."),
						"Survival probability of cancer patients in 32 TCGA cancer types.",
						"Each cohort is dichotomized into positively and negatively correlated with your signature.",
						"The survival probabilities are computed and tested using the Kaplan-Meier method."
					),
					br(),
					plotWidget(ns("tcga_tcgasurv"))
				)
			)
		)
	)
}
