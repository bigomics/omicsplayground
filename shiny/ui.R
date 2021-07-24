help.tabs <- navbarMenu(
	"Help",
	tabPanel(title=HTML("<a href='https://omicsplayground.readthedocs.io' target='_blank'>Documentation")),
	tabPanel(title=HTML("<a href='https://www.youtube.com/watch?v=_Q2LJmb2ihU&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-' target='_blank'>Video tutorials</a>")),
	tabPanel(title=HTML("<a href='https://github.com/bigomics/omicsplayground' target='_blank'>GitHub")),
	tabPanel(title=HTML("<a href='https://hub.docker.com/r/bigomics/omicsplayground' target='_blank'>Docker")),
	tabPanel(title=HTML("<a href='https://groups.google.com/d/forum/omicsplayground' target='_blank'>Google groups"))
)

TABVIEWS <- list(
	"load"   = tabView("Home",LoadingInputs("load"),LoadingUI("load")),
	"view"   = tabView("DataView",DataViewInputs("view"),DataViewUI("view")),
	"clust"  = tabView("Clustering",ClusteringInputs("clust"),ClusteringUI("clust")),
	"wgcna"  = tabView("WGCNA (beta)",WgcnaInputs("wgcna"),WgcnaUI("wgcna")),
	"expr"   = tabView("Differential expression",ExpressionInputs("expr"),ExpressionUI("expr")),
	"cor"    = tabView("Correlation analysis", CorrelationInputs("cor"), CorrelationUI("cor")),
	"enrich" = tabView("Geneset enrichment",EnrichmentInputs("enrich"), EnrichmentUI("enrich")),
	"func"   = tabView("Pathway analysis", FunctionalInputs("func"), FunctionalUI("func")),
	"word"   = tabView("Word cloud", WordCloudInputs("word"), WordCloudUI("word")),
	"drug"   = tabView("Drug connectivity", DrugConnectivityInputs("drug"), DrugConnectivityUI("drug")),
	"isect"  = tabView("Compare signatures", IntersectionInputs("isect"), IntersectionUI("isect")),
	"sig"    = tabView("Test signatures", SignatureInputs("sig"), SignatureUI("sig")),
	"bio"    = tabView("Find biomarkers", BiomarkerInputs("bio"), BiomarkerUI("bio")),
	"cmap"   = tabView("Similar experiments", ConnectivityInputs("cmap"), ConnectivityUI("cmap")),
	"scell"  = tabView("CellProfiling", SingleCellInputs("scell"), SingleCellUI("scell")),
	"tcga"   = tabView("TCGA survival (beta)", TcgaInputs("tcga"), TcgaUI("tcga")),
	"comp"   = tabView("Compare datasets (beta)", CompareInputs("comp"), CompareUI("comp"))
)

if(DEV) {
	TABVIEWS$corsa  = tabView("CORSA (dev)",CorsaInputs("corsa"),CorsaUI("corsa"))
	TABVIEWS$system = tabView("Systems analysis (dev)",SystemInputs("system"),SystemUI("system"))
	TABVIEWS$multi  = tabView("Multi-level (dev)", MultiLevelInputs("multi"), MultiLevelUI("multi"))
}

names(TABVIEWS)
TABVIEWS <- TABVIEWS[names(TABVIEWS) %in% names(which(ENABLED))]
names(TABVIEWS)

logout.tab <- tabPanel(title=HTML("<a id='logout' href='/logout'>Logout"))

createUI <- function(tabs)
{
	message("\n======================================================")
	message("======================= UI ===========================")
	message("======================================================\n")
	
	version <- scan("../VERSION", character())[1]
	TITLE = paste(opt$TITLE,version)
	LOGO = div(img(src="bigomics-logo-white-48px.png", height="48px"),
						 TITLE, id="navbar-logo", style="margin-top:-13px;")    
	title = tagList(LOGO)
	windowTitle = TITLE
	theme = shinythemes::shinytheme("cerulean")
	id = "maintabs"
	##selected = "Home"    
	header = tagList(
		tags$head(tags$link(rel = "stylesheet", href = "playground.css")),
		tags$head(tags$link(rel="shortcut icon", href="favicon.ico")),
		shinyjs::useShinyjs(),
		firebase::useFirebase(),
		TAGS.JSSCRIPT,
		tags$script(async=NA, src="https://platform.twitter.com/widgets.js"),
		div(textOutput("current_dataset"),class='current-data')
		##QuestionBoard_UI("qa")
	)
	names(header) <- NULL
	
	busy.img = sample(dir("www/busy",pattern=".gif$",full.name=TRUE))[1]
	busy.img
	busy.img = "www/busy.gif"
	footer.gif = tagList(
		busy_start_up(
			text = "\nPrepping your Omics Playground...", mode = "auto",
			background="#2780e3", color="#ffffff",
			loader = img(src=base64enc::dataURI(file=busy.img))
		)
	)
	## if(runif(1) < 0.1)
	footer = footer.gif ## every now and then show easter egg..
	
	##-------------------------------------
	## create TAB list
	##-------------------------------------
	createNavbarMenu <- function(title, tabs, icon=NULL) {
		tablist <- TABVIEWS[tabs]
		names(tablist) <- NULL
		do.call( navbarMenu, c(tablist, title=title, icon=icon) )
	}
	##tablist <- TABVIEWS[tabs]
	tablist <- list()
	i=1
	for(i in 1:length(tabs)) {
		itab <- tabs[[i]]
		itab <- itab[which(ENABLED[itab])] ## only enabled
		if(length(itab)>1) {
			message("[MAIN] creating menu items for: ",paste(itab,collapse=" "))
			m <- createNavbarMenu( names(tabs)[i], itab )
			tablist[[i]] <- m
		} else if(length(itab)==1) {
			message("[MAIN] creating menu item for: ",itab)
			tablist[[i]] <- TABVIEWS[[itab]] 
		} else {
			
		}
	}
	tablist <- tablist[!sapply(tablist,is.null)]
	
	## add help menu
	tablist[["helpmenu"]] <- help.tabs
	if(SHINYPROXY) {
		tablist[["logout"]] <- logout.tab
	}
	
	##-------------------------------------
	## create navbarPage
	##-------------------------------------
	selected = "Home"    
	names(tablist) <- NULL
	do.call( navbarPage, c(tablist,
												 title=title, id=id,
												 selected=selected,
												 windowTitle = windowTitle,
												 header = tagList(header),
												 footer = tagList(footer),
												 theme = theme))
}

tabs = list(
	"Home" = c("load"),
	"DataView" = "view",
	"Clustering" = c("clust","wgcna"),
	"Expression" = c("expr","cor"),
	"Enrichment" = c("enrich","func","word","drug"),
	"Signature" = c("isect","comp","sig","bio","cmap","tcga"),
	"CellProfiling" = "scell",
	"DEV" = c("corsa","system","multi")
)
createUI(tabs)