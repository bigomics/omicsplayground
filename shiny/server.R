function(input, output, session) {
	
	message("\n========================================================")
	message("===================== SERVER ===========================")
	message("========================================================\n")
	
	message("[SERVER] calling boards...")
	message("[SERVER] USER_MODE = ", USER_MODE)
	server.start_time <- Sys.time()
	
	require(firebase)
	firebase=firebase2=NULL
	if(AUTHENTICATION=="firebase") {
		firebase  <- FirebaseEmailPassword$new()
		firebase2 <- FirebaseSocial$new()
	}
	
	## firebase <- NULL
	max.limits <- c("samples" = opt$MAX_SAMPLES,
									"comparisons" = opt$MAX_COMPARISONS,
									"genes" = opt$MAX_GENES,
									"genesets" = opt$MAX_GENESETS)
	env <- list()  ## communication "environment"
	env[["load"]]  <- callModule(
		LoadingBoard, "load", max.limits = max.limits,
		authentication = AUTHENTICATION, enable_delete = opt$ENABLE_DELETE,
		enable_save = opt$ENABLE_SAVE, firebase=firebase, firebase2=firebase2)
	
	##auth <- env[["load"]][['auth']]
	
	## load other modules if 
	if(ENABLED["view"])   env[["view"]]   <- callModule( DataViewBoard, "view", 
																											 inputData = env[["load"]][["inputData"]])
	if(ENABLED["clust"])  env[["clust"]]  <- callModule( ClusteringBoard, "clust", env)
	if(ENABLED["expr"])   env[["expr"]]   <- callModule( ExpressionBoard, "expr", env)
	if(ENABLED["enrich"]) env[["enrich"]] <- callModule( EnrichmentBoard, "enrich", env)
	if(ENABLED["func"])   env[["func"]]   <- callModule( FunctionalBoard, "func", env)
	if(ENABLED["word"])   env[["word"]]   <- callModule( WordCloudBoard, "word", env)
	if(ENABLED["drug"])   env[["drug"]]   <- callModule( DrugConnectivityBoard, "drug", env)
	if(ENABLED["isect"])  env[["isect"]]  <- callModule( IntersectionBoard, "isect", env)
	if(ENABLED["sig"])    env[["sig"]]    <- callModule( SignatureBoard, "sig", env)
	if(ENABLED["cor"])    env[["cor"]]    <- callModule( CorrelationBoard, "cor", env)
	if(ENABLED["bio"])    env[["bio"]]    <- callModule( BiomarkerBoard, "bio", env)
	if(ENABLED["cmap"])   env[["cmap"]]   <- callModule( ConnectivityBoard, "cmap", env)
	if(ENABLED["scell"])  env[["scell"]]  <- callModule( SingleCellBoard, "scell", env)
	if(ENABLED["tcga"])   env[["tcga"]]   <- callModule( TcgaBoard, "tcga", env)
	if(ENABLED["wgcna"])  env[["wgcna"]]  <- callModule( WgcnaBoard, "wgcna", env)
	if(ENABLED["comp"])   env[["comp"]]   <- callModule( CompareBoard, "comp", env)
	if(DEV) {            
		if(ENABLED["corsa"])  env[["corsa"]]  <- callModule( CorsaBoard, "corsa", env)
		if(ENABLED["system"]) env[["system"]] <- callModule( SystemBoard, "system", env)
		if(ENABLED["multi"])  env[["multi"]]  <- callModule( MultiLevelBoard, "multi", env)
		env[["qa"]] <- callModule( QuestionBoard, "qa", lapse = -1)
	}
	
	## message("[SERVER] all boards called:",paste(names(env),collapse=" "))
	message("[SERVER] boards enabled:",paste(names(which(ENABLED)),collapse=" "))
	## outputOptions(output, "clust", suspendWhenHidden=FALSE) ## important!!!
	
	output$current_dataset <- renderText({
		pgx <- env[["load"]][["inputData"]]()
		name <- gsub(".*\\/|[.]pgx$","",pgx$name)
		if(length(name)==0) name = "(no data)"
		name
	})
	
	## Dynamicall hide/show certain sections depending on USERMODE/object
	observe({
		pgx <- env[["load"]][["inputData"]]() ## trigger on change dataset
		
		## hide all main tabs until we have an object
		if(is.null(pgx)) {
			lapply(MAINTABS, function(m) hideTab("maintabs",m))
			if(!opt$ENABLE_UPLOAD)  hideTab("load-tabs","Upload data")
			if(is.null(ACCESS.LOG)) hideTab("load-tabs","Visitors map")            
			return(NULL)
		}
		
		message("[SERVER] dataset changed. reconfiguring menu...")
		## show all main tabs
		lapply(MAINTABS, function(m) showTab("maintabs",m))
		
		if(USER_MODE == "basic") {
			hideTab("maintabs","CellProfiling")
			hideTab("maintabs","DEV")
			hideTab("enrich-tabs1","GeneMap")
			hideTab("clust-tabs2","Feature ranking")
			hideTab("expr-tabs1","Volcano (methods)")
			hideTab("expr-tabs2","FDR table")
			hideTab("enrich-tabs1","Volcano (methods)")
			hideTab("enrich-tabs2","FDR table")
			hideTab("cor-tabs","Functional")    ## too slow
			hideTab("cor-tabs","Differential")  ## too complex
		}
		
		## hideTab("cor-tabs","Functional")       
		if(USER_MODE == "dev" || DEV) {
			showTab("maintabs","DEV")
			showTab("view-tabs","Resource info")
			showTab("enrich-tabs1","GeneMap")
			showTab("scell-tabs1","CNV")  ## DEV only
			showTab("scell-tabs1","Monocle") ## DEV only
			showTab("cor-tabs","Functional")
		} else {
			hideTab("maintabs","DEV")
			hideTab("view-tabs","Resource info")
			hideTab("enrich-tabs1","GeneMap")
			hideTab("scell-tabs1","CNV")  ## DEV only
			hideTab("scell-tabs1","Monocle") ## DEV only
			hideTab("cor-tabs","Functional")            
		}
		
		## Dynamically show upon availability in pgx object
		if(opt$ENABLE_UPLOAD) showTab("load-tabs","Upload data")            
		tabRequire(pgx, "connectivity", "maintabs", "Similar experiments")
		tabRequire(pgx, "drugs", "maintabs", "Drug connectivity")
		tabRequire(pgx, "wordcloud", "maintabs", "Word cloud")
		tabRequire(pgx, "deconv", "maintabs", "CellProfiling")
		fileRequire("tcga_matrix.h5", "maintabs", "TCGA survival (beta)")         
		if(!is.null(ACCESS.LOG)) showTab("load-tabs","Visitors map")                    
		
		message("[SERVER] reconfiguring menu done.")        
	})
	
	server.init_time <- round(Sys.time() - server.start_time, digits=4)    
	message("[SERVER] server.init_time = ",server.init_time," ",attr(server.init_time,"units"))
	
}