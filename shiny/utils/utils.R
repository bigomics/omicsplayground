envcat <- function(var) {
	message(var," = ",Sys.getenv(var))
}

# TODO: this is a version that respects the DEBUG flag, but it's not being used
dbg <- function(...) {
	if(DEBUG) {
		msg <- list(...)
		msg <- paste(sapply(msg, function(s) paste(s,collapse=" ")),collapse=" ")
		message(msg)
	}
}

dbg <- function(...) {
	##msg = paste0(ifelse(is.null(module),"",paste0("<",module,"> ")),msg)
	msg = sapply( list(...),paste,collapse=" ")
	message(paste0("DBG ",sub("\n$","",paste(msg,collapse=" "))))
}

tipify2 <- function(...) {
	shinyBS::tipify(..., placement="top", options = list(container = "body"))
}
tipifyL <- function(...) {
	shinyBS::tipify(..., placement="left", options = list(container = "body"))
}
tipifyR <- function(...) {
	shinyBS::tipify(..., placement="right", options = list(container = "body"))
}
tipifyT <- function(...) {
	shinyBS::tipify(..., placement="top", options = list(container = "body"))
}
tipifyB <- function(...) {
	shinyBS::tipify(..., placement="bottom", options = list(container = "body"))
}

# TODO: this function isn't being used
premium.feature <- function(...) {
	message("[premium.feature] USER_MODE = ",USER_MODE)
	message("[premium.feature] DEV = ",DEV)        
	el <- list(...)
	if(USER_MODE %in% c("pro","premium","dev")) return(el)
	tipify(disabled(...),
				 "This is a Premium feature. Upgrade to enable this feature."
	)    
	
}

# TODO: this function can be a variable
in.shinyproxy <- function() {
	## Determine if we are in ShinyProxy
	##
	vars <- c("SHINYPROXY_USERNAME","SHINYPROXY_USERGROUPS",
						"PLAYGROUND_USERID","PLAYGROUND_LEVEL")
	vars <- c("SHINYPROXY_USERNAME","SHINYPROXY_USERGROUPS")
	vals <- sapply(vars,Sys.getenv)
	all(vals!="")
}

tabRequire <- function(pgx, slot, tabname, subtab) {
	if(!slot %in% names(pgx)) {
		cat(paste("[MAIN] object has no ",slot," results. hiding tab.\n"))
		hideTab(tabname, subtab)
	} else {
		showTab(tabname, subtab)
	}
}

fileRequire <- function(file, tabname, subtab) {
	file1 <- search_path(c(FILES,FILESX),file)
	has.file <- !is.null(file1) && file.exists(file1)
	if(!has.file) {
		message(paste("[MAIN] file ",file," not found. Hiding",subtab,"\n"))
		hideTab(tabname, subtab)
	} else {
		message(paste("[MAIN] file ",file," available. Showing",subtab,"\n"))        
		showTab(tabname, subtab)
	}
}

## From https://github.com/plotly/plotly.js/blob/master/src/components/modebar/buttons.js
all.plotly.buttons = c(
	"toImage",
	"senDataToCloud","editInChartStudio","zoom2d","pan2d","select2d",
	"lasso2d","drawclosedpath","drawopenpath","drawline","drawrect",
	"drawcircle","eraseshape","zoomIn2d","zoomOut2d",
	"autoScale2d","resetScale2d","zoom3d","pan3d",
	"orbitRotation","tableRotation","resetCameraDefault3d",
	"resetCameraLastSave3d","hoverClosest3d","zoomInGeo",
	"zoomOutGeo","resetGeo","hoverClosestGeo","hoverClosestGl2d",
	"hoverClosestPie","resetViewSankey","toggleHover",
	"hoverClosestCartesian","hoverCompareCartesian",
	"resetViews","toggleSpikelines",
	"resetViewMapbox","zoomInMapbox","zoomOutMapbox"
)