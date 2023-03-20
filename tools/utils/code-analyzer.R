##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


setwd("~/Playground/omicsplayground")
rfiles <- dir("components", recursive=TRUE, patter=".*[.][rR]$",full.names=TRUE)

## scan all function declarations
func.defined <- c()
f <- rfiles[1]
for(f in rfiles) {
  src <- readLines(f)
  ## detect allw function declaration
  here.func <- grep("=[ ]*function|<-[ ]*function",src)
  these.func  <- trimws(sub("<-.*|=.*","",src[here.func]))
  these.func  <- gsub(".*[ ]|#","",these.func)
  if(length(these.func)) {
    ff <- cbind(these.func, f)
    func.defined <- rbind( func.defined, ff)
  }
}
colnames(func.defined) <- c("function.name","file")

## create regexpression for all functions
all.func <- func.defined[,"function"]
func.rexp <- paste0("[\\^ =-\\(]",all.func,"[\\(@,]")  ## NEED RETHINK
names(func.rexp) <- all.func

## detect which functions are used in each file
func.used <- c()
f <- rfiles[10]
for(f in rfiles) {
  ## detect functions used
  src <- readLines(f)
  using.func <- grep(paste(func.rexp,collapse="|"), src)
  using.func
  if(length(using.func)==0) next
  n=using.func[1]
  func.calls <- c()
  for(n in using.func){
    s <- src[n]
    f1 <- names(which(sapply(func.rexp, function(f) grepl(f,s))))
    func.calls <- c(func.calls, f1)
  }
  if(length(func.calls)>0) {
    tt <- table(func.calls)
    rr <- data.frame(f, "function.called"=names(tt), "nfreq"=as.integer(tt))
    func.used <- rbind(func.used, rr)
  }
}
colnames(func.used) <- c("file","function.name","nfreq")

head(func.defined,20)
head(func.used,20)

## detect multiple defined functions
ndefined <- table(func.defined[,"function.name"])
multiple.defined <- names(which(ndefined > 1))
ww <- tapply( func.defined[,"file"], func.defined[,"function.name"],
  function(w) paste(gsub(".*/","",sort(unique(w))),collapse=', '))
df1 <- data.frame( 'function.name'=names(ndefined), n.defined=as.integer(ndefined), where.defined=ww)
head(df1)
head(df1[which(df1$n.defined>1),],20)

## detect not used functions
head(func.defined,20)
head(func.used,20)

nused <- table(func.used[,"function.name"])
uu <- tapply( func.used[,"file"], func.used[,"function.name"],
  function(w) paste(gsub(".*/","",sort(unique(w))),collapse=', '))
df2 <- data.frame( 'function.name'=names(nused), n.used=as.integer(nused), where.used=uu)
head(df2)

df2 <- df2[match(df1$function.name, df2$function.name),]
df2$n.used[is.na(df2$n.used)] <- 0
df2$function.name <- NULL
df <- cbind(df1, df2)
rownames(df) <- NULL

df <- df[,c("function.name","n.defined","n.used","where.defined","where.used")]

## write.csv(df, file="code-analyzer-output.csv")


## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

library(shiny)
library(DT)

ui = fluidPage(
  h2("Code analytics: function declaration and calls"),
  div( class="row",
    actionButton("show_all","all"),
    actionButton("show_multi","multiple-defined"),
    actionButton("show_notused","not used")
  ),
  DT::DTOutput("dt")
)

server = function(input, output, session) {

  filtered_df <- reactiveVal(df)

  observeEvent( input$show_all, {
    filtered_df( df )
  })

  observeEvent( input$show_notused, {
    df1 <- df[ df$n.used==0,] 
    df1 <- df1[order(-df1$n.defined),]    
    filtered_df( df1)
  })
    
  observeEvent( input$show_multi, {
    df1 <- df[ df$n.defined > 1,] 
    df1 <- df1[order(-df1$n.defined),]    
    filtered_df(df1)
  })

  output$dt <- DT::renderDT(
    DT::datatable( filtered_df(),
      extensions = 'Scroller',
      options = list(
        deferRender = FALSE,
        dom = 't',
        ##columnDefs = list(list(className = 'dt-center', targets = 5)),
        scrollY = 800,
        scroller = TRUE,
        scrollX = TRUE,
        pageLength = 80)
    ) %>%
      formatStyle(
        'n.defined',
        backgroundColor = styleInterval(1, c('white', 'yellow'))
      ) %>%
      formatStyle(
        'n.used',
        backgroundColor = styleInterval(0, c('salmon', 'white'))
      ) 
  )
}

shiny::shinyApp(ui, server=server, options=list(launch.browser=TRUE))
