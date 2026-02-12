##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


if(0) {

  install.packages("xkcd")
  library(extrafont)
  tryCatch({
    temp_font <- file.path(tempdir(), "xkcd.ttf")
    download.file("https://toledoem.github.io/img/xkcd.ttf",
      destfile = temp_font, mode = "wb", timeout = 60)
  }, error = function(e) {
    warning("Failed to download xkcd font. See https://github.com/ipython/xkcd-font")
  })
  
  # If downloaded, copy into the user's fonts directory for registration (not packaged)
  fonts_dir <- path.expand("~/.fonts")
  if (!dir.exists(fonts_dir)) dir.create(fonts_dir, recursive = TRUE)
  if (exists("temp_font") && file.exists(temp_font)) {
    file.copy(temp_font, file.path(fonts_dir, "xkcd.ttf"), overwrite = TRUE)
  }
  
  # Register fonts (import only when needed)
  font_import(pattern = "[X/x]kcd", prompt = FALSE)
  fonts()
  fonttable()
  if (.Platform$OS.type != "unix") {
    loadfonts(device = "win")
  } else {
    loadfonts()
  }

}



#' The application server-side logic
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @export
prism_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    get_dataframe <- function() {
      within(mtcars, {
        vs <- factor(vs)
        am <- factor(am)
        cyl  <- factor(cyl)
        gear <- factor(gear)
      })      
    }
        
    last_plotcode <- ""

    ##get_plotcode <- eventReactive( input$chartbot_user_input, {
    get_plotcode <- eventReactive( input$chartbot_send, {        

      dbg("chartbot_user_input reacted!") 
      
      msg <- "Plot weight versus mpg, color by cylinder. Apply ggprism theme. Use large font"
      msg <- "Make the dots extra large"
      msg <- "Replace dots with large triangles"
      msg <- "Label all japanese cars"      

      msg <- input$chartbot_user_input
      if(msg == "") return(NULL)
      if(msg == "reset") {
        last_plotcode <<- ""
        shiny::updateTextInput(session, "chartbot_user_input", value = "",
          placeholder = "What do you want to plot?")
        empty <- "ggplot() + theme_void()"
        return(empty)
      }
      
      data <- get_dataframe()
      vars <- paste(colnames(data),collapse=", ")
      rows <- paste(rownames(data),collapse=", ")      

      pointsize=4;fontsize=12;theme="classic"
      pointsize <- input$pointsize
      fontsize <- input$fontsize
      theme <- input$theme
      
      msg2 <- paste(
        "You are asked the following request about modifying a ggplot figure.", 
        "Just give the raw plotting code. Do not give explanations. Prepend package namespace if the function is not from ggplot2.",
        "\nThis is the request of the user: ", msg,        
        "\nThese are the variables in the dataframe called 'data': ", vars,
        "\nThese are the rownames of 'data': ", rows,        
        "\nThis is the last plotting code of the graph: ", last_plotcode,
        "\nDefault point size unless asked by user: ", pointsize,
        "\nDefault font size unless asked by user: ", fontsize,
        "\nDefault theme unless asked by user: ", theme                
      )

      plotcode <- playbase::ai.ask(msg2, "groq:openai/gpt-oss-20b")
      plotcode <- paste0(plotcode,"\n")
      plotcode <- gsub("```[rR]|```","",plotcode)
      ##plotcode <- sub(".*ggplot\\(","ggplot(",plotcode)      
            
      last_plotcode <<- plotcode      
      shiny::updateTextInput(session, "chartbot_user_input", value = "",
        placeholder = "Any edits?")
      ##shinychat::chat_clear("chartbot")

      return(plotcode)      
    })

    
    output$plot1 <- renderPlot({
      plotcode <- get_plotcode()
      shiny::req(plotcode)
      data <- get_dataframe()
      dbg("WARNING: performing eval(): plotcode = \n", plotcode)
      require(ggplot2)
      eval(parse(text=plotcode))  
    })

    output$plotcode <- renderUI({
      plotcode <- get_plotcode()
      shiny::req(plotcode)
      plotcode <- sub("ggplot\\(","ggplot(<br>&nbsp;&nbsp;",plotcode)
      plotcode <- gsub("[+]","+<br>&nbsp;",plotcode)      
      HTML(plotcode)      
    })
    
    
  })
}
