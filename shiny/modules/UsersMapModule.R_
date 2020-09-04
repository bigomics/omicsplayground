##
## This file is part of the Omics Playground project.
## Copyright (c) 2020 BigOmics Analytics Sagl. All rights reserved.
##

UsersMapInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList()
}

UsersMapUI <- function(id) {
    ns <- NS(id)  ## namespace
    fillCol(
        flex = c(1),
        height = 780,
        tabsetPanel(
            id = ns("tabs"),
            tabPanel("UsersMap",uiOutput(ns("usersmap_UI")))
        )
    )
}

UsersMapModule <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE

    fullH = 750
    rowH = 660  ## row height of panel
    tabH = 200  ## row height of panel
    tabH = '70vh'  ## row height of panel    
    
    description = "<b>UsersMap</b>. <br> Geo location of Omics Playground users."
    output$description <- renderUI(HTML(description))

    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================

    output$inputsUI <- renderUI({
        tagList()
    })
    ##outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!

    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    ##---------------------------------------------------------------
    ##------------- modules for UsersMap ---------------------------
    ##---------------------------------------------------------------
    geolocate.RENDER <- reactive({

        require(rgeolocate)

        accessfile = file.path(FILESX,"access-ncov2019.log")
        if(!file.exists(accessfile)) return(NULL)
                                            
        acc <- read.table(accessfile)
        ip <- as.character(acc[,1])
        ##loc <- ip_api(unique(ip))
        ip <- unique(ip)
        
        ##file <- system.file("extdata","GeoLite2-Country.mmdb", package = "rgeolocate")
        ##loc <- maxmind(ip, file, "country_code")
        file <- file.path(FILESX,"GeoLite2-City.mmdb")
        loc <- maxmind(ip, file, c("country_code", "country_name", "city_name"))
        country_code <- unique(loc$country_code)
        names(country_code) <- loc[match(country_code,loc$country_code),"country_name"]
        sort(table(results$country_code))
        sort(table(results$city_name))

        tt <- table(loc$country_name)
        df <- data.frame( country_name = names(tt),
                         country_code = country_code[names(tt)],
                         frequency = (as.integer(tt)))
        
        library(rworldmap)
        data(countryExData)
        head(countryExData)
        sPDF <- getMap()  
        mapCountryData(sPDF, nameColumnToPlot='continent')

        sPDF <- joinCountryData2Map(
            df,
            joinCode = "ISO2",
            nameJoinColumn = "country_code")
        
        par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
        mapCountryData(
            sPDF, nameColumnToPlot="frequency",
            mapTitle = "Number of unique IPs",
            numCats=10, catMethod="logFixedWidth")   
        
    })
    
    geolocate_info = "<strong>geolocate.</strong>."
    
    callModule(
        plotModule,
        id = "geolocate", label="a", 
        plotlib="baseplot", func=geolocate.RENDER, 
        info.text = geolocate_info,
        ##options = geolocate_options,
        pdf.width=12, pdf.height=7, pdf.pointsize=13,
        height = 0.5*rowH, res=72,
        ##datacsv = enrich_getWordFreq,
        title = "User's map"
    )


    ##---------------------------------------------------------------
    ##-------------- UI Layout for UsersMap ------------------------
    ##---------------------------------------------------------------

    geolocate_caption = "<b>(a)</b> <b>Geo locate.</b>"
    
    output$geolocate_UI <- renderUI({
        fillCol(
            height = fullH,
            flex = c(NA,1),
            div(HTML(geolocate_caption,"<br>"),class="caption"),
            plotWidget(ns("geolocate_map"))
        )
    })



}
