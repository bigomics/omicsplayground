setwd('components/app/R')
source('global.R')
source('ui.R')
source('server.R')
shinyApp(
    ui = app_ui,
    server = app_server,
    uiPattern = '.*',
    options = list(
        launch.browser = FALSE,
        host = "0.0.0.0",
        port = 3838
    )
)
