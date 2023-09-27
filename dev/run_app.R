setwd('components/app/R')
source('global.R')
source('ui.R')
source('server.R')
shinyApp(
    ui = app_ui,
    server = app_server,
    uiPattern = '.*',
    options = list(
        launch.browser = TRUE,
        port = 3838
    )
)
