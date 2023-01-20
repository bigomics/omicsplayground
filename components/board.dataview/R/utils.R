plotly_default1 <- function(e) {
    e %>%
        plotly::layout(
            #xaxis = list(fixedrange=TRUE),
            #yaxis = list(fixedrange=TRUE),
            font = list(family = "Lato")
            ## title = pd$gene                    
        ) %>%
        ##plotly::config(displayModeBar = FALSE) %>%
        plotly::config(displaylogo = FALSE) %>%
        plotly::config(
            modeBarButtons = list(list("toImage","zoom2d","resetScale2d")),
            toImageButtonOptions = list(format='svg', height=500, width=900)
        )
}
