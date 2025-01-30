##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

UserProfileBoard <- function(id, auth, nav_count) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    dbg("[UserProfileBoard] >>> initializing UserBoard...")

    shiny::observeEvent(input$board_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>User Profile</strong>"),
        shiny::HTML(
          "The User Profile page allows you to view and change your
                subscription plan and to view the latest news about application
                development."
        ),
        easyClose = TRUE, size = "l"
      ))
    })

    output$plan <- renderUI({
      p(
        span(paste0("Subscription level: ", tools::toTitleCase(auth$level)))
      )
    })

    output$userdata <- renderTable(
      {
        dbg("[UserBoard::userdata]  renderDataTable")
        cl <- "badge badge-info"
        values <- c(
          Name = auth$username,
          Email = auth$email,
          End = auth$expiry,
          Datasets = auth$options$MAX_DATASETS
        )
        values[which(values == "")] <- "(not set)"
        data.frame(" " = names(values), "  " = values, check.names = FALSE)
      },
      width = "100%",
      striped = TRUE
    )

    render.baseplot <- function() {
      x <- unlist(nav_count())
      x <- x[grep("welcome|userprofile|setting", names(x), invert = TRUE)]
      x <- x[order(names(x))]
      score <- 1 * (length(x) > 12) + 1 * (sum(x) > 30) + 1 * (sum(x) > 100) + 1 * (max(x) > 20)
      level <- c("youngling", "padawan", "knight", "master", "grand master")[1 + score]
      if (length(x) < 20) x <- head(c(x, rep(0, 20)), 20)
      names(x) <- gsub("-tab", "", names(x))
      par(mar = c(8, 4, 2, 1))
      barplot(x, las = 3, ylim = c(0, max(3, 1.2 * max(x))), ylab = "module visits")
      title(paste("Your level: ", level), cex.main = 1)
    }

    render.plotly <- function() {
      x <- unlist(nav_count())
      x <- x[grep("welcome|userprofile|setting", names(x), invert = TRUE)]
      x <- x[order(names(x))]
      score <- 1 * (length(x) > 12) + 1 * (sum(x) > 30) + 1 * (sum(x) > 100) + 1 * (max(x) > 20)
      level <- c("Youngling", "Padawan", "Knight", "Master", "Grand Master")[1 + score]
      title <- paste("Your level is ", level)
      if (score > 2) title <- paste("Congratulations!", title)

      if (length(x) < 20) x <- head(c(x, rep(0, 20)), 20)
      names(x) <- gsub("-tab", "", names(x))
      df <- data.frame(x = names(x), y = x)
      ymax <- c(0, max(3, 1.2 * max(x)))

      plotly::plot_ly(
        data = df,
        y = ~y,
        x = ~x,
        type = "bar",
        ## orientation = "h",
        ## color = ~color, ## TODO: use variable that encodes grouping
        colors = omics_pal_d()(length(unique(df$color))),
        hovertemplate = "%{y}: %{x}<extra></extra>"
      ) %>%
        plotly_default() %>%
        plotly::layout(
          title = title,
          xaxis = list(title = FALSE),
          yaxis = list(
            title = "module visits",
            range = list(0, ymax)
          ),
          font = list(family = "Lato"),
          showlegend = FALSE,
          bargap = .4,
          margin = list(l = 5, r = 10, b = 10, t = 25)
        )
    }


    PlotModuleServer(
      "usage",
      func = render.plotly,
      plotlib = "plotly",
      res = 85,
      add.watermark = FALSE
    )
  })
}
