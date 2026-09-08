##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

app_ui <- function(x) {
  if (identical("/cookie", x$PATH_INFO)) {
    value <- x$HTTP_HEADER_USER_COOKIE
    return(cookies::set_cookie_response(
      cookie_name = "persistentOPG",
      cookie_value = value,
      http_only = TRUE,
      secure_only = TRUE,
      redirect = "/close",
      same_site = "Strict",
      path = "/"
    ))
  } else if (identical("/cookie_nonce", x$PATH_INFO)) {
    value <- x$HTTP_HEADER_USER_COOKIE
    return(cookies::set_cookie_response(
      cookie_name = "persistentOPG_nonce",
      cookie_value = value,
      http_only = TRUE,
      secure_only = TRUE,
      redirect = "/close",
      same_site = "Strict",
      path = "/"
    ))
  } else if (identical("/cookie_remove", x$PATH_INFO)) {
    return(cookies::set_cookie_response(
      cookie_name = "persistentOPG",
      cookie_value = "",
      expiration = -1,
      http_only = TRUE,
      secure_only = TRUE,
      redirect = "/close",
      same_site = "Strict",
      path = "/"
    ))
  } else if (identical("/close", x$PATH_INFO)) {} else if (identical("/", x$PATH_INFO)) {

    theme <- bslib::bs_add_variables(bslib::bs_theme(),
      "grid-breakpoints" = # here e.g. with lg: 800px;
        "(xs: 0, sm: 576px, md: 768px, lg: 1200px, xl: 1600px, xxl: 2000px)",
      .where = "declarations"
    )

    header <- shiny::tagList(
      shiny::tags$head(htmltools::includeHTML(file.path(APPDIR,"assets/hubspot-embed.html"))),
      ##    gtag2, ## Google Tag Manager???
      shiny::tags$head(shiny::tags$script(src = "custom/temp.js")),
      shiny::tags$script(src = "custom/close-message.js"),
      shiny::tags$head(shiny::tags$script(src = "static/shared-badges.js")),
      shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "custom/styles.min.css")),
      shiny::tags$head(shiny::tags$link(rel = "shortcut icon", href = "custom/favicon.ico")),
      ## visnetwork must be inside a proper tagList, NOT named, otherwise
      ## all subsequent items (shinyalert, shinyjs, waiter, etc.) are
      ## silently dropped from the DOM (R's tagList treats named args as
      ## attributes, not children).
      shiny::tagList(
        shinyjs::useShinyjs(),
        shinyalert::useShinyalert(force = TRUE),
        waiter::use_waiter(),
        sever::useSever(),
        bigLoaders::addBigLoaderDeps(),
        shinybrowser::detect(),
        ## If you need to load visNetwork JS library, add:
        ## visNetwork::visNetworkOutput("a", height = "0px")
        ## BUT ensure it is NOT a named argument, or it will silently
        ## drop all subsequent items from the DOM (ask me how I know).
      ),
      shinybusy::busy_start_up(
        text = tags$h2("\nPrepping your personal playground..."), mode = "auto",
        background = "#2780e3", color = "#ffffff",
        loader = shinybusy::spin_epic("hollow-dots", color = "#FFF")
      )
    )

    nav_weblink <- function(title, href, onClick = NULL) {
      bslib::nav_item(NULL, shiny::tags$a(
        href = href, target = "_blank", onClick = onClick,
        HTML(paste(title,"<i class='fa-solid fa-arrow-up-right-from-square weblink' style='font-size: 13px;'></i>"))
      ))
    }

    nav_signout <- function(title, href, onClick = NULL) {
      bslib::nav_item(NULL, shiny::tags$a(
        href = href, target = "_blank", onClick = onClick,
        HTML(paste(title,"<i class='fa-solid fa-right-from-bracket weblink' style='font-size: 15px;'></i>"))
      ))
    }


    signout_link <- nav_signout(
      "Sign out", href = NULL, onClick = "logoutInApp(); setTimeout(() => window.location.reload(), 200);"
    )
    if (opt$AUTHENTICATION == "shinyproxy") {
      ## For ShinyProxy we need to redirect to /logout for clean session
      ## logout. Then we need a redirect to the /login page.
      signout_link <- nav_signout("Sign out", onClick = "shinyproxy_logout();", href = "/login")
    } else if (opt$AUTHENTICATION == "apache-cookie") {
      ## For apache SSO we need to redirect to /mellon/logout for SSO logout
      signout_link <- nav_signout("Sign out", onClick = NULL,
        href = paste0(opt$APACHE_COOKIE_PATH, "mellon/logout?ReturnTo=#"))
    }

    ## Plain bslib page rather than an outer bigdash::bigPage(): the
    ## "Dashboard" and "Qsee" nav_panels below embed opg_ui()/qsee_ui(),
    ## each their own bigdash::bigPage(). 
    ui <- shiny::bootstrapPage(
      title = "BigOmics",
      theme = bigdash::big_theme(),
      bigdash::dependencies(),
      div(
        id = "app-shell",
        class = "d-flex",
        style = "min-height:100vh;",
        div(
          class = "flex-grow-1 p-0 w-100",
          header,
          bslib::navset_pill_list(
            id = "app-sidebar",
            ##widths = c("50px","calc(100% - 50px)"),
            widths = c(1,11),
            selected = "Home",
            well = TRUE,
            ## bslib::nav_panel(
            ##   title = "Home",
            ##   icon = icon("home"),
            ##   omicspanel(WelcomeBoardUI("welcome2"))
            ## ),
            bslib::nav_panel(
              title = "Home",
              icon = icon("home"),
              #icon = icon("app-store-ios", style="font-size: 38px;"),
              launcher_ui("apps")
            ),
            bslib::nav_panel(
              title = "Library",
              icon=icon("book"),
              omicspanel(LoadingUI("load"))
            ),
            bslib::nav_panel(
              title = "Dashboard",
              icon = icon("chart-line"),
              ## The sidebar DOM is built once per session; the MultiOmics
              ## view's flat-boards layout (opg_promote_multiomics()) has to
              ## be baked in here, not left to runtime show/hide -- see
              ## opg_promote_multiomics()'s docstring.
              opg_ui("app", menu_tree = opg_promote_multiomics(opg_menu_tree()))
            ),
            ## AI tabs render only when the deployment licenses AI (opt$ENABLE_AI).
            ## The runtime "Enable AI" switch further shows/hides them per session
            ## via bigdash.toggleTab in appsettings_server.R.
            if (isTRUE(opt$ENABLE_AI)) {
              bslib::nav_panel(
                title = HTML("Studio"),
                value = "Studio",
                icon = icon("clapperboard"),
                omicspanel(StudioUI("studio"))
              )
            },
            if (isTRUE(opt$ENABLE_AI) && copilot_packages_ok()) {
              bslib::nav_panel(
                #title = HTML("AI&nbsp;Copilot"),
                title = shiny::tagList(icon("robot"), tags$br(), HTML("Obi&nbsp;AI")),
                value = "Copilot",
                omicspanel(CopilotBoardUI("copilot2"))
              )
            },
            if(isTRUE(opt$DEVMODE)) {
              bslib::nav_panel(
                title = "Runs", icon = icon("person-running"),
                omicspanel(RunMonitorUI("runmonitor"))
              )
            },
            ## Hidden panels
            bslib::nav_panel_hidden("Upload",
              omicspanel(UploadUI("upload"))
            ),
            bslib::nav_panel_hidden("UserProfile",
              omicspanel(UserProfileUI("user_profile"))
            ),
            if (isTRUE(opt$ENABLE_ADMIN)) {
              bslib::nav_panel_hidden("AdminPanel",
                omicspanel(AdminPanelUI("admin_panel"))
              )
            },
            ## Tools
            if(isTRUE(opt$DEVMODE)) {
              bslib::nav_panel_hidden("Prism",
                omicspanel(prism_ui("prism"))
              )
            },
            if(isTRUE(opt$DEVMODE)) {
              bslib::nav_panel_hidden("IDconvert",
                omicspanel(idconvert_ui("idconvert"))
              )
            },
            ## if(isTRUE(opt$DEVMODE)) {
            ##   bslib::nav_panel_hidden(
            ##     value = "AcrossDatasets",
            ##     omicspanel(AcrossUI("across"))
            ##   )
            ## },
            ## lower settings buttons
            bslib::nav_spacer(),
            bslib::nav_panel("Settings", icon=icon("cog"),
              omicspanel(AppSettingsUI("app_settings"))
            ),
            bslib::nav_menu(
              title = "Help",
              icon = icon("circle-question"),
              bslib::nav_item(NULL, actionLink("navbar_about", "About")),
              nav_weblink("Documentation", href="https://omicsplayground.readthedocs.io/"),
              nav_weblink("Video tutorials", href="https://bigomics.ch/tutorials/"),
              nav_weblink("Google forum", href="https://groups.google.com/d/forum/omicsplayground/"),
              nav_weblink("Reddit r/omicsplayground", href="https://www.reddit.com/r/omicsplayground"),
              nav_weblink("Submit a support ticket", href="https://share-eu1.hsforms.com/1glP7Cm6GQrWIGXgZrC0qrweva7t"),
              nav_weblink("Github issues", href="https://github.com/bigomics/omicsplayground/issues/"),
              nav_weblink("Case studies", href="https://bigomics.ch/case-studies/")
            ),
            bslib::nav_menu(
              title = "",
              icon = icon("user"),
              bslib::nav_item(NULL, actionLink("my_profile", "My profile")),
              if (isTRUE(opt$ENABLE_ADMIN)) {
                bslib::nav_item(NULL, actionLink("show_admin", "Admin panel"))
              },
              bslib::nav_item(NULL, InviteFriendUI("invite", type="link")),
              nav_weblink("Pricing &amp; Features", href="https://bigomics.ch/pricing/"),
              nav_weblink("Buy us coffee", href="https://buymeacoffee.com/bigomics"),
              signout_link
              )
          )
        )
      )
    )

    return(ui)
  }
}
