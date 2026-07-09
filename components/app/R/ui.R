##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
      shiny::tags$head(htmltools::includeHTML("www/hubspot-embed.html")),
      ##    gtag2, ## Google Tag Manager???
      shiny::tags$head(shiny::tags$script(src = "custom/temp.js")),
      shiny::tags$head(shiny::tags$script(src = "static/copy-info-helper.js")),
      shiny::tags$script(src = "custom/close-message.js"),
      shiny::tags$head(shiny::tags$script(src = "static/add-tick-helper.js")),
      shiny::tags$head(shiny::tags$script(src = "static/shared-badges.js")),
      shiny::tags$head(shiny::tags$script(src = "custom/dropdown-helper.js")),
      shiny::tags$head(shiny::tags$link(rel = "stylesheet", href = "custom/styles.min.css")),
      shiny::tags$head(shiny::tags$link(rel = "shortcut icon", href = "custom/favicon.ico")),
      visnetwork = visNetwork::visNetworkOutput("a", height = "0px"),
      shinyjs::useShinyjs(),
      waiter::use_waiter(),
      sever::useSever(),
      bigLoaders::addBigLoaderDeps(),
      firebase::useFirebase(firestore = TRUE, analytics = TRUE),
      shinybrowser::detect(),
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
    
    ## new multi-app UI
    nav_page <- function(...) {
      bslib::page_fluid(
        theme = bigdash::big_theme(),
        title = NULL,
        style = "padding: 0px;",
        bigdash::navbar(
          title = tags$img(
            src = "assets/img/bigomics.png",
            width = "110"
          ),
          center = tags$div(
              "title in navbar",
              style = 'text-align:center;width: 100%;'
          ),
          left = NULL,
          NULL
        ),
        ... 
      )
    }
    
    nav_page <- function(p) {p}  ## dummy
    
    ui <- bigdash::bigPage(
      header,
      navbar = NULL,
      bslib::navset_pill_list(
        id = "app-sidebar",
        ##widths = c("50px","calc(100% - 50px)"),
        widths = c(1,11),
        selected = "Home",
        well = TRUE,
        bslib::nav_panel(
          title = "Home",
          icon=icon("home"),
          nav_page(WelcomeBoardUI("welcome2"))
        ),
        bslib::nav_panel(
          title = "Library",
          icon=icon("book"),
          nav_page(
            div(LoadingUI("load"), class = "px-4 py-0")
          )
        ),
        bslib::nav_panel(
          title = "Dashboard",
          icon = icon("chart-line"),
          opg_ui()
        ),
        if (isTRUE(opt$ENABLE_ACROSS)) {
          bslib::nav_panel(
            title = HTML("Across&nbsp;datasets"),
            value = "AcrossDatasets",
            icon = icon("layer-group"),
            div(AcrossUI("across"), class = "px-4 py-0")
          )
        },
        bslib::nav_panel(
          title = HTML("AI&nbsp;Studio"),
          value="Studio",
          icon = icon("clapperboard"),
          div(StudioUI("studio"), class = "px-4 py-0")          
        ),
        if (copilot_packages_ok()) {
          bslib::nav_panel(
            #title = HTML("AI&nbsp;Copilot"),
            title = tagList(icon("robot"), tags$br(), "Obi"),
            value = "Copilot",
            div(CopilotBoardUI("copilot2"), class = "px-4 py-0")
          )
        },
        if(isTRUE(opt$DEVMODE)) {
          bslib::nav_panel(
            title = "Runs", icon=icon("person-running"),
            div( class = "px-4 py-0",
              ##shiny::div(id = "navheader-current-section", HTML("Runs")),
              ##p("Monitor and inspect the details of computation runs"),
              RunMonitorUI("runmonitor")
            )
          )
        },
        if(isTRUE(opt$DEVMODE)) {
          bslib::nav_panel(title = "Tools", icon = icon("tools"),
            tools_ui("tools")
          )
        },
        ## Hidden panels (e.g. tools)
        bslib::nav_panel_hidden("Prism",
          div(prism_ui("prism"), class='px-4 py-0')
        ),
        bslib::nav_panel_hidden("Upload",
          div(UploadUI("upload"), class='px-4 py-0')           
        ),
        bslib::nav_panel_hidden("UserProfile",
          div(UserProfileUI("user_profile"), class='px-4 py-0')
        ),
        if (isTRUE(opt$ENABLE_ADMIN)) {
          bslib::nav_panel_hidden("UserProfile",
            div(AdminPanelUI("admin_panel"), class='px-4 py-0')
          )
        },
        
        ## lower settings buttons
        bslib::nav_spacer(),
        bslib::nav_panel("Settings", icon=icon("cog"),
          div(AppSettingsUI("app_settings"), class='px-4 py-0') 
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

    return(ui)
  }
}
