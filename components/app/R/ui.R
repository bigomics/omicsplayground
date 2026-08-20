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

    ## Cache-bust on mtime. These responses carry no Cache-Control and no ETag,
    ## so browsers fall back to heuristic caching and keep a stylesheet or a
    ## script from an earlier deploy - which shows up as the app rendering with
    ## half its CSS, or a JS fix that appears not to have landed. Both "custom"
    ## and "static" are addResourcePath'd to www, so one mtime lookup serves.
    bust <- function(path) {
      f <- file.path("www", basename(path))
      if (!file.exists(f)) return(path)
      paste0(path, "?v=", as.integer(file.mtime(f)))
    }

    header <- shiny::tagList(
      shiny::tags$head(htmltools::includeHTML("www/hubspot-embed.html")),
      ##    gtag2, ## Google Tag Manager???
      shiny::tags$head(shiny::tags$script(src = bust("custom/temp.js"))),
      shiny::tags$head(shiny::tags$script(src = bust("static/copy-info-helper.js"))),
      shiny::tags$script(src = bust("custom/close-message.js")),
      shiny::tags$head(shiny::tags$script(src = bust("static/add-tick-helper.js"))),
      shiny::tags$head(shiny::tags$script(src = bust("static/shared-badges.js"))),
      shiny::tags$head(shiny::tags$script(src = bust("custom/dropdown-helper.js"))),
      shiny::tags$head(shiny::tags$link(rel = "stylesheet",
                                        href = bust("custom/styles.min.css"))),
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
    
    omicspanel <- function(p) {
      div(p, style="margin: 0 8px 0 8px;", class = "omicspanel")
    }  
    
    ##ui <- bslib::page_fillable(
    ui <- bigdash::bigPage(
      header,
      navbar = NULL,
      bslib::navset_pill_list(
        id = "app-sidebar",
        ##widths = c("50px","calc(100% - 50px)"),
        widths = c(1,11),
        selected = ifelse( isTRUE(opt$DEVMODE), "Apps", "Home"),
        well = TRUE,
        bslib::nav_panel(
          title = "Home",
          icon = icon("home"),
          omicspanel(WelcomeBoardUI("welcome2"))
        ),
        if(isTRUE(opt$DEVMODE)) {
          bslib::nav_panel(
            title = "Apps",
            icon = icon("app-store-ios", style="font-size: 38px;"),
            launcher_ui("apps")
          )
        },        
        bslib::nav_panel(
          title = "Library",
          icon=icon("book"),
          omicspanel(LoadingUI("load"))
        ),
        bslib::nav_panel(
          title = "Dashboard",
          icon = icon("chart-line"),
          opg_ui()
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
            title = shiny::tagList(icon("robot"), tags$br(), "Obi AI"),
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
        ## Hidden panels (e.g. tools)
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
        bslib::nav_panel_hidden("Prism",
          omicspanel(prism_ui("prism"))
        ),
        bslib::nav_panel_hidden("IDconvert",
          omicspanel(idconvert_ui("idconvert"))
        ),
        bslib::nav_panel_hidden("Qsee",
          omicspanel(qsee_ui("qsee"))
        ),
        ## Not DEVMODE-gated: a methylation dataset bypasses the Dashboard and
        ## lands here (see opg_server.R), so this panel has to exist on every
        ## build, not only the ones that show the Apps launcher.
        bslib::nav_panel_hidden("Methylome",
          omicspanel(methylome_ui("methylome"))
        ),
        bslib::nav_panel_hidden(
          value = "AcrossDatasets",
          #title = HTML("Across&nbsp;datasets"),
          #icon = icon("layer-group"),
          omicspanel(AcrossUI("across"))
        ),
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
      ),
      ## The "Dashboard" nav_panel above embeds opg_ui(), which builds its
      ## own bigdash::bigPage() (with its own sidebar/settings drawer) as a
      ## tab of this outer bigdash::bigPage(). Both default to the same
      ## "app" id, so both end up as `.bigdash-app` roots with
      ## id/data-bigdash-id="app" -- an invalid duplicate DOM id. bigdash's
      ## own JS (settings.js) iterates every `.bigdash-app` root once to
      ## bind the settings-lock click handler, so with two roots resolving
      ## to the same id it binds the handler *twice* on the same lock icon;
      ## the two bound handlers toggle the locked/unlocked state back and
      ## forth on every click, so the lock appears completely unresponsive
      ## (this is the reported bug -- see qsee_ui.R/qsee_server.R for a
      ## module that avoids this by using a distinct, fully-scoped id).
      ## This outer shell never uses its own bigdash sidebar/settings (no
      ## `sidebar`/`settings` args are passed above), so it's safe to strip
      ## its `.bigdash-app` marker before bigdash's init code runs, leaving
      ## opg_ui()'s nested page as the sole "app"-id root. Runs as a plain
      ## synchronous <script> (not $(document).ready(...)) so it executes
      ## while the page is still parsing, before bigdash's own
      ## `$(function(){ ... })` init handlers fire.
      tags$script(HTML(
        "(function() {
          var roots = document.querySelectorAll('.bigdash-app');
          if (roots.length > 1) {
            roots[0].classList.remove('bigdash-app');
            roots[0].removeAttribute('data-bigdash-id');
          }
        })();"
      ))
    )

    return(ui)
  }
}
