##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


read_user_options <- function(user_dir) {
  user_opt_file <- file.path(user_dir, "OPTIONS")
  new_opt <- opt ## opt from global
  if (!file.exists(user_opt_file)) {
    ## IK: no need
    ## file.copy(from = opt.file, to = user_opt_file)
  } else {
    user_opt <- playbase::pgx.readOptions(file = user_opt_file)
    for (opt_name in names(user_opt)) {
      new_opt[[opt_name]] <- user_opt[[opt_name]]
    }
  }
  new_opt
}

upgrade.dialog <- function(ns, current.plan) {
  btn_basic <- "Go Basic!"
  btn_starter <- "Get Starter!"
  btn_premium <- "Get Premium!"
  if (current.plan == "free") btn_basic <- "Current Plan"
  if (current.plan == "starter") btn_starter <- "Current Plan"
  if (current.plan == "premium") btn_premium <- "Current Plan"

  modalDialog(
    title = h3("Find the right OmicsPlayground plan for you"),
    size = "m",
    div(
      class = "row",
      style = "padding-left:4rem;padding-right:4rem;text-align:center;",
      div(
        class = "col-md-4",
        style = "background:#F2FAFF;",
        HTML("<h4><b>Basic</b></h4>"),
        p("Try for free"),
        h3("Free!"),
        tags$ul(
          class = "list-unstyled",
          tags$li("Host up to 3 datasets"),
          tags$li("45 minutes time limit"),
          tags$li("Up to 25 samples / dataset"),
          tags$li("Up to 5 comparisons")
        ),
        shiny::actionButton(ns("get_basic"), btn_basic),
        br()
      ),
      div(
        class = "col-md-4",
        style = "background:#E8F8FF;",
        h4(HTML("<b>Starter</b>")),
        p("Great to start"),
        h3("Soon!"),
        tags$ul(
          class = "list-unstyled",
          tags$li("Host up to 10 datasets"),
          tags$li("3 hours time limit"),
          tags$li("Up to 100 samples / dataset"),
          tags$li("Up to 10 comparisons")
        ),
        shiny::actionButton(ns("get_starter"), btn_starter),
        br()
      ),
      div(
        class = "col-md-4",
        style = "background:#E2F4FF;",
        HTML("<h4><b>Premium</b></h4>"),
        p("For power users or small groups"),
        h3("Soon!"),
        tags$ul(
          class = "list-unstyled",
          tags$li("Host up to 100 datasets"),
          tags$li("8 hours time limit"),
          tags$li("Up to 2000 samples / dataset"),
          tags$li("Up to 100 comparisons")
        ),
        shiny::actionButton(ns("get_premium"), btn_premium),
        br()
      )
    ), ## content div
    div(
      style = "margin-top:3rem;text-align:center;",
      HTML("Looking for OmicsPlayground for <b>Enterprise</b>? <a href='mailto:info@bigomics.com'>Contact sales for info and pricing</a>.")
    ),
    footer = tagList(
      fillRow(
        flex = c(NA, 0.03, NA, 1, NA, NA),
        tags$label(
          class = "radio-inline",
          tags$input(
            id = "yearlyCheck",
            type = "radio",
            name = "yearly",
            onclick = "priceChange(name)",
            checked = TRUE
          ),
          "Billed yearly"
        ),
        br(),
        tags$label(
          class = "radio-inline",
          tags$input(
            id = "monthlyCheck",
            type = "radio",
            name = "monthly",
            onclick = "priceChange(name)"
          ),
          "Billed monthly"
        ),
        br(),
        shiny::actionButton(ns("manage"), "Manage Subscription"),
        modalButton("Dismiss")
      )
    )
  ) ## modalDialog
}

js.emailFeedbackMessage <- function(session, msg, type = "error") {
  session$sendCustomMessage(
    "email-feedback",
    list(
      type = type,
      msg = msg
    )
  )
}

checkAuthorizedDomain <- function(email, domain) {
  if (is.null(domain) || domain == "" || domain == "*") {
    return(TRUE)
  }
  domain1 <- strsplit(domain, split = "\\|")[[1]]
  domain1 <- paste0(paste0("@", domain1, "$"), collapse = "|")
  authorized <- grepl(domain1, email)
  authorized
}

checkAuthorizedUser <- function(email, credentials_file = NULL) {
  if (is.null(credentials_file) || credentials_file == FALSE) {
    return(TRUE)
  }
  if (!file.exists(credentials_file)) {
    return(TRUE)
  }
  CREDENTIALS <- read.csv(credentials_file, colClasses = "character")
  valid_user <- email %in% CREDENTIALS$email
  if (!valid_user) {
    return(FALSE)
  }
  sel <- match(email, CREDENTIALS$email)
  valid_date <- as.Date(CREDENTIALS$expiry[sel]) > as.Date(Sys.time())
  authorized <- valid_user && valid_date
  authorized
}

checkValidEmailFormat <- function(email) {
  if(is.null(email)) return(FALSE)
  if(email %in% c(""," ",NA)) return(FALSE)  
  grepl(".*@.*[.].*", email)
}

checkPersonalEmail <- function(email) {
  grepl("gmail|ymail|outlook|yahoo|hotmail|mail.com$|icloud|msn.com$", email)
}

checkEmail <- function(email, domain = NULL, credentials_file = NULL, check.personal = TRUE) {
  dbg("[AuthenticationModule:checkEmail]")
  chk <- list()
  if (!checkValidEmailFormat(email)) {
    return(list(valid = FALSE, msg = "not a valid email"))
  }
  if (!checkAuthorizedDomain(email, domain)) {
    return(list(valid = FALSE, msg = "domain not authorized"))
  }
  if (!checkAuthorizedUser(email, credentials_file)) {
    return(list(valid = FALSE, msg = "not authorized user"))
  }
  if (check.personal) {
    if (checkPersonalEmail(email)) {
      return(list(valid = FALSE, msg = "No personal email allowed. Please use your business, academic or institutional email."))
    }
  }
  list(valid = TRUE, "email ok")
}


