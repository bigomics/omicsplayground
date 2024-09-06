##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


create_user_dir_if_needed <- function(user_dir, pgxdir) {
  if (!dir.exists(user_dir)) {
    dir.create(user_dir)
    example_file <- file.path(pgxdir, "example-data.pgx")
    if (file.exists(example_file)) {
      file.copy(example_file, user_dir)
    }
  }
}

check_user_options_db <- function(email, user_database = NULL) {
  if (is.null(user_database)) {
    return(FALSE)
  }
  connection <- connect_db(user_database)
  user_opt <- query_by_email(email, connection)
  if (is.null(user_opt)) {
    return(FALSE) # user NOT in db
  } else {
    return(TRUE) # user IN db
  }
}

read_user_options <- function(user_dir) {
  user_opt_file <- file.path(user_dir, "OPTIONS")
  new_opt <- opt ## opt from global
  if (file.exists(user_opt_file) && !dir.exists(user_opt_file)) {
    user_opt <- playbase::pgx.readOptions(file = user_opt_file)
    ## restrict user options only to these options.
    ALLOWED_USER_OPTS <- c(
      "ENABLE_CHIRP", "ENABLE_DELETE", "ENABLE_PGX_DOWNLOAD",
      "ENABLE_PUBLIC_SHARE", "ENABLE_UPLOAD", "ENABLE_USER_SHARE",
      "MAX_DATASETS", "MAX_SAMPLES", "MAX_COMPARISONS",
      "MAX_GENES", "MAX_GENESETS", "MAX_SHARED_QUEUE",
      "TIMEOUT", "WATERMARK"
    )
    dbg("[read_user_options] 1 : names(user_opt) = ", names(user_opt))
    user_opt <- user_opt[which(names(user_opt) %in% ALLOWED_USER_OPTS)]
    dbg("[read_user_options] 2 : names(user_opt) = ", names(user_opt))
    for (opt_name in names(user_opt)) {
      new_opt[[opt_name]] <- user_opt[[opt_name]]
    }
  }
  # add user dir to opt file (IK: this is not an option!)
  ##  new_opt$user_dir <- user_dir
  new_opt
}

read_user_options_db <- function(email, user_database = NULL) {
  connection <- connect_db(user_database)
  user_opt <- query_by_email(email, connection)
  new_opt <- opt ## opt from global
  disconnect_db(connection)
  if (!is.null(user_opt)) {
    ## restrict user options only to these options.
    ALLOWED_USER_OPTS <- c(
      "ENABLE_CHIRP", "ENABLE_DELETE", "ENABLE_PGX_DOWNLOAD",
      "ENABLE_PUBLIC_SHARE", "ENABLE_UPLOAD", "ENABLE_USER_SHARE",
      "MAX_DATASETS", "MAX_SAMPLES", "MAX_COMPARISONS",
      "MAX_GENES", "MAX_GENESETS", "MAX_SHARED_QUEUE",
      "TIMEOUT", "WATERMARK"
    )
    dbg("[read_user_options] 1 : names(user_opt) = ", names(user_opt))
    user_opt <- user_opt[which(names(user_opt) %in% ALLOWED_USER_OPTS)]
    user_opt <- user_opt %>%
      dplyr::mutate(
        ENABLE_CHIRP = as.logical(ENABLE_CHIRP),
        ENABLE_DELETE = as.logical(ENABLE_DELETE),
        ENABLE_PGX_DOWNLOAD = as.logical(ENABLE_PGX_DOWNLOAD),
        ENABLE_PUBLIC_SHARE = as.logical(ENABLE_PUBLIC_SHARE),
        ENABLE_UPLOAD = as.logical(ENABLE_UPLOAD),
        ENABLE_USER_SHARE = as.logical(ENABLE_USER_SHARE),
        WATERMARK = as.logical(WATERMARK)
      )
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

## -----------------------------------------------------------------
## Email checking functions
## -----------------------------------------------------------------

checkAuthorizedDomain <- function(email, domain) {
  if (is.null(domain) || domain == "" || domain == "*") {
    return(TRUE)
  }
  domain1 <- strsplit(domain, split = "\\|")[[1]]
  domain1 <- paste0(paste0(domain1, "$"), collapse = "|")
  authorized <- grepl(tolower(domain1), tolower(email))
  authorized
}

checkBlockedDomain <- function(email, domain) {
  if (is.null(domain) || domain == "" || domain == "*") {
    return(FALSE)
  }
  domain1 <- strsplit(domain, split = "\\|")[[1]]
  domain1 <- paste0(paste0(domain1, "$"), collapse = "|")
  blocked <- grepl(tolower(domain1), tolower(email))
  blocked
}

checkAuthorizedUser <- function(email, credentials_file = NULL) {
  if (is.null(credentials_file) || credentials_file == FALSE) {
    return(TRUE)
  }
  if (!file.exists(credentials_file)) {
    return(TRUE)
  }
  CREDENTIALS <- read.csv(credentials_file, colClasses = "character")
  valid_user <- tolower(email) %in% tolower(CREDENTIALS$email)
  if (!valid_user) {
    return(FALSE)
  }
  sel <- match(tolower(email), tolower(CREDENTIALS$email))
  valid_date <- as.Date(CREDENTIALS$expiry[sel]) > as.Date(Sys.time())
  authorized <- valid_user && valid_date
  authorized
}

checkExpiredUser <- function(email, user_database) {
  if (is.null(user_database) || user_database == FALSE) {
    return(TRUE)
  }
  if (!file.exists(user_database)) {
    return(TRUE)
  }
  connection <- connect_db(user_database)
  # connection <- DBI::dbConnect(RSQLite::SQLite(), dbname = user_database)
  query_result <- DBI::dbGetQuery(connection, paste0("
    SELECT expiry
    FROM users
    WHERE email = '", email, "'
  "))
  # DBI::dbDisconnect(connection)
  disconnect_db(connection)
  if (nrow(query_result) == 0) {
    return(TRUE)
  } else {
    valid_date <- as.Date(query_result[1, 1]) > as.Date(Sys.time())
    return(valid_date)
  }
}

checkValidEmailFormat <- function(email) {
  if (is.null(email)) {
    return(FALSE)
  }
  if (email %in% c("", " ", NA)) {
    return(FALSE)
  }
  grepl(".*@.*[.].*", email)
}

checkPersonalEmail <- function(email) {
  grepl("gmail|ymail|outlook|yahoo|hotmail|mail.com$|icloud|msn.com$|qq.com", tolower(email))
}

checkMissingEmail <- function(email) {
  (is.null(email) || is.na(email) || email %in% c("", " ", NA))
}

## PGX.DIR="~/Playground/omicsplayground/data/"
checkExistUserFolder <- function(email) {
  user_dirs <- list.dirs(PGX.DIR, full.names = FALSE, recursive = FALSE)
  user_dirs <- grep("@", user_dirs, value = TRUE)
  tolower(email) %in% tolower(user_dirs)
}

checkEmail <- function(email, domain = NULL, blocked_domain = NULL, user_database = NULL,
                       check.personal = TRUE, check.existing = FALSE, credentials_file = NULL) {
  chk <- list()
  if (checkMissingEmail(email)) {
    return(list(valid = FALSE, msg = "missing email"))
  }
  if (!checkValidEmailFormat(email)) {
    return(list(valid = FALSE, msg = "not a valid email"))
  }
  if (!checkAuthorizedDomain(email, domain)) {
    return(list(valid = FALSE, msg = "domain not authorized"))
  }
  if (checkBlockedDomain(email, blocked_domain)) {
    return(list(valid = FALSE, msg = "domain blocked"))
  }
  if (!checkExpiredUser(email, user_database)) {
    return(list(valid = FALSE, msg = "user expired"))
  }
  if (check.personal) {
    if (checkPersonalEmail(email)) {
      return(list(valid = FALSE, msg = "No personal email allowed. Please use your business, academic or institutional email."))
    }
  }
  if (check.existing) {
    if (!checkExistUserFolder(email)) {
      return(list(valid = FALSE, msg = "username does not exist"))
    }
  }
  list(valid = TRUE, "email ok")
}

checkHubspot <- function(user_email) {
  dbg("[HubSpotCheckModule] checking hubspot for user data for email: ", user_email)
  search_body <- create_hubspot_search(user_email)
  credential <- file.path(ETC, "hubspot_creds")
  if (!file.exists(credential)) {
    stop("[HubSpotCheckModule] Error: Hubspot check is TRUE but no /etc/hubspot_creds file is on the deploy.")
  }
  HUBSPOT_SECRET <- readLines(credential)
  response <- hubspot_post(
      url = "https://api.hubapi.com/crm/v3/objects/contacts/search",
      body = search_body,
      access_token = HUBSPOT_SECRET
  )
  if (httr::content(response)$total == 1) {
    lastname <- httr::content(response)$results[[1]]$properties$lastname
    firstname <- httr::content(response)$results[[1]]$properties$firstname
    jobtitle <- httr::content(response)$results[[1]]$properties$jobtitle
    background <- httr::content(response)$results[[1]]$properties$background
    bioinfos_team <- httr::content(response)$results[[1]]$properties$bioinformaticians_in_the_team
    id <- httr::content(response)$results[[1]]$id
    is_null_or_empty <- function(x) {
        is.null(x) || (is.character(x) && x == "")
    }
    is_any_null_or_empty <- is_null_or_empty(lastname) |
            is_null_or_empty(firstname) |
            is_null_or_empty(jobtitle) |
            is_null_or_empty(background) |
            is_null_or_empty(bioinfos_team)
    time_since_update <- as.numeric(Sys.Date() - as.Date(httr::content(response)$results[[1]]$updatedAt))
    if (time_since_update < 1) {
      dbg("[HubSpotCheckModule] user is less than one day old, skipping data complete check")
      is_any_null_or_empty <- FALSE
    }
    if(is_any_null_or_empty) { # Missing info, redirect to auth
      dbg("[HubSpotCheckModule] missing data for ", user_email, " redirecting to auth")
      shinyjs::runjs(
        paste0(
          "window.location.replace('https://auth.bigomics.ch/#!/login?email=",
          user_email,
          "');"
        )
      )
    }
  } else {
    warning("[HubSpotCheckModule] Error: User ", user_email, " not on HubSpot")
  }
}

## ================================================================================
## ================================= END OF FILE ==================================
## ================================================================================
