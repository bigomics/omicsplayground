#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # list functions in global
  
  board = options()$board
  authentication = options()$authentication

  server_fn_name <- glue::glue("{board}board")
  board_server <- grep(server_fn_name, ls(envir = .GlobalEnv), value = TRUE, ignore.case = TRUE)
  length <- nchar(board) + nchar("board")
    
  board_server <- board_server[which(length == nchar(board_server))]

  board_server_fn <- get(board_server)

  # authentication

  if (!is.null(authentication) &&  authentication != "password") {
      credentials_file <- NULL
  }

  if (!is.null(authentication) && authentication == "password") {
      auth <- PasswordAuthenticationModule(
      id = "auth",
      credentials_file = credentials_file,
      allow_personal = opt$ALLOW_PERSONAL_EMAIL,
      domain = opt$DOMAIN
      )
  } else if (!is.null(authentication) && authentication == "firebase") {
      auth <- FirebaseAuthenticationModule(
      id = "auth",
      domain = opt$DOMAIN,
      firebase.rds = "firebase.rds",
      credentials_file = credentials_file,
      allow_personal = opt$ALLOW_PERSONAL_EMAIL,
      allow_new_users = opt$ALLOW_NEW_USERS
      )
  } else if (!is.null(authentication) && authentication == "email-link") {
      auth <- EmailLinkAuthenticationModule(
      id = "auth",
      pgx_dir = PGX.DIR,
      domain = opt$DOMAIN,
      firebase.rds = "firebase.rds",
      credentials_file = credentials_file,
      allow_personal = opt$ALLOW_PERSONAL_EMAIL,
      allow_new_users = opt$ALLOW_NEW_USERS
      )
  } else if (!is.null(authentication) && authentication == "login-code") {
      auth <- LoginCodeAuthenticationModule(
      id = "auth",
      mail_creds = file.path(ETC, "gmail_creds"),
      domain = opt$DOMAIN,
      user_database = user_database,
      allow_personal = opt$ALLOW_PERSONAL_EMAIL,
      allow_new_users = opt$ALLOW_NEW_USERS
      )
  } else if (!is.null(authentication) && authentication == "login-code-no-mail") {
      auth <- LoginCodeNoEmailAuthenticationModule(
      id = "auth",
      mail_creds = file.path(ETC, "gmail_creds"),
      domain = opt$DOMAIN,
      credentials_file = credentials_file,
      allow_personal = opt$ALLOW_PERSONAL_EMAIL,
      allow_new_users = opt$ALLOW_NEW_USERS
      )
  } else if (!is.null(authentication) && authentication == "shinyproxy") {
      username <- Sys.getenv("SHINYPROXY_USERNAME")
      auth <- NoAuthenticationModule(
      id = "auth",
      show_modal = TRUE,
      username = username,
      email = username
      )
  } else if (!is.null(authentication) && authentication %in% c("none2","none")) {
      ## no authentication but also not showing main modal (enter)
      auth <- NoAuthenticationModule(id = "auth", show_modal = FALSE)
  } else {
      auth <- NoAuthenticationModule(id = "auth", show_modal = FALSE)
  }

  trigger_server <- reactive({
        req(input$pgx_path)
        pgx <- playbase::pgx.load(input$pgx_path)
        pgx <- playbase::pgx.initialize(pgx)
        DATATYPEPGX <<- tolower(pgx$datatype)
        server <- board_server_fn(board, pgx = pgx)

  })
  
  observeEvent(input$pgx_path, {
    trigger_server()
  })
}
