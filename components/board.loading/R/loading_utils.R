## -------------------------------------------------------------------
## util functions
## -------------------------------------------------------------------

sendShareMessage <- function(pgxname, sender, share_user, path_to_creds = "gmail_creds") {
  if (!file.exists(path_to_creds)) {
    info("[sendShareMessage] WARNING : mail not sent. cannot get mail creds =", path_to_creds)
    return(NULL)
  }

  share_user <- trimws(share_user)
  sender <- trimws(sender)

  blastula::smtp_send(
    blastula::compose_email(
      body = blastula::md(
        glue::glue(
          "Hello,

          The user <strong>{sender}</strong> shared a dataset with you on Omics Playground.
          
          If you are a new user, please register here: https://auth.bigomics.ch/#!/register

          If you already use Omics Playground, simply login here: https://auth.bigomics.ch/#!/login

          You can find the shared dataset in the 'Load dataset/Sharing' tab.

          Yours,
          BigOmics Team"
        )
      ),
      footer = blastula::md(
        glue::glue("Email sent on {blastula::add_readable_time()}.")
      )
    ),
    from = "bigomics.app@gmail.com",
    to = share_user,
    subject = paste("Omics Playground: You received a new dataset from", sender),
    credentials = blastula::creds_file(path_to_creds)
  )
}


get_coworkers <- function(pgxdir, email) {
  domain <- sub(".*@", "", email)
  if (email == "" || domain == "") {
    return(NULL)
  }
  cow <- dir(pgxdir, pattern = paste0(domain, "$"))
  cow <- setdiff(cow, email)
  cow
}

is_valid_email <- function(email) {
  is_personal <- grepl("gmail|ymail|outlook|yahoo|hotmail|mail.com$|icloud|msn.com$", email)
  valid_email <- grepl(".*@.*[.].*", email)
  valid_email <- valid_email && !grepl("[*/\\}{]", email) ## no special chars
  return(!is_personal && valid_email)
}

makebuttonInputs2 <- function(FUN, len, id, tooltip = NULL, ...) {
  inputs <- character(length(len))
  for (i in seq_along(len)) {
    if (is.null(tooltip)) {
      inputs[i] <- as.character(FUN(paste0(id, len[i]), ...))
    } else {
      inputs[i] <- as.character(
        withTooltip(FUN(paste0(id, len[i]), ...), tooltip)
      )
    }
  }
  inputs
}
