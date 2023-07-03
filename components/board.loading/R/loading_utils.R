
##-------------------------------------------------------------------
## util functions
##-------------------------------------------------------------------

sendShareMessage <- function(pgxname, sender, share_user, path_to_creds='gmail_creds') {
  
  dbg("[sendShareMessage] pgxname = ", pgxname)
  dbg("[sendShareMessage] sender = ", sender)
  dbg("[sendShareMessage] share_user = ", share_user)
  
  if(!file.exists(path_to_creds)) return(NULL)
  
  blastula::smtp_send(
    blastula::compose_email(
      body = blastula::md(
        glue::glue(
          "Hello, {sender} shared a dataset with you on OmicsPlayground! ",
          "Login to accept the new dataset.")
      ),
      footer = blastula::md(
        glue::glue("Email sent on {blastula::add_readable_time()}.")
      )
    ),
    from = "app@bigomics.ch",
    to = share_user,
    subject = paste("Your friend",sender,"shared data on OmicsPlayground"),
    credentials = blastula::creds_file(path_to_creds)
  )
}


get_coworkers <- function(pgxdir, email) {
  domain <- sub(".*@","",email)
  if(email=="" || domain=="") return(NULL)
  cow <- dir(pgxdir, pattern=paste0(domain,"$"))
  cow <- setdiff(cow, email)
  cow
}

is_valid_email <- function(email) {
  is_personal <- grepl("gmail|ymail|outlook|yahoo|hotmail|mail.com$|icloud|msn.com$",email)
  valid_email <- grepl(".*@.*[.].*",email)
  valid_email <- valid_email && !grepl("[*/\\}{]",email) ## no special chars
  return(!is_personal && valid_email)
}

