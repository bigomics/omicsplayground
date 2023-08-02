##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


sendErrorMessageToCustomerSuport <- function(user, pgx_name, error, path_to_creds = "gmail_creds") {
  if (!file.exists(path_to_creds)) {
    info("[sendShareMessage] WARNING : mail not sent. cannot get mail creds =", path_to_creds)
    return(NULL)
  }

  user <- trimws(user)

  blastula::smtp_send(
    blastula::compose_email(
      body = blastula::md(
        glue::glue(
          "Hello, The user {user} had a dataset that failed to compute. ",
          "Please find here the log and data path.",
          "The error is : {error}",
          "The ds name is : {pgx_name}"
        )
      ),
      footer = blastula::md(
        glue::glue("Email sent on {blastula::add_readable_time()}.")
      )
    ),
    from = "bigomics.app@gmail.com",
    to = "mauro.masiero@bigomics.ch",
    subject = paste("Problem with dataset from", user),
    credentials = blastula::creds_file(path_to_creds)
  )
}

