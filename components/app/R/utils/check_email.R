##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

check_personal_email <- function(auth, pgx_dir, title=NULL, text=NULL) {
  email <- auth$email()
  is_personal_email <- grepl("gmail|ymail|outlook|yahoo|mail.com$|icloud",email)
  existing_user_dirs <- basename(list.dirs(pgx_dir))
  user_exists <- (email %in% existing_user_dirs)

  if(is.null(text)) {
    text <- "You are using a personal email adress. Please provide your business, academic or institutional email. Your login name will be changed and your data will be copied. After that please login again with your new email. Are you sure?"
  }
  if(is.null(title)) {
    title <- "Please change your email"
  }
  
  if(is_personal_email && user_exists) {
    shinyalert::shinyalert(
      inputId = "new_email",                        
      title = title,
      text = text,
      type = "input",
      callbackR = function(new_email) {
        old_dir_exists <- (email %in% existing_user_dirs)
        new_dir_exists <- (new_email %in% existing_user_dirs)
        ## copy old data to new data
        if(old_dir_exists && !new_dir_exists) {
          dbg("@@@@@ dry run @@@@@ copying data from",email,"to",new_email)
          old_dir <- file.path( pgx_dir, email)
          new_dir <- file.path( pgx_dir, new_email)
          base::file.rename(old_dir, new_dir)
          
          shinyalert::shinyalert(title="", text="Your login name has been changed and your data have been moved. Please login again with your new email.")          
          shinyjs::runjs("logoutInApp()")
        } else if(!old_dir_exists && !new_dir_exists) {
          shinyalert::shinyalert(title="", text="Your login name has been changed. Please login again with your new email.")          
          shinyjs::runjs("logoutInApp()")          
        } else if(new_dir_exists) {
          dbg("@@@@@ dry run @@@@@ ERROR!!",new_email,"exists")
          title = "Email already exists"
          text <- "This email already exists. Please provide a different business, academic or institutional email."
          check_personal_email(auth, pgx_dir, title=title, text=text)
        }
      }
    )
  } ## end if user_exists

}
