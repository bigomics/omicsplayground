extract_cookie_value <- function(cookie_string, cookie_name) {
    cookies <- unlist(strsplit(cookie_string, "; "))
    user_cookie <- grep(paste0("^", cookie_name, "="), cookies, value = TRUE)
    if(length(user_cookie) > 0) {
        return(unlist(strsplit(user_cookie, "="))[2])
    } else {
        return(NULL)
    }
}
