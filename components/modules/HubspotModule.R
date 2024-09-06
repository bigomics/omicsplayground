hubspot_post <- function(url, body, access_token, type = "POST") {
    auth <- hubspot_auth_header(access_token)
    if (type == "POST") {
      response <- httr::POST(
          url,
          body = body,
          httr::content_type("application/json"),
          auth
      )
    } else if (type == "PATCH") {
      response <- httr::PATCH(
          url,
          body = body,
          httr::content_type("application/json"),
          auth
      )
    }
    return(response)
}

create_hubspot_search <- function(email) {
    body <- list(
    filterGroups = list(
      list(
        filters = list(
          list(
            value = tolower(email),
            propertyName = "email",
            operator = "EQ"
          )
        )
      )
    ),
    properties = list("email", "firstname", "lastname", "jobtitle", "background", "bioinformaticians_in_the_team")
  ) |> jsonlite::toJSON(auto_unbox = TRUE)
}

hubspot_auth_header <- function(access_token) {
  httr::add_headers(Authorization = paste0("Bearer ", access_token))
}
