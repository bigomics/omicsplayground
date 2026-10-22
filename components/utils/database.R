connect_db <- function(database_path) {
  connection <- DBI::dbConnect(RSQLite::SQLite(), dbname = database_path)
  return(connection)
}

disconnect_db <- function(connection) {
  DBI::dbDisconnect(connection)
}

query_by_email <- function(email, connection) {
  query_result <- DBI::dbGetQuery(connection, paste0("
    SELECT *
    FROM users
    WHERE email = '", email, "'
  "))
  if (nrow(query_result) == 0) {
    return(NULL)
  } else {
    return(query_result)
  }
}
