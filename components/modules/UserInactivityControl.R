start_inactivityControl <- function(session, timeout, inactivityCounter) {
  reactive({
    invalidateLater(timeout / 6 * 1000)
    ia_counts <- isolate(inactivityCounter())
    shiny::isolate(inactivityCounter(ia_counts + 1)) ## increase counter
    ## If >30 min inactivity, close session
    if (ia_counts == 6) {
      dbg("[SERVER:userLogout] >>> USER session ended due to inactivity")
      sever::sever(sever_ciao(), bg_color = "#004c7d")
      session$close()
    }
  })
}
