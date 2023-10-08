start_inactivityControl <- function(session, delta, inactivityCounter) {
    reactive({
        invalidateLater(delta * 1000)
        # Add 1 to inactivity counter
        shiny::isolate(inactivityCounter(inactivityCounter() + 1))
        # If >30 min inactivity, close session
        if(inactivityCounter() == 6){
        dbg("[SERVER:userLogout] >>> USER session ended due to inactivity")
        sever::sever(sever_ciao(), bg_color = "#004c7d")
        session$close()
        }
    })
}