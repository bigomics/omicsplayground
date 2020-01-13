

alertDataLoaded <- function(session, ngs) {
    if(!is.null(ngs)) return()
    sendSweetAlert(
        session = session,
        ##title = "No dataset loaded",
        title = NULL,
        text = "Please first load a dataset"
    )
}
