

alertDataLoaded <- function(session, ngs) {
    if(!is.null(ngs)) return()
    sendSweetAlert(
        session = session,
        ##title = "No dataset loaded",
        title = NULL,
        text = "Please first load a dataset"
    )
}

pgx.randomCartoon <- function() {
    cartoon_list <- list(
        list(slogan="Visual analytics. See and understand", img="data-graph-wisdom.jpg"),
        list(slogan="Fasten your seat belts. Accelerated discovery", img="cartoon-speedup.jpg"),
        list(slogan="Analytics anywhere. Anytime.", img="cartoon-cloudservice.jpg"),
        ## list(slogan="Do-it-yourself. Yes you can.", img="bigomics-rockstar3.jpg"),
        list(slogan="Analyze with confidence. Be a rockstar", img="bigomics-rockstar3.jpg"),
        list(slogan="Fast track your Bioinformatics", img="selfservice-checkout2.png"),
        list(slogan="Integrate more. Dig deeper", img="cartoon-integration.jpg"),
        list(slogan="Your analysis doesn't take coffee breaks", img="gone-for-coffee.png"),
        list(slogan="Too much data? Help yourself", img="cartoon-datahelp2.jpg"),    
        list(slogan="Big Friendly Omics", img="big-friendly-omics1.jpg"),
        list(slogan="Big Data meets Biology", img="bigdata-meets.png")
    )
    ##randomCartoon <- reactive({
    ##invalidateLater(20000)
    cartoon <- sample(cartoon_list,1)[[1]]
    cartoon$img2 = file.path("cartoons",cartoon$img)
    cartoon$img  = file.path("www/cartoons",cartoon$img)
    cartoon
}
