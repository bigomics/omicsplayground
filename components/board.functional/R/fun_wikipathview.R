##    devtools::install_github("m-jahn/fluctuator")
library(fluctuator)

wikipathview <- function(wp, val, dir) {
    
  require(fluctuator)
  get_pathway_svg <- function(wp) {
    fn <- grep(paste0("_",wp,"_"), dir(dir, pattern="svg$",
      full.names=TRUE), ignore.case=TRUE, value=TRUE)
    if(length(fn)) {
        svg <- fluctuator::read_svg(fn)
    } else {
        svg <- NULL
    }
    svg
  }  
  get_labels <- function(svg) {
    nodes <- fluctuator::get_attributes(svg, node="text", node_attr="node_set", attr=NULL)
    text <- sapply(nodes$id, function(id)
        fluctuator::get_values(svg, node=id, node_attr="id")) ## SLOW!
    names(text) <- nodes$id
    text
  }
  color_all_nodes <- function(svg, node_ids, colors) {
    aa  <- fluctuator::get_attributes(svg, node=node_ids, node_attr="id")  
    path_url <- unlist(aa['clip-path'])
    path_ids <- rep(NA,length(path_url))
    i=1
    for(i in 1:length(path_url)) {
        pp <- fluctuator::get_attributes(svg, node=path_url[i], node_attr="clip-path")  
        j <- which(pp$node_set=="path" & grepl("stroke:black",pp$style))
        if(length(j)==0) next()
        path_ids[i] <- pp$id[j]
    }
    
    fluctuator::set_attributes(
        svg, node=path_ids, node_attr='id', attr = "style",
        pattern = "fill:none", replacement = paste0("fill:",colors))
  }
    
  svg <- get_pathway_svg(wp)
  if(is.null(svg)) {
    return(NULL)
  }
    
  labels <- get_labels(svg)  ## slow...
  if(is.null(val)) {
    node_ids <- names(labels)        
    svg2 <- color_all_nodes(svg, node_ids, "#00ff0033")  
  } else {
    sum(names(val) %in% labels)
    if(sum(names(val) %in% labels)==0) {
      return(NULL)
    }
    labels <- labels[labels %in% names(val)]
    labels
    val <- val[labels]
    val
    node_ids <- names(labels)
    rr <- round(66*pmin(1,abs(val/2.0))**0.5)  
    colors <- ifelse(val>0, paste0("#ff0000",rr), paste0("#0055ff",rr))
    svg2 <- color_all_nodes(svg, node_ids, colors)  
  }
  svg2
}

if(0) {
    wp="WP558"
    wp="WP179"
    svg2 <- wikipathview(wp="WP179", val=NULL, dir="svg2")
    fluctuator::write_svg(svg2, file = "WP179-green.svg")
}


