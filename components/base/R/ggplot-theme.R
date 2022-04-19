library(systemfonts)

#' The plot theme to be used for figures in the OmicsPlayground app.
#'
#' @param style (string) Overall color style of text labels.
#' Either "default" or "light".
#' @param base_size (integer) Base point size
#' @param grid (string) Grid lines. Options include  "none" or
#' any combination of "X", "Y", "x" and "y".
#' @param axistitle (string) Axis titles. Options include "none" or
#' any combination of "X", "Y", "x" and "y".
#' @param axistext (string) Axis text labels for values or groups.
#' Options include "none" or any combination of "X", "Y", "x" and "y".
#' @param axis_num (string) Should axis text be formatted as monospaced? Set 
#' to  x and y, respectively, in case numeric values are displayed. Options 
#' include "none" or any combination of "X", "Y", "x" and "y". 
#' @param legend_num (logical) Should legend text be formatted as monospaced?
#' Defaults to FALSE (no monospace font). Set to TRUE in case of numeric values.
#' #' @param panelborder (logical) Should a panel border be drawn?
#' Defaults to FALSE (no border). If TRUE it also adds tick marks to both axes.
#' @param margin (numeric) Should a margin of x be added to the plot?
#' Defaults to 0 (no margin by default).
#' @param ... Other arguments passed to ggplot methods.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mpg, aes(class)) + geom_bar() + theme_omics()
#' ggplot(mpg, aes(class)) + geom_bar() +
#'   theme_omics(style = "light", grid = "xy", margin = 20)
#' ggplot(mpg, aes(class)) + geom_bar() +
#'   theme_omics(grid = "none", axistext = "x", 
#'               axistitle = "none", panelborder = TRUE)
#' }
#'
#' @export
theme_omics <- function(style = "default", base_size = 15,
                        grid = "xy", axistitle = "xy", 
                        axistext = "xy", axis_num = "none",
                        legend_num = FALSE, panelborder = FALSE, 
                        margin = 0, ...) {
  if(!style %in% c("default", "light")) stop('style must be either "default" or "light"')
  if(!is.character(grid)) stop('grid must be a character: "none" or any combination of "X", "Y", "x" and "y"')
  if(!is.character(axistitle)) stop('axistitle must be a character: "none" or any combination of "X", "Y", "x" and "y"')
  if(!is.character(axistext)) stop('axistext must be a character: "none" or any combination of "X", "Y", "x" and "y"')
  if(!is.character(axis_num)) stop('axis_num must be a character: "none" or any combination of "X", "Y", "x" and "y"')
  if(!is.logical(legend_num)) stop('legend_num must be a logical variable')
  if(!is.logical(panelborder)) stop('panelborder must be a logical variable')
  if(!is.numeric(margin)) stop('margin must be a numeric value')
  
  fontfamily <- "" ## Clear Sans
  fontfamily_mono <- "" ## Fira Code
  
  if (style == "default") { base_col <- "black"; light_col <- "grey15" }
  if (style == "light") { base_col <- "grey30"; light_col <- "grey45" }
  
  out <-
    ggplot2::theme_minimal(base_size = base_size, base_family = fontfamily) +
    ggplot2::theme(
      text = ggplot2::element_text(
        color = base_col
      ),
      line = ggplot2::element_line(
        color = light_col
      ),
      rect = ggplot2::element_rect(
        color = light_col, 
        fill = "transparent"
      ),
      plot.title = ggtext::element_textbox_simple(
        size = base_size * 1.7,
        face = "bold",
        lineheight = .8, 
        box.color = NA,
        margin = ggplot2::margin(t = 0, b = base_size * .67)
      ),
      plot.subtitle = ggtext::element_textbox_simple(
        color = "grey40",
        size = base_size, 
        lineheight = 1.2,
        margin = ggplot2::margin(t = 0, b = base_size * 1.5)
      ),
      plot.caption = ggtext::element_textbox_simple(
        color = "grey40",
        size = base_size / 2, 
        lineheight = 1.2, 
        margin = ggplot2::margin(t = base_size * 1.5, b = 0)
      ),
      axis.title.x = ggplot2::element_text(
        size = base_size * .8,
        face = "bold",
        margin = ggplot2::margin(t = base_size / 3, r = 3, b = 3, l = 3)
      ),
      axis.title.y = ggplot2::element_text(
        size = base_size * .8,
        face = "bold",
        margin = ggplot2::margin(t = 3, r = 3, b = base_size / 3, l = 3)
      ),
      axis.title.y.left = ggplot2::element_text(
        size = base_size * .8,
        face = "bold",
        margin = ggplot2::margin(t = 3, r = 3, b = base_size / 3, l = 3)
      ),
      axis.title.y.right = ggplot2::element_text(
        size = base_size * .8,
        face = "bold",
        margin = ggplot2::margin(t = 3, r = 3, b = base_size / 3, l = 3)
      ),
      axis.text.x = ggplot2::element_text(
        color = light_col,
        size = base_size * .7,
        margin = ggplot2::margin(t = base_size / 4, r = 1, b = 1, l = 1)
      ),
      axis.text.y = ggplot2::element_text(
        color = light_col,
        size = base_size * .7,
        margin = ggplot2::margin(t = 1, r = base_size / 4, b = 1, l = 1)
      ),
      axis.ticks.length = grid::unit(.35, "lines"),
      strip.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "grey87", size = .25
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        color = "transparent", 
        fill = "transparent"
      ),
      plot.background = ggplot2::element_rect(
        color = "white", 
        fill = "white"
      ),
      plot.margin = ggplot2::margin(
        t = margin, r = margin, l = margin, b = margin * .67
      ),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.position = "right",
      legend.title = ggplot2::element_text(
        color = light_col,
        size = base_size * .8,
        face = "bold",
        margin = margin(b = 10)
      ),
      legend.text = ggplot2::element_text(
        color = light_col,
        size = base_size * .7
      )
    )
  
  if (grid != "none") {
    if (!stringr::str_detect(grid, "X|x")) {
      out <- out +
        ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
    }
    if (!stringr::str_detect(grid, "Y|y")) {
      out <- out +
        ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
    }
  } else {
    out <- out + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  }
  
  if (axistitle != "none") {
    if (!stringr::str_detect(axistitle, "X|x")) {
      out <- out +
        ggplot2::theme(axis.title.x = ggplot2::element_blank())
    }
    if (!stringr::str_detect(axistitle, "Y|y")) {
      out <- out +
        ggplot2::theme(axis.title.y = ggplot2::element_blank())
    }
  } else {
    out <- out +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank())
  }
  
  if (axis_num != "none") {
    if (stringr::str_detect(axis_num, "X|x")) {
      out <- out +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
          family = fontfamily_mono,
          color = light_col,
          size = base_size * .7,
          margin = ggplot2::margin(t = base_size / 4, r = 1, b = 1, l = 1)
        ))
    }
    if (stringr::str_detect(axis_num, "Y|y")) {
      out <- out +
        ggplot2::theme(axis.text.y = ggplot2::element_text(
          family = fontfamily_mono,
          color = light_col,
          size = base_size * .7,
          margin = ggplot2::margin(t = base_size / 4, r = 1, b = 1, l = 1)
        ))
    }
  }
  
  if (panelborder == TRUE) {
    out <- out +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(
          color = "grey70", 
          fill = "transparent",
          size = .5
        ),
        axis.ticks = ggplot2::element_line(
          size = .5, 
          color = "grey70"
        ) 
      )
  }
  
  if (axistext != "none") {
    if (!stringr::str_detect(axistext, "X|x")) {
      out <- out +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
    }
    if (!stringr::str_detect(axistext, "Y|y")) {
      out <- out +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())
    }
  } else {
    out <- out +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
  }
  
  if (legend_num == TRUE) {
    out <- out +
      ggplot2::theme(legend.text = ggplot2::element_text(
        family = fontfamily_mono,
        color = light_col,
        size = base_size * .75
      ))
  }
  
  return(out)
}


#' Discrete guide
#'
#' A nicely styled legend for categorical aesthetics (color, fill, shape,
#' size and alpha)
#'
#' @param aes (string) Aesthetic of the legend that should be modified.
#' Options include "color", "fill", "shape", "size" and "alpha".
#' @param reverse (boolean) If TRUE legend items are shown in reversed order.
#' Defaults to FALSE.
#' @param ... Other arguments passed to ggplot methods.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mpg, aes(cty, displ)) + geom_point(aes(color = manufacturer)) +
#'   theme_omics() + guide_discrete()
#' ggplot(mpg, aes(class)) + geom_bar(aes(fill = manufacturer)) +
#'   theme_omics() + guide_discrete(aes = "fill")
#' }
#'
#' @export
guide_discrete <- function(aes = "color", reverse = FALSE, ...) {
  if(!aes %in% c("color", "fill", "shape", "size", "alpha")) stop('aes must be one of "color", "fill", "shape", "size" or "alpha"')
  
  if (aes == "color") {
    out <- ggplot2::guides(
      color = ggplot2::guide_legend(
        reverse = reverse,
        title.position = "top",
        keywidth = grid::unit(1.4, "lines"),
        keyheight = grid::unit(1.4, "lines")
      ))
  }
  if (aes == "fill") {
    out <- ggplot2::guides(
      fill = ggplot2::guide_legend(
        reverse = reverse,
        title.position = "top",
        keywidth = grid::unit(1.4, "lines"),
        keyheight = grid::unit(1.4, "lines")
      ))
  }
  if (aes == "shape") {
    out <- ggplot2::guides(
      shape = ggplot2::guide_legend(
        reverse = reverse,
        title.position = "top",
        keywidth = grid::unit(1.4, "lines"),
        keyheight = grid::unit(1.4, "lines")
      ))
  }
  if (aes == "size") {
    out <- ggplot2::guides(
      size = ggplot2::guide_legend(
        reverse = reverse,
        title.position = "top",
        keywidth = grid::unit(1.4, "lines"),
        keyheight = grid::unit(1.4, "lines")
      ))
  }
  if (aes == "alpha") {
    out <- ggplot2::guides(
      alpha = ggplot2::guide_legend(
        reverse = reverse,
        title.position = "top",
        keywidth = grid::unit(1.4, "lines"),
        keyheight = grid::unit(1.4, "lines")
      ))
  }
  
  return(out)
}


#' A nicely styled colorbar for continuous aesthetics (color, fill and alpha)
#'
#' @param aes (string) Aesthetic of the legend that should be modified.
#' Options include "color", "fill", "shape" and "size".
#' @param aes (type) Type of the guide. Options include "bar" and "steps".
#' @param width (number) Width of the color bar.
#' Options include "color", "fill", "shape" and "size".
#' @param ... Other arguments passed to ggplot methods.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mpg, aes(cty, hwy)) + geom_point(aes(color = displ)) +
#'   theme_omics() + guide_continuous()
#' ggplot(mpg, aes(class)) + geom_bar(aes(fill = manufacturer)) +
#'   theme_omics() + guide_continuous(aes = "fill")
#' }
#'
#' @export
guide_continuous <- function(aes = "color", type = "bar", width = .4, ...) {
  if(!aes %in% c("color", "fill", "alpha")) stop('aes must be one of "color", "fill" or "alpha"')
  if(!type %in% c("bar", "steps")) stop('type must be either "bar" or "steps"')
  
  if (aes == "color") {
    if (type == "bar") {
      out <- ggplot2::guides(
        color = ggplot2::guide_colorbar(
          title.position = "top", 
          label.hjust = 1,
          barwidth = grid::unit(width, "lines")
        ))
    }
    if (type == "steps") {
      out <- ggplot2::guides(
        color = ggplot2::guide_colorsteps(
          title.position = "top", 
          label.hjust = 1,
          barwidth = grid::unit(width, "lines"),
          ticks.colour = "transparent",
          show.limits = TRUE
        ))
    }
  }
  
  
  if (aes == "fill") {
    if (type == "bar") {
        out <- ggplot2::guides(
          fill = ggplot2::guide_colorbar(
            title.position = "top", 
            label.hjust = 1,
            barwidth = grid::unit(width, "lines")
          ))
    }
    if (type == "steps") {
      out <- ggplot2::guides(
        fill = ggplot2::guide_colorbar(
          title.position = "top", 
          label.hjust = 1,
          barwidth = grid::unit(width, "lines"),
          ticks.colour = "transparent",
          show.limits = TRUE
        ))
    }
  }
  
  if (aes == "alpha") {
    if (type == "bar") {
      out <- ggplot2::guides(
        alpha = ggplot2::guide_colorbar(
          title.position = "top", 
          label.hjust = 1,
          barwidth = grid::unit(width, "lines")
        )) 
    }
    if (type == "steps") {
      out <- ggplot2::guides(
        alpha = ggplot2::guide_colorsteps(
          title.position = "top", 
          label.hjust = 1,
          barwidth = grid::unit(width, "lines"),
          ticks.colour = "transparent",
          show.limits = TRUE
        )) 
    }
  }
  
  return(out)
}
