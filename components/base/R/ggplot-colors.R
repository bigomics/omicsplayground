#' Function to generate color palettes
#'
#' @param hex Hex code for main color as string
#' @param discrete Reorder for categorical data to place main color first
#' @return A set of three colors
#'
#' @examples 
#' generate_colors("#2EC09C")
#' generate_colors("#BE34EF")
generate_colors <- function(hex, discrete = FALSE) {
  color_set <- c(
    colorspace::desaturate(colorspace::lighten(hex, .6), .2),
    hex,
    colorspace::darken(hex, .6, space = "HLS")
  )
  
  if (discrete == TRUE) color_set <- color_set[c(2, 1, 3)]
  
  return(color_set)
}

#' Function to extract omics colors as hex codes
#'
#' @param ... Character names of colors
#'
#' @examples
#' omics_colors()
#' omics_colors("brand_blue")
#'
#' @export
omics_colors <- function(...) {
  omics_cols <- c(
    ## data colours
    `brand_blue`      = "#004ca7",
    `light_blue`      = "#b8ddff",
    `turquoise`       = "#009c9f",
    `green`           = "#99d4a9",
    `light_green`     = "#c7dd03",
    `orange`          = "#ff9c00",
    `red`             = "#f23451",
    `yellow`          = "#ffd90e",
    `purple`          = "#c7a4ff",
    ## gradient + sequential colours
    `mid_blue`        = "#2780e3", ## plus red for gradient
    `bright_blue`     = "#2fb5e3", ## both plus turquoise, green, light_green and yellow 
    ## grey
    `light_grey`      = "#f8f8f8",
    `mid_grey`        = "#d8d8d8",
    `dark_grey`       = "#989898",
    `super_dark_grey` = "#363535"
  )
  
  cols <- c(...)
  
  if (is.null(cols))
    return (omics_cols)
  
  omics_cols[cols]
}


#' Return function to interpolate a continuous omics color palette
#'
#' @param palette Character name of palette in omics_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments to pass to [grDevices::colorRampPalette()]
#'
#' @examples
#' omics_pal_c()(10)
#' omics_pal_c("red")(3)
#' omics_pal_c("bLUE_rED", reverse = TRUE)(5)
#'
#' @export
omics_pal_c <- function(palette = "brand_blue", reverse = FALSE, ...) {
  
  palette <- stringr::str_to_lower(palette)
  
  if (!palette %in% c("brand_blue", "bright_blue", "red", "turquoise", "yellow", "orange", "purple", "green", "light_green", "grey", "blue_red", "midblue_orange", "turq_yellow", "purple_orange", "blue_red_grey", "turq_yellow_grey", "purple_orange_grey")) stop('palette should be one of "brand_blue", "bright_blue", "red", "turquoise", "yellow", "orange", "purple", "green", "light_green", "grey", "blue_red", "midblue_orange", "turq_yellow", "purple_orange", "blue_red_grey", "turq_yellow_grey", or "purple_orange_grey".')
  
  omics_palettes <- list(
    ## sequential palettes based on gradients of unique data colours
    `brand_blue`         = unname(generate_colors(omics_colors("brand_blue"))),
    `mid_blue`           = unname(generate_colors(omics_colors("mid_blue"))),
    `bright_blue`        = unname(generate_colors(omics_colors("bright_blue"))),
    `turquoise`          = unname(generate_colors(omics_colors("turquoise"))),
    `red`                = unname(generate_colors(omics_colors("red"))),
    `yellow`             = unname(generate_colors(omics_colors("orange"))),
    `orange`             = unname(generate_colors(omics_colors("yellow"))),
    `purple`             = unname(generate_colors(omics_colors("purple"))),
    `green`              = unname(generate_colors(omics_colors("green"))),
    `light_green`        = unname(generate_colors(omics_colors("light_green"))),
    `grey`               = unname(omics_colors(c("mid_grey", "super_dark_grey"))),
    ## diverging palettes based on multiple data colours
    `blue_red`           = unname(omics_colors(c("red", "brand_blue"))),
    `midblue_orange`     = unname(omics_colors(c("orange", "mid_blue"))),
    `turq_yellow`        = unname(omics_colors(c("yellow", "grey", "turquoise"))),
    `purple_orange`      = unname(omics_colors(c("orange", "purple"))),
    `blue_red_grey`      = unname(omics_colors(c("red", "mid_grey", "brand_blue"))),
    `turq_yellow_grey`   = unname(omics_colors(c("yellow", "mid_grey", "turquoise"))),
    `purple_orange_grey` = unname(omics_colors(c("orange", "mid_grey", "purple")))
  )
  
  pal <- omics_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  grDevices::colorRampPalette(pal, ...)
}


#' Return function to interpolate a discrete omics color palette
#'
#' @param palette Character name of palette in omics_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#'
#' @examples
#' omics_pal_d()
#' omics_pal_d("light")
#' omics_pal_d("mUTED", reverse = TRUE)
#'
#' @export
omics_pal_d <- function(palette = "default", reverse = FALSE) {
  
  palette <- stringr::str_to_lower(palette)
  
  if (!palette %in% c("default", "light", "dark", "super_light", "super_dark", "muted", "muted_light", "expanded")) stop('palette should be one of "default", "light", "dark", "super_light", "super_dark", "muted", "muted_light" or "expanded".')
  if (!is.logical(reverse)) stop('reverse should be logical.')
  
  cat_colors <- unname(omics_colors("brand_blue", "red", "turquoise", "yellow", "purple", "green", "orange", "light_green", "bright_blue", "dark_grey"))
  
  omics_palettes <- list(
    `default`     = cat_colors,
    `light`       = colorspace::lighten(cat_colors, .3), 
    `dark`        = colorspace::darken(cat_colors, .2), 
    `super_light` = colorspace::lighten(cat_colors, .6), 
    `super_dark`  = colorspace::darken(cat_colors, .4),
    `muted`       = colorspace::desaturate(cat_colors, .5),
    `muted_light` = colorspace::desaturate(colorspace::lighten(cat_colors, .6), .3), 
    `expanded`    = c(cat_colors, colorspace::darken(cat_colors, .4), colorspace::desaturate(cat_colors, .5), colorspace::lighten(cat_colors, .6))
  )
  
  pal <- omics_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)  
  
  scales::manual_pal(pal)
}


#' Color scale constructor for continuous omics color palettes
#'
#' @param palette Character name of palette in omics_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE
#'            or FALSE
#' @examples
#' library(ggplot2)
#' ggplot(mpg, aes(hwy, cty, color = displ)) + 
#'   geom_point(size = 4) + 
#'   scale_color_omics_c()
#' ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Petal.Width)) + 
#'   geom_point(size = 4) + 
#'   scale_color_omics_c("TURQ_yellow_GREY", reverse = TRUE)
#' @export
scale_color_omics_c <- function(palette = "brand_blue", reverse = FALSE, ...) {
  
  palette <- stringr::str_to_lower(palette)
  
  if (!palette %in% c("brand_blue", "bright_blue", "red", "turquoise", "yellow", "orange", "purple", "green", "light_green", "grey", "blue_red", "midblue_orange", "turq_yellow", "purple_orange", "blue_red_grey", "turq_yellow_grey", "purple_orange_grey")) stop('palette should be one of "brand_blue", "bright_blue", "red", "turquoise", "yellow", "orange", "purple", "green", "light_green", "grey", "blue_red", "turq_yellow", "purple_orange", "blue_red_grey", "turq_yellow_grey", or "purple_orange_grey".')
  
  pal <- omics_pal_c(palette = palette, reverse = reverse)
  
  ggplot2::scale_color_gradientn(colours = pal(256), ...)
}

#' Color scale constructor for discrete omics colors
#'
#' @param palette Character name of palette in omics_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE
#'            or FALSE
#'
#' @examples
#' library(ggplot2)
#' ggplot(mpg, aes(displ, cty, color = class)) +
#'   geom_point(size = 4) + 
#'   scale_color_omics_d()
#' ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Species)) +
#'   geom_point(size = 4) + 
#'   scale_color_omics_d("mUTED", reverse = TRUE)
#'   
#' @export
scale_color_omics_d <- function(palette = "default", reverse = FALSE, ...) {
  
  palette <- stringr::str_to_lower(palette)
  
  if (!palette %in% c("default", "light", "dark", "super_light", "super_dark", "muted", "muted_light", "expanded")) stop('palette should be one of "default", "light", "super_light", "super_dark", "muted", "muted_light" or "expanded".')
  if (!is.logical(reverse)) stop('reverse should be logical.')
  
  pal <- omics_pal_d(palette = palette, reverse = reverse)
  
  ggplot2::discrete_scale("colour", paste0("omics_", palette), palette = pal, ...)
}

#' Fill scale constructor for continuous omics color palettes
#'
#' @param palette Character name of palette in omics_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE
#'            or FALSE
#'
#' @examples
#' library(ggplot2)
#' ggplot(mpg, aes(displ, cty, fill = hwy)) +
#'   geom_point(shape = 21, size = 4, stroke = 1) + 
#'   scale_fill_omics_c()
#' ggplot(iris, aes(Sepal.Width, Sepal.Length, fill = Petal.Width)) +
#'   geom_point(shape = 21, size = 4, stroke = 1) +
#'   scale_fill_omics_c("BLUE_RED_GREY", reverse = TRUE)
#' 
#' @export
scale_fill_omics_c <- function(palette = "brand_blue", reverse = FALSE, ...) {
  
  palette <- stringr::str_to_lower(palette)
  
  if (!palette %in% c("brand_blue", "bright_blue", "red", "turquoise", "yellow", "orange", "purple", "green", "light_green", "grey", "blue_red", "midblue_orange", "turq_yellow", "purple_orange", "blue_red_grey", "turq_yellow_grey", "purple_orange_grey")) stop('palette should be one of "brand_blue", "bright_blue", "red", "turquoise", "yellow", "orange", "purple", "green", "light_green", "grey", "blue_red", "midblue_orange", "turq_yellow", "purple_orange", "blue_red_grey", "turq_yellow_grey", or "purple_orange_grey".')
  
  pal <- omics_pal_c(palette = palette, reverse = reverse)
  
  ggplot2::scale_fill_gradientn(colours = pal(256), ...)
}



#' Fill scale constructor for discrete omics color palettes
#'
#' @param palette Character name of palette in omics_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE
#'            or FALSE
#'
#' @examples
#' library(ggplot2)
#' ggplot(mpg, aes(class, fill = class)) + 
#'   geom_bar() + 
#'   scale_fill_omics_d()
#' ggplot(mpg, aes(class, fill = class)) + 
#'   geom_bar() +
#'   scale_fill_omics_d("mUTED", reverse = TRUE)
#' ggplot(iris, aes(Species, Sepal.Width, fill = Species)) + 
#'   geom_jitter(shape = 21, size = 4, stroke = 1) +
#'   scale_fill_omics_d("super_light")
#'
#' @export
scale_fill_omics_d <- function(palette = "default", reverse = FALSE, ...) {
  
  palette <- stringr::str_to_lower(palette)
  
  if (!palette %in% c("default", "light", "dark", "super_light", "super_dark", "muted", "muted_light", "expanded")) stop('palette should be one of "default", "light", "dark", "super_light", "super_dark", "muted", "muted_light" or "expanded".')
  if (!is.logical(reverse)) stop('reverse should be logical.')
  
  pal <- omics_pal_d(palette = palette, reverse = reverse)
  
  ggplot2::discrete_scale("fill", paste0("omics_", palette), palette = pal, ...)
}
