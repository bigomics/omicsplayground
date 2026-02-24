##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


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
    ## main colors
    `brand_blue`      = "#3181de", ## BS color (primary)
    `dark_orange`     = "#BC7408",
    `turquoise`       = "#009c9f", ## BS color (success)
    `terra_cotta`     = "#bf616a", ## BS color (danger)
    `purple`          = "#A06EC5",
    `green`           = "#5B9B5B",
    `pastel_blue`     = "#81a1c1", ## BS color (secondary)
    `dark_grey`       = "#989898",
    ## back-up colors
    `red`             = "#f23451",
    `orange`          = "#e3a45a", ## BS color (warning)
    `bright_blue`     = "#2fb5e3",
    ## grey shades
    `super_dark_grey` = "#3b4252", ## BS color (dark)
    `mid_grey`        = "#d8d8d8",
    `grey`            = "#eeeeee"
  )

  cols <- c(...)

  if (is.null(cols)) {
    return(omics_cols)
  }

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

  if (!palette %in% c("brand_blue", "dark_orange", "turquoise", "terra_cotta", "purple", "green", "pastel_blue", "grey", "red", "bright_blue", "orange", "blue_red", "blue_orange", "pastel_blue_orange", "turq_orange", "purple_orange", "blue_red_grey", "blue_orange_grey", "pastel_blue_orange_grey", "turq_orange_grey", "purple_orange_grey")) stop('palette should be one of "brand_blue", "dark_orange", "turquoise", "terra_cotta", "purple", "green", "pastel_blue", "grey", "red", "bright_blue", "orange", "blue_red", "blue_orange", "pastel_blue_orange", "turq_orange", "purple_orange", "blue_red_grey", "blue_orange_grey", "pastel_blue_orange_grey", "turq_orange_grey" or "purple_orange_grey".')
  if (!is.logical(reverse)) stop("reverse should be logical.")

  omics_palettes <- list(
    ## sequential palettes based on gradients of unique data colours
    `brand_blue` = unname(generate_colors(omics_colors("brand_blue"))),
    `dark_orange` = unname(generate_colors(omics_colors("dark_orange"))),
    `turquoise` = unname(generate_colors(omics_colors("turquoise"))),
    `terra_cotta` = unname(generate_colors(omics_colors("terra_cotta"))),
    `purple` = unname(generate_colors(omics_colors("purple"))),
    `green` = unname(generate_colors(omics_colors("green"))),
    `pastel_blue` = unname(generate_colors(omics_colors("pastel_blue"))),
    `red` = unname(generate_colors(omics_colors("red"))),
    `bright_blue` = unname(generate_colors(omics_colors("bright_blue"))),
    `orange` = unname(generate_colors(omics_colors("orange"))),
    `grey` = unname(omics_colors(c("grey", "super_dark_grey"))),
    ## diverging palettes based on multiple data colours
    `blue_red` = unname(omics_colors(c("red", "brand_blue"))),
    `blue_orange` = unname(omics_colors(c("dark_orange", "bright_blue"))),
    `pastel_blue_orange` = unname(omics_colors(c("orange", "pastel_blue"))),
    `turq_orange` = unname(omics_colors(c("dark_orange", "turquoise"))),
    `blue_red_grey` = unname(omics_colors(c("red", "grey", "brand_blue"))),
    `blue_orange_grey` = unname(omics_colors(c("dark_orange", "grey", "brand_blue"))),
    `pastel_blue_orange_grey` = unname(omics_colors(c("orange", "grey", "pastel_blue"))),
    `turq_orange_grey` = unname(omics_colors(c("dark_orange", "grey", "turquoise")))
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

  ## Custom gradient palette: interpolate 3 user-defined colours

  if (palette == "custom_gradient") {
    ct <- get_color_theme()
    cols <- c(ct$palette_c1, ct$palette_c2, ct$palette_c3)
    if (reverse) cols <- rev(cols)
    return(function(n) {
      grDevices::colorRampPalette(cols)(n)
    })
  }

  if (!palette %in% c("default", "light", "dark", "super_light", "super_dark", "muted", "muted_light", "expanded", "highlight_blue", "highlight_red", "highlight_orange")) stop('palette should be one of "default", "light", "dark", "super_light", "super_dark", "muted", "muted_light", "expanded", "highlight_blue", "highlight_red" or "highlight_orange".')
  if (!is.logical(reverse)) stop("reverse should be logical.")

  cat_colors <- unname(omics_colors("brand_blue", "dark_orange", "turquoise", "terra_cotta", "purple", "green", "pastel_blue", "dark_grey"))
  hl1_colors <- c(unname(omics_colors("brand_blue")), colorRampPalette(omics_colors("grey", "dark_grey"))(5))
  hl2_colors <- c(unname(omics_colors("red")), colorRampPalette(omics_colors("grey", "dark_grey"))(5))
  hl3_colors <- c(unname(omics_colors("dark_orange")), colorRampPalette(omics_colors("grey", "dark_grey"))(5))

  omics_palettes <- list(
    `default`          = cat_colors,
    `light`            = colorspace::lighten(cat_colors, .3),
    `dark`             = colorspace::darken(cat_colors, .2),
    `super_light`      = colorspace::lighten(cat_colors, .6),
    `super_dark`       = colorspace::darken(cat_colors, .4),
    `muted`            = colorspace::desaturate(cat_colors, .5),
    `muted_light`      = colorspace::desaturate(colorspace::lighten(cat_colors, .6), .3),
    `expanded`         = c(cat_colors, colorspace::darken(cat_colors, .4), colorspace::lighten(cat_colors, .6), colorspace::desaturate(cat_colors, .5)),
    `highlight_blue`   = hl1_colors,
    `highlight_red`    = hl2_colors,
    `highlight_orange` = hl3_colors
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
#'   scale_color_omics_c("turq_orange_GREY", reverse = TRUE)
#' @export
scale_color_omics_c <- function(palette = "brand_blue", reverse = FALSE, ...) {
  palette <- stringr::str_to_lower(palette)

  if (!palette %in% c("brand_blue", "dark_orange", "turquoise", "terra_cotta", "purple", "green", "pastel_blue", "grey", "red", "bright_blue", "orange", "blue_red", "blue_orange", "pastel_blue_orange", "turq_orange", "purple_orange", "blue_red_grey", "blue_orange_grey", "pastel_blue_orange_grey", "turq_orange_grey", "purple_orange_grey")) stop('palette should be one of "brand_blue", "dark_orange", "turquoise", "terra_cotta", "purple", "green", "pastel_blue", "grey", "red", "bright_blue", "orange", "blue_red", "blue_orange", "pastel_blue_orange", "turq_orange", "purple_orange", "blue_red_grey", "blue_orange_grey", "pastel_blue_orange_grey", "turq_orange_grey" or "purple_orange_grey".')
  if (!is.logical(reverse)) stop("reverse should be logical.")

  pal <- omics_pal_c(palette = palette, reverse = reverse)

  ggplot2::scale_color_gradientn(colours = pal(256), ...)
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
#'   scale_fill_omics_c("blue_red_grey", reverse = TRUE)
#'
#' @export
scale_fill_omics_c <- function(palette = "brand_blue", reverse = FALSE, ...) {
  palette <- stringr::str_to_lower(palette)

  if (!palette %in% c("brand_blue", "dark_orange", "turquoise", "terra_cotta", "purple", "green", "pastel_blue", "grey", "red", "bright_blue", "orange", "blue_red", "blue_orange", "pastel_blue_orange", "turq_orange", "purple_orange", "blue_red_grey", "blue_orange_grey", "pastel_blue_orange_grey", "turq_orange_grey", "purple_orange_grey")) stop('palette should be one of "brand_blue", "dark_orange", "turquoise", "terra_cotta", "purple", "green", "pastel_blue", "grey", "red", "bright_blue", "orange", "blue_red", "blue_orange", "pastel_blue_orange", "turq_orange", "purple_orange", "blue_red_grey", "blue_orange_grey", "pastel_blue_orange_grey", "turq_orange_grey" or "purple_orange_grey".')
  if (!is.logical(reverse)) stop("reverse should be logical.")

  pal <- omics_pal_c(palette = palette, reverse = reverse)

  ggplot2::scale_fill_gradientn(colours = pal(256), ...)
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

  if (!palette %in% c("default", "light", "dark", "super_light", "super_dark", "muted", "muted_light", "expanded", "highlight_blue", "highlight_red", "highlight_orange")) stop('palette should be one of "default", "light", "dark", "super_light", "super_dark", "muted", "muted_light", "expanded", "highlight_blue", "highlight_red" or "highlight_orange".')
  if (!is.logical(reverse)) stop("reverse should be logical.")

  pal <- omics_pal_d(palette = palette, reverse = reverse)

  ggplot2::discrete_scale("colour", paste0("omics_", palette), palette = pal, ...)
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

  if (!palette %in% c("default", "light", "dark", "super_light", "super_dark", "muted", "muted_light", "expanded", "highlight_blue", "highlight_red", "highlight_orange")) stop('palette should be one of "default", "light", "dark", "super_light", "super_dark", "muted", "muted_light", "expanded", "highlight_blue", "highlight_red" or "highlight_orange".')
  if (!is.logical(reverse)) stop("reverse should be logical.")

  pal <- omics_pal_d(palette = palette, reverse = reverse)

  ggplot2::discrete_scale("fill", paste0("omics_", palette), palette = pal, ...)
}

#' @export
color_from_middle <- function(data, color1, color2) {
  ## from https://stackoverflow.com/questions/33521828/
  max_val <- max(abs(data), na.rm = TRUE)
  DT::JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'", max_val, color1, max_val, color1, color2, color2, max_val, max_val))
}
