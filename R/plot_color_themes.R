#' plot_color_themes
#'
#' standard plot color themes
#'
#' @section
#'
#' @examples p + scale_fill_plot(palette = "main")
#' p + add.scatter
#' @name plot_color_themes
NULL


#' plotColors
#'makes a list of named colors for the theme
#' @export
#' @rdname plot_color_themes

plotColors <- c(
  `darkGrey` = "#827E7E",
  `red` = "#C93030",
  `green` = "#2F8A34",
  `blue` = "#484CC7",
  `midGrey` = "#B8B1B1",
  `lightGrey` = "#E6DDDD",
  `brown` = "#824B4B",
  #`offWhite` = "gray94",
  `light-blue` = "#A3A3FA"
)

#' plotColorGetter
#' use name of index to pull a vector of hex color codes
#' @export
#' @rdname plot_color_themes

plotColorGetter <- function(...) {
  cols <- c(...)

  if(is.null(cols))
    return(plotColors)

  plotColors[cols]
}

#' PlotColor_palettes
#' @param ... character list set up named themes
#' @export
#' @rdname plot_color_themes

plotColor_palettes <- list(
  `main` = plotColorGetter("darkGrey", "red", "green", "blue", "midGrey", "lightGrey", "brown"),

  `attract` = plotColorGetter("darkGrey", "green", "blue", "midGrey", "lightGrey", "brown"),

  `grey-blue` = plotColorGetter("darkGrey", "blue"),

  `multi-control` = plotColorGetter("darkGrey", "blue", "midGrey", "light-blue"),

  `2-Ps` = plotColorGetter("darkGrey", "blue", "light-blue"),

  `2-each` = plotColorGetter("darkGrey", "lightGrey", "blue", "light-blue"),

  `grey-blue-light` = plotColorGetter("darkGrey", "blue", "light-blue")
)


#' paletteGetter
#' @param palette character for palette selection
#' @param reverse logical to reverse color orders
#' @param ... additional arguments to colorRampPalette (like alpha etc...)
#' @export
#' @rdname plot_color_themes

paletteGetter <- function(palette = "main", reverse = FALSE, ...) {
  pal <- plotColor_palettes[[palette]]

  if (reverse) pal <- rev(pal)

  colorRampPalette(pal, ...)
}

#' scale_color_plot
#' @param palette character for palette selection
#' @param discrete discrete or continuous palette
#' @param reverse logical to reverse color orders
#' @param ... additional arguments to discrete_scale or scale_color_gradientn
#' @export
#' @rdname plot_color_themes

scale_color_plot <- function(palette = "main", discrete = TRUE, reverse = FALSE, drop = drop, ...) {
  pal <- paletteGetter(palette = palette, reverse = reverse)

  if (discrete) {
    discrete_scale("colour", paste0("plotColor_", palette), palette = pal, drop = drop, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}

#' scale_fill_plot
#' @param palette character for palette selection
#' @param discrete discrete or continuous palette
#' @param reverse logical to reverse color orders
#' @param ... additional arguments to discrete_scale or scale_color_gradientn
#' @export
#' @rdname plot_color_themes

scale_fill_plot <- function(palette = "main", discrete = TRUE, reverse = FALSE, drop = drop, ...) {
  pal <- paletteGetter(palette = palette, reverse = reverse)

  if (discrete) {
    discrete_scale("fill", paste0("plotColor_", palette), palette = pal, drop = drop, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}
