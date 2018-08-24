#' plot_plasticityIndex
#'
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples data %>% plot_plasticityIndex()

plot_plasticityIndex <- function(data, xvar = strain, dot_color = strain, palette = main, ...) {


  xvar <- quo_name(enquo(xvar))
  colors <- quo_name(enquo(palette))
  p <- ggplot(data, aes_string(x = xvar))

  # if(colorvar) {
  #   dot_color <- quo_name(enquo(colorvar))
  # } else {
  #   dot_color <- quo_name(enquo(xvar))
  # }

  # if(dot_color) {
  #   dot_color <- quo_name(enquo(dot_color))
  # } else {
  #   dot_color <- xvar
  # }

  dot_color <- quo_name(enquo(dot_color))

  if(dot_color == xvar) {
    p <- p + ggbeeswarm::geom_quasirandom(aes_string(y = 'rel.Logit', colour = xvar), width = 0.1)
  } else {
    p <- p + ggbeeswarm::geom_quasirandom(aes_string(y = 'rel.Logit', colour = dot_color), width = 0.1)
  }

    p + theme_classic() +
    add.median('rel.Logit', width = 0.5) +
    add.quartiles('rel.Logit') +
    scale_x_discrete(drop = FALSE) +
    add.n.categorical(!!xvar, ypos = -3) +
    coord_cartesian(ylim = c(-3,3)) +
    labs(y = "plasticity index") +
    figure.axes(no.x = FALSE)

}
