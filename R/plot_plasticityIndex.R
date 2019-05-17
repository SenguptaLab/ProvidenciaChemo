#' plot_plasticityIndex
#'
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples data %>% plot_plasticityIndex()

plot_plasticityIndex <- function(data,
                                 xvar = strain,
                                 dot_color = strain,
                                 palette = main,
                                 BayesFit = FALSE,
                                 width = 0.5,
                                 ...) {

  xvar <- quo_name(enquo(xvar))
  colors <- quo_name(enquo(palette))
  p <- ggplot(data, aes_string(x = xvar))
  width = width

  dot_color <- quo_name(enquo(dot_color))

  if(dot_color == xvar) {
    p <- p + ggbeeswarm::geom_quasirandom(aes_string(y = 'rel.Logit', colour = xvar), width = 0.1)
  } else {
    p <- p + ggbeeswarm::geom_quasirandom(aes_string(y = 'rel.Logit', colour = dot_color), width = 0.1)
  }

    p <- p + theme_classic() +
    add.median('rel.Logit', width = width) +
    add.quartiles('rel.Logit', width = 0.3*width) +
    scale_x_discrete(drop = FALSE) +
    #add.n.categorical(!!xvar, ypos = -3) +
    add.n(!!xvar, y.pos = -3) +
    coord_cartesian(ylim = c(-3,3)) +
    labs(y = "plasticity index") #+
    #figure.axes(no.x = FALSE)

    if(BayesFit) {
      p +
        stat_pointinterval(aes(y=rel.Logit, x = 1.2),
                           data = fitted, fatten_point = 0,
                           size_range = c(0.3, 1), colour = "grey") +
        stat_summary(data = fitted,
                     aes(y=rel.Logit, x = 1.2),
                     fun.y = median,
                     fun.ymin = median,
                     fun.ymax = median,
                     geom = "crossbar",
                     width = 0.05,
                     lwd = 0.35,
                     colour = "grey")
    } else {
      p
    }

}
