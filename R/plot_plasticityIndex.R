#' plot_plasticityIndex
#'
#' This function is used to plot a modulation, ie a relative change in
#' log odds of a chemotaxis response. In principle, this function could
#' be used to plot other relative proportion data. The plot is relative to
#' a common control, and expects output from format_AvoidData().
#'
#' @param data dataset used to plot, must be proportion data
#' @param xvar x-axis before facetting. This should be the variable that will
#' divide the nearest data categories. Usually this is "strain" - the food source.
#' @param dot_color factor by which to order the color categories
#' @param palette see eg using plot_color_themes
#' @param width width of horizontal elements, median, error bars.
#' @param BayesFit optional to include a Bayes credile interval based on
#' posterior draws. Must have a "fitted" object which is derived from
#' recenter_fittedValues() function. If using this, the xvar is locked to
#' "data_type" column to align the raw data and fitted values. Additional
#' orienting of the data needs to be done by facetting.
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

  if(BayesFit) {
    xvar = "data_type"
  } else {
    xvar <- quo_name(enquo(xvar))
  }

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
    add.n(!!xvar, y.pos = -3) +
    labs(y = "plasticity index") #+
    #figure.axes(no.x = FALSE)

    if(BayesFit == TRUE & fitted[1,"interval_type"] == "fitted_draws") {
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
      if(BayesFit == TRUE & fitted[1,"interval_type"] == "HDI") {
        p + geom_crossbar(data = fitted,
                          aes(x = 1.2,
                              y = median,
                              ymin = median,
                              ymax = median),
                          width = 0.05,
                          colour = "grey") +
          geom_pointrange(data = fitted,
                        aes(x = 1.2,
                            y = median,
                            ymin = lower.2.5,
                            ymax = upper.97.5),
                            size = 0.5,
                        fatten = 0,
                        color = "grey") +
          geom_pointrange(data = fitted,
                          aes(x = 1.2,
                              y = median,
                              ymin = lower.17,
                              ymax = upper.83),
                          size = 1, fatten = 0,
                          colour = "grey")
      } else {
      p
      }
    }
}
