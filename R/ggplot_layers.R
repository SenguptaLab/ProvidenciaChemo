#' ggplot_layers
#'
#' adds layers to ggplots
#'
#' @section
#'
#' @param width optional width for median lines
#' @examples p <- ggplot(aes(x=genotype, y=pct)
#' p + add.scatter
#' @name ggplot_layers
NULL

#' @export
#' @rdname ggplot_layers
#'

add.scatter <- function(yval) {
  yval <- quo_name(enquo(yval))
  ggbeeswarm::geom_quasirandom(aes_string(y = yval),colour = "black", cex=1,
                   width = 0.075,size=0.3, alpha = 0.75,
                   method = 'smiley')
  }
#' @export
#' @rdname ggplot_layers


add.median <- function(yval, width, colour = "black", dodge.width = 0, group = NA) {
  yval <- quo_name(enquo(yval))
  group <- quo_name(enquo(group))
  if(missing(width)) {
    width = 0.25
  } else {
    width = width
  }
  stat_summary(aes_string(y=yval, group = group),fun.y = median,
               fun.ymin = median,
               fun.ymax = median,
               geom = "crossbar",
               width = width,
               lwd = 0.35,
               colour = colour,
               position = position_dodge(width = dodge.width))
}

#' @export
#' @rdname ggplot_layers


add.mean <- function(yval, width, colour = black, dodge.width = 0, group = NA) {
  color <- quo_name(enquo(colour))
  yval <- quo_name(enquo(yval))
  group <- quo_name(enquo(group))
  if(missing(width)) {
    width = 0.25
  } else {
    width = width
  }
  stat_summary(aes_string(y=yval, group = group),fun.y = mean,
               fun.ymin = mean,
               fun.ymax = mean,
               geom = "crossbar",
               width = width,
               lwd = 0.35,
               colour = color,
               position = position_dodge(width = dodge.width))
}

#' @export
#' @rdname ggplot_layers

add.quartiles <- function(yval, width = 0.15, dodge.width = 0, group = NA) {
  yval <- quo_name(enquo(yval))
  group <- quo_name(enquo(group))
  width = width

  stat_summary(aes_string(y = yval, group = group),
               fun.y = median,
               fun.ymin = function(z) {quantile(z,0.25)},
               fun.ymax = function(z) {quantile(z,0.75)},
               geom = "errorbar",
               width = width,
               lwd = 0.4,
               position = position_dodge(width = dodge.width))
}

#' @export
#' @rdname ggplot_layers

figure.axes <- function(no.x = TRUE, ...) {
  if(no.x){
  list(theme(axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 15),
        strip.text.x = ggplot2::element_blank()),
    labs(subtitle = NULL))
  } else {
    list(theme(axis.title.x = ggplot2::element_blank(),
               axis.text.x = ggplot2::element_text(size = 15),
               axis.text.y = ggplot2::element_text(size = 15),
               strip.text.x = ggplot2::element_blank()),
         labs(subtitle = NULL))
    }
}

#' @export
#' @rdname ggplot_layers

add.n.categorical <- function(group, dark, ypos = 0, ...) {
  group <- quo_name(enquo(group))
  ypos <- quo_name(enquo(ypos))
  if(missing(dark)) {
    stat_summary(aes_string(x = as.numeric(as.factor(group)) + 0.3, y = ypos),
                 fun.data = fun_length, geom = "text", size = 3)
  } else {
    stat_summary(aes_string(x = as.numeric(as.factor(group)) + 0.3, y= ypos),
                 fun.data = fun_length, geom = "text", size = 3, color = "white")
  }

}

#' @export
#' @rdname ggplot_layers

add.n <- function(xval, y.pos) {
  if(missing(y.pos)) {
    y.pos <- 0
  } else {
    y.pos <- y.pos
  }
  xval <- quo_name(enquo(xval))
  stat_summary(aes_string(x = xval, y=y.pos),
               fun.data = fun_length, geom = "text", size = 3)
}

#' @export
#' @rdname ggplot_layers

add.Bayes.CI <- function(modsum, emmeans=FALSE, xvar, ...) {
  if(emmeans) {
    xvar <- quo_name(enquo(xvar))

    list(geom_errorbar(data=modsum, aes_(x= xvar,
                                         y=~estimate,
                                         ymin=~lower.HPD,
                                         ymax=~lower.HPD),
                       width=0,colour ="grey", lwd=0.15))
         # geom_errorbar(data=modsum, aes(x=x.pos,y=mean, ymin = lower.25, ymax = upper.75),
         #               width=0,colour = "darkgrey", lwd = 0.15+0.7),
         # geom_segment(data = modsum, aes_(x = xvar,#x.pos-(0.009*nrow(mixed)),
         #                                  y=~estimate,
         #                                xend = xvar,#x.pos+(0.009*nrow(mixed)),
         #                                yend = ~estimate),colour = "darkgrey"))
  } else {
  list(geom_errorbar(data=mixed, aes(x=x.pos,y=mean, ymin=lower.CL, ymax=upper.CL),
                     width=0,colour ="grey", lwd=0.15),
       geom_errorbar(data=mixed, aes(x=x.pos,y=mean, ymin = lower.25, ymax = upper.75),
                     width=0,colour = "darkgrey", lwd = 0.15+0.7),
       geom_segment(data = mixed, aes(x = x.pos-(0.009*nrow(mixed)),
                                      y = mean, xend = x.pos+(0.009*nrow(mixed)),
                                      yend = mean), colour = "darkgrey"))
  }
}

#' @export
#' @rdname ggplot_layers
#alt to theme classic
theme_my_classic <- ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1, size=12),legend.key = ggplot2::element_blank())

#' @export
#' @rdname ggplot_layers
#alt to theme classic
theme_my_ppt <- ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x=ggplot2::element_text(angle=45, hjust=1, size=12),
    legend.key = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    legend.background = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "#E0D8D6", colour = NA))

#' @export
#' @rdname ggplot_layers
#annotate plots with sample size, useful for box and scatter plots
fun_length <- function (x) {
  # annotate plot with sample size
  return(data.frame(y=min(x), label = paste0("(", length(x),")")))
}



#' @export
#' @rdname ggplot_layers
#plotting theme I use for most plots
theme_my <- ggplot2::theme_bw() + ggplot2::theme(
  axis.line        = ggplot2::element_line(colour = "black"),
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  panel.border     = ggplot2::element_blank(),
  strip.background = ggplot2::element_blank(),
  legend.key       = ggplot2::element_blank(),
  axis.text.x=ggplot2::element_text(angle=45, hjust=1, size=12)
)

#' @export
#' @rdname ggplot_layers
#grabs the legend info for separate plotting
g_legend <- function(a.gplot){
  tmp <- ggplotGrob(a.gplot)
  #leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  #legend <- tmp$grobs[[leg]]
  legend <- gtable::gtable_filter(tmp, 'guide-box', fixed = TRUE)
  return(legend)
}

#' @export
#' @rdname ggplot_layers
#plot mean bar, sem errorbar and quasirandom scatter (req. ggbeeswarm)
geom_bardots <- function(fillvar, dotvar, ...) {
  list(stat_summary(geom = "bar", fun.y = mean, aes(fill = {{ fillvar }}), width = 0.5, alpha = 0.75),
       ggbeeswarm::geom_quasirandom(aes(colour = {{ dotvar }}), width = 0.1, alpha = 0.5),
       stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.3))
}

#' @export
#' @rdname ggplot_layers
#plot boxplot and scatter, plus Bayes cred interval from emmeans contrasts, named 'fitted'
geom_relLatency <- function(fitted, fillvar, dotvar, yvar) {
  list(ggbeeswarm::geom_quasirandom(aes(colour = {{ dotvar }},
                                        y = {{ yvar }}),
                                    width = 0.2,
                                    alpha = 0.5,
                                    shape = 16),
       geom_boxplot(aes(fill = {{ fillvar }},
                        y = {{ yvar}}),
                    alpha = 0.7,
                    outlier.shape = NA),
       geom_linerange(data = fitted,
                      aes(x = 1.5,
                          ymin = ll,
                          ymax = hh),
                      lwd = 0.35,
                      colour = "grey34"),
         geom_crossbar(data = fitted,
                       aes(y = m,
                           x = 1.5,
                           ymin = l,
                           ymax = h),
                       width = 0.05,
                       lwd = 0.35,
                       alpha = 0.5,
                       colour = "grey69",
                       fill = "grey69"))
       }


theme_black = function(base_size = 12, base_family = "") {

  theme_classic(base_size = base_size, base_family = base_family) %+replace%

    theme(
      # Specify axis options
      axis.line = element_line(colour = "white"),
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),
      axis.ticks = element_line(color = "white", size  =  0.2),
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),
      axis.ticks.length = unit(0.3, "lines"),
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),
      legend.key = element_rect(color = "white",  fill = "black"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = base_size*0.8, color = "white"),
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),
      legend.position = "right",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      # Specify facetting options
      strip.background = element_blank(), #(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = base_size*0.8, color = "white", face = "italic"),
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),
      plot.title = element_text(size = base_size*1.2, color = "white"),
      plot.margin = unit(rep(1, 4), "lines")

    )

}

#' @export
#' @rdname ggplot_layers
#space out facets in a ggplot object that has multiple x-facets - n_facets is the number of major groups - only works for minor facet grouping
# of n = 1, ie facet_grid(.~genotype + food), where food has 2 levels and genotype has n_facet levels
space_facets <- function(plot, n_facets = 2, major_div = 5, minor_div = 0.3) {
  gt <- ggplot_gtable(ggplot_build(plot))
  gt$widths[seq(1:(n_facets-1))*4 + 4] = major_div*gt$widths[seq(1:(n_facets-1))*4 + 4]
  gt$widths[seq(1:n_facets)*4-2] = minor_div*gt$widths[seq(1:n_facets)*4-2]
  grid::grid.draw(gt)
}
