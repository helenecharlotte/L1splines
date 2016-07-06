##--- packages

library(nleqslv)
library(ggplot2)
library(grid)
library(Rmisc)
library(gridExtra)
library(reshape2)
library(fda)
#library(microbenchmark)
library(xtable)
library(quadprog)
library(matrixcalc)
library(Matrix)
library(rgl)
library(wesanderson)
library(MASS)


##--- function that puts theme on plot

#' ggplot2 theme.
#'
#' User-friendly ggplot2 theme.
#'
#' @param p The ggplot2 object.
#' @param ls Legend size (default = 20).
#' @param ts Title size (default = 15).
#' @param as Axis text size (default = 15).
theme.p = function(p, ls=20, ts=15, as=15) p +
  theme(legend.text     = element_text(size = ls),
        plot.title      = element_text(size = ts),
        axis.text       = element_text(size = as, color="black"),
        axis.title      = element_text(size = as),
        legend.position = "bottom",
        legend.key.size  = unit(0.6, "in"),
        legend.key       = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "gray96"),
        #  panel.grid.minor = element_blank(),
        #  panel.grid.major = element_blank(),
        #   plot.background  = element_rect(fill = "transparent",colour = "white"),
        #  axis.ticks       = element_line(color="black",size = 0.2),
        panel.border     = element_rect(fill=NA, color="black"))

##--- specify path for output

path = "C:/Users/Helene/Documents/Uni/Speciale/latex/graphics/"

##--- function that extracts legend from plot

#' Getlegend.
#'
#' Function that extracts legend from ggplot object.
#'
#' @param p The ggplot2 object.
getlegend.p = function(p){
  tmp = ggplot_gtable(ggplot_build(p))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}

##--- profile function

profileFun = function(FUN, ..., rep=100) {
  Rprof(tmp <- tempfile(), line.profiling=T)
  for (i in 1:rep) FUN(...)
  Rprof(NULL)
  summaryRprof(tmp)
}

##--- bmplot

bmplot = function(p)
  autoplot(p) +
    theme(plot.title       = element_text(face="bold", size=12),
          axis.text        = element_text(size = 12, color="black"),
          axis.title       = element_text(size = 12),
          legend.text      = element_text(size = 12),
          legend.key.size  = unit(0.4, "in"),
          plot.margin      = unit(c(0.1,0.1,0.1,0.1), "in"),
          panel.grid.major = element_line(size = 0.7),
          panel.grid.minor = element_line(size = 0.55),
          strip.text.y     = element_text(size = 16),
          axis.title.y     = element_text(margin=margin(0,15,0,0)),
          axis.title.x     = element_text(margin=margin(15,0,0,0)))

##--- mean squared error

#' Mean squared error.
#'
#' Function that computes the mean squared error of two vectors.
#'
#' @param y Numeric vector.
#' @param yhat Numeric vector.
MSEfun = function(y, yhat)
  return(mean((y - yhat)^2))


##--- color vector

#colvec = wes_palette("Moonrise2", 4, type = "discrete")[c(1,2,4)]
colvec = c("#CC6666", "#66CC99", "#9999CC")

##--- normalizing function

t.norm.fun = function(t)
  return((t - t[1]) / max(t - t[1]))
