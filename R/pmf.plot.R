#' Plot a discrete cumulative distribution function
#'
#' Plots a discrete/categorical cumulative distribution function.
#' @param x numeric vector giving the possible values of the discrete random variable.
#' @param y numeric vector giving the cumulative probabilities corresponding to `x`.
#' @param sort logical; if `FALSE` (default), the `x` values are assumed to be in
#' ascending order.
#' @param xlim,ylim the x and y limits of the plot.
#' @param xlab,ylab the axis titles.
#' @param main an overall title for the plot.
#' @param add logical; if `TRUE` add to an already existing plot, else generate a new plot.
#' @param verticals logical; if `TRUE`, draw vertical lines at steps.
#' @param do.points logical; if `TRUE` (default), draw points at the knot locations.
#' @param pch point character or symbol if `do.points`.
#' @param col default color of all points and lines.
#' @param col.points color of points if `do.points`.
#' @param cex.points symbol expansion factor if `do.points`.
#' @param col.vert color of vertical lines.
#' @param lty,lwd line type and thickness for horizontal lines.
#' @param lty.vert,lwd.vert line type and thickness for vertical lines.
#' @param ... further graphical arguments (passed to `plot.default()` or to `segments()`).
#' @returns
#' Returns a numeric vector with the random deviates and an attribute `ngen`
#' with the required number of uniform generations.
#' @seealso [plot.stepfun()]
#' @examples
#' n <- 5
#' p <- 0.5
#' x <- 0:n
#' y <- pbinom(x, size = n, prob = p)
#' cdf.plot(x, y,
#'     main = paste0("B(", n,", ", p, ") cumulative distribution function"))
#' @export
# Based on stats::plot.stepfun()
cdf.plot <- function (x, y, sort = FALSE, xlim, ylim = range(c(0, y)),
                  xlab = "x", ylab = "F(x)", main = "", add = FALSE,
                  verticals = TRUE, do.points = TRUE, pch = 20, col = 1,
                  col.points = col, cex.points = 1, col.vert = col,
                  lty = 1, lwd = 2, lty.vert = 3, lwd.vert = 1, ...) {
  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))
  if (sort) {
    i <- order(x)
    x <- x[i]
    y <- y[i]
  }
  if (missing(xlim)) {
    rx <- range(x)
    dr <- if (length(x) > 1L)
              max(0.08 * diff(rx), 0.6*median(diff(x))) else
              abs(x)/16
    xlim <- rx + dr * c(-1, 1)
  } else
    dr <- diff(xlim)
  ti <- c(xlim[1L] - dr, x, xlim[2L] + dr)
  yi <- c(0, y)
  ti.l <- ti[-length(ti)]
  ti.r <- ti[-1L]
  dev.hold()
  on.exit(dev.flush())
  if (add) {
    segments(ti.l, yi, ti.r, yi, col = col, lty = lty, lwd = lwd, ...)
  } else {
    plot(NA, NA, type = "n", xlim = xlim, ylim = ylim, xlab = xlab,
         ylab = ylab, main = main, ...)
    segments(ti.l, yi, ti.r, yi, col = col, lty = lty, lwd = lwd)
  }
  if (do.points)
    points(x, y, pch = pch, col = col.points, cex = cex.points)
  if (verticals)
    segments(x, yi[-length(yi)], x, y, col = col.vert, lty = lty.vert,
             lwd = lwd.vert)
}

