#' Plot methods for a vlstarpred object
#'
#' Plot method for objects with class attribute \code{vlstarpred}.
#'
#' @param col.ci Character vector, colors for the interval forecast.
#' @param col.obs Character vector, colors for the observed values.
#' @param col.pred Character vector, colors for the predicted values.
#' @param col.vline Character vector, colors for the vertical line.
#' @param lty.ci Vector, lty for the interval forecast.
#' @param lty.obs Vector, lty for the plot of the observed values.
#' @param lty.pred Vector, lty for the plot of the predicted values.
#' @param lty.vline Vector, lty for the vertical line.
#' @param lwd.ci Vector, lwd for the interval forecast.
#' @param lwd.obs Vector, lwd for the plot of the observed values.
#' @param lwd.pred Vector, lwd for the plot of the predicted values.
#' @param lwd.vline Vector, lwd for the vertical line.
#' @param main Character vector, the titles of the plot.
#' @param mar Setting of margins.
#' @param names Character vector, the variables names to be plotted. If left \code{NULL}, all variables are plotted.
#' @param oma Setting of outer margins.
#' @param type Character, if \code{multiple} all plots are drawn in a single device, otherwise the plots are shown consecutively.
#' @param x An object of class \sQuote{\code{vlstarpred}}.
#' @param xlab Character vector signifying the labels for the x-axis.
#' @param ylab Character vector signifying the labels for the y-axis.
#' @param ylim Vector, the limits of the y-axis.
#' @param \dots Passed to internal plot function.
#' @export
#' @author Andrea Bucci
#' @seealso \code{\link{predict.VLSTAR}}

plot.vlstarpred <- function(x, type = c('single', 'multiple'), names = NULL,
                            main = NULL, xlab = NULL, ylab = NULL,
                            lty.obs = 2,lty.pred = 1, lty.ci = 3, lty.vline = 1, lwd.obs = 1, lwd.pred = 1,
                            lwd.ci = 1, lwd.vline = 1, col.obs = NULL, col.pred = NULL, col.ci = NULL,
                            col.vline = NULL, ylim = NULL, mar = par("mar"), oma = par("oma"),...){
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  type <- match.arg(type)
  ifelse(is.null(names), ynames <- colnames(x$y), ynames <- names)
  nrowy <- nrow(x$y)
  if(is.null(col.obs)){
    col.obs = 'black'
  }
  if(is.null(col.ci)){
    col.ci = 'red'
  }
  if(length(col.ci) != 1){
    stop('Please provide a unique color for the confidence interval')
  }
  col.ci = rep(col.ci, 2)
  if(is.null(col.pred)){
    col.pred = 'blue'
  }
  if(is.null(col.vline)){
    col.vline = 'gray'
  }
  col = c(col.pred, col.obs, col.ci, col.vline)
  lty = c(lty.pred, lty.obs, rep(lty.ci,2), lty.vline)
  lwd = c(lwd.pred, lwd.obs, rep(lwd.ci,2), lwd.vline)
  ncoly <- length(ynames)
  nc <- ifelse(ncoly > 4, 2, 1)
  ifelse(is.null(main), main <- paste("Forecast of series", ynames), main <- rep(main, ncoly)[1:ncoly])
  ifelse(is.null(ylab), ylab <- rep("", ncoly), ylab <- rep(ylab, ncoly)[1:ncoly])
  ifelse(is.null(xlab), xlab <- rep("", ncoly), xlab <- rep(xlab, ncoly)[1:ncoly])
  plotprediction <- function(x, main, name, ylab, xlab, col, lty, lwd, ...){
  predy <- c(rep(NA, nrowy - 1), x$y[nrowy, name], x$forecasts[[name]][, 1])
  predl <- c(rep(NA, nrowy - 1), x$y[nrowy, name], x$forecasts[[name]][, 2])
  predu <- c(rep(NA, nrowy - 1), x$y[nrowy, name], x$forecasts[[name]][, 3])
  yobs <- c(x$y[, name], rep(NA, length(x$forecasts[[name]][, 1])))
  if(is.null(ylim)){
    min.y <- min(na.omit(c(predy, predl, predu, yobs)))
    max.y <- max(na.omit(c(predy, predl, predu, yobs)))
    ylim <- c(min.y, max.y)
  }
  plot.ts(predy, main = main, ylim = ylim, col = col[1], lwd=lwd[1], lty = lty[1], xlab = xlab, ylab = ylab, ...)
  lines(yobs, col = col[2], lty = lty[2], lwd = lwd[2])
  lines(predl, col = col[3], lty = lty[3], lwd = lwd[3])
  lines(predu, col = col[4], lty = lty[4], lwd = lwd[4])
  abline(v = nrowy, col = col[5], lty = lty[5], lwd = lwd[5])
  }
  if(type == "single"){
    par(mar = mar, oma = oma)
    if(ncoly > 1) par(ask = TRUE)
    for(i in 1:ncoly){
      plotprediction(x = x, name = ynames[i], main = main[i], col = col, lty = lty, lwd = lwd, ylab = ylab[i], xlab = xlab[i], ...)
    }
  } else if(type == "multiple"){
    nr <- ceiling(ncoly / nc)
    par(mfcol = c(nr, nc), mar = mar, oma = oma)
    for(i in 1:ncoly){
      plotprediction(x = x, name = ynames[i], main = main[i], col = col, lty = lty, lwd = lwd, ylab = ylab[i], xlab = xlab[i], ...)
    }
  }
}

