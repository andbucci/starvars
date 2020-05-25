\name{plot}

\alias{plot}

\alias{plot.VLSTAR}


\encoding{latin1}

\title{Plot methods for a VLSTAR object}

\description{
  Plot method for objects with class attribute \code{VLSTAR}.
}

\usage{

\method{plot}{VLSTAR}(x, names = NULL, main.fit = NULL, main.acf = NULL, main.pacf = NULL, main.logi = NULL,
ylim.fit = NULL, ylim.resid = NULL, lty.fit = NULL, lty.resid = NULL, lty.logi = NULL,
lwd.fit = NULL, lwd.resid = NULL, lwd.logi = NULL, lag.acf = NULL, lag.pacf = NULL, col.fit = NULL,
col.resid = NULL, col.logi = NULL,  ylab.fit = NULL, ylab.resid = NULL, ylab.acf = NULL, ylab.pacf = NULL,
ylab.logi = NULL, xlab.fit = NULL, xlab.resid = NULL, xlab.logi = NULL, mar = par("mar"), oma = par("oma"),
adj.mtext = NA, padj.mtext = NA, col.mtext = NA,...)
}

\arguments{
  \item{adj.mtext}{Adjustment for \code{mtext()}.}
  \item{col.fit}{Character vector, colors for diagram of fit.}
  \item{col.mtext}{Character, color for \code{mtext()}, only applicable.}
  \item{col.resid}{Character vector, colors for residual plot.}
  \item{col.logi}{Character vector, colors for logistic function plot.}
  \item{lag.acf}{Integer, lag.max for ACF of residuals.}
  \item{lag.pacf}{Integer, lag.max for PACF of residuals.}
  \item{lty.fit}{Vector, lty for diagram of fit.}
  \item{lty.resid}{Vector, lty for residual plot.}
  \item{lty.logi}{Vector, lty for the plot of the logistic function.}
  \item{lwd.fit}{Vector, lwd for diagram of fit.}
  \item{lwd.resid}{Vector, lwd for residual plot.}
  \item{lwd.logi}{Vector, lwd for the plot of the logistic function.}
  \item{main.acf}{Character vector, main for residuals' ACF.}
  \item{main.fit}{Character vector, main for diagram of fit.}
  \item{main.pacf}{Character vector, main for residuals' PACF.}
  \item{main.logi}{Character vector, main for the plot of the logistic function.}
  \item{mar}{Setting of margins.}
  \item{names}{Character vector, the variables names to be plotted. If
    left \code{NULL}, all variables are plotted.}
  \item{oma}{Setting of outer margins.}
  \item{padj.mtext}{Adjustment for \code{mtext()}.}
  \item{x}{An object of class \sQuote{\code{VLSTAR}} obtained through \command{VLSTAR()}.}
  \item{xlab.fit}{Character vector, xlab for diagram of fit.}
  \item{xlab.resid}{Character vector, xlab for residual plot.}
  \item{xlab.logi}{Character vector, xlab for the plot of the logistic function.}
  \item{ylab.acf}{Character, ylab for ACF.}
  \item{ylab.fit}{Character vector, ylab for diagram of fit.}
  \item{ylab.pacf}{Character, ylab for PACF}
  \item{ylab.resid}{Character vector, ylab for residual plot.}
  \item{ylab.logi}{Character vector, ylab for the plot of the logistic function.}
  \item{ylim}{Vector, the limits of the y-axis.}
  \item{ylim.fit}{Vector, ylim for diagram of fit.}
  \item{ylim.resid}{Vector, ylim for residual plot.}
  \item{\dots}{Passed to internal plot function.}
}

\details{
 The values of the logistic function, given the estimated values of gamma and c through \code{VLSTAR}, are reported.
}

\author{Andrea Bucci}

\references{
  Hamilton, J. (1994), \emph{Time Series Analysis}, Princeton
  University Press, Princeton.

Terasvirta T. and Yang Y. (2014), Specification, Estimation and Evaluation of Vector Smooth Transition Autoregressive Models with Applications. \emph{CREATES Research Paper 2014-8}
}


\seealso{
  \code{\link{VLSTAR}}
}

\examples{
##
##See 'VLSTAR' examples
##
}

\keyword{VLSTAR}
