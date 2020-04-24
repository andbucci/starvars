\name{rcov}
\alias{rcov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Realized Covariance
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function returns the vectorization of the lowest triangular of the Realized Covariance matrices for different frequency.
}
\usage{
rcov(data, freq, ret = FALSE,
    cholesky = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  a (T x N) xts object containing the N price/return series over period \emph{T}
%%     ~~Describe \code{y1} here~~
}
  \item{freq}{
  a string defining the desired frequency for the Realized Covariance matrices between "daily", "monthly", "quarterly" or "yearly"
%%     ~~Describe \code{x1} here~~
}
  \item{ret}{
  boolean, in case it is TRUE the data are converted to returns, FALSE by default
%%     ~~Describe \code{p} here~~
}
  \item{cholesky}{
  boolean, in case it is TRUE the Cholesky factors of the Realized Covariance matrices are calculated, FALSE by default
%%     ~~Describe \code{exo} here~~
}
}
\value{
\item{Realized Covariances}{a \eqn{M \times N(N+1)/2} matrix of realized covariances, where \emph{M} is the number of lower frequency data}
\item{Cholesky Factors (optional)}{a \eqn{M \times N(N+1)/2} matrix of Cholesky factors of the realized covariance matrices, where \emph{M} is the number of lower frequency data}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Andersen T.G., Bollerslev T., Diebold F.X. and Labys P. (2003), Modeling and Forecasting Realized Volatility. \emph{Econometrica}. 71: 579-625

Barndorff-Nielsen O.E. and Shephard  N. (2002), Econometric analysis of realised volatility and its use in estimating stochastic volatility models. \emph{Journal of the Royal Statistical Society}. 64(2): 253-280
%% ~put references to the literature/web site here ~
}
\author{
The code was written by Andrea Bucci
%%  ~~who you are~~
}

\examples{
require(highfrequency)

data(sample_5minprices_jumps)

rc <- rcov(sample_5minprices_jumps, freq = 'daily', cholesky = TRUE, ret = TRUE)
rc
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~VLSTAR }% use one of  RShowDoc("KEYWORDS")
