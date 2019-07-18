\name{lrvarbart}
\alias{lrvarbart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Long-run variance using Bartlett kernel
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function returns the long-run variance of a time series, relying on the Bartlett kernel. The window size of the kernel is the cube root of the sample size.
}
\usage{
lrvarbrt(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a (T x 1) vector containing the time series over period \emph{T}
%%     ~~Describe \code{y1} here~~
}
}
\value{
\item{lrv}{long-run variance}
\item{bandwidth}{size of the window}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Hamilton J. D. (1994), Time Series Analysis. \emph{Princeton University Press}

Tsay R.S. (2005), Analysis of Financial Time Series. \emph{John Wiley & SONS}


%% ~put references to the literature/web site here ~
}
\author{
The code was written by Andrea Bucci.
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(Realized)

lrvarbart(Realized[,1])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~VLSTAR }% use one of  RShowDoc("KEYWORDS")
