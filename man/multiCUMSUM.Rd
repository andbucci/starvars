\name{multiCUMSUM}
\alias{multiCUMSUM}
\alias{print.multiCUMSUM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate CUMSUM test
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function returns the test statistics for the presence of co-breaks in a set of multivariate time series.
}
\usage{
multiCUMSUM(data, alpha = 0.95)
\method{print}{multiCUMSUM}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  a (T x N) matrix/data.frame containing the N time series over period \emph{T}
}
\item{alpha}{
Confidence level. By default set to 0.95
}
\item{x}{\sQuote{multiCUMSUM} object}
\item{\dots}{further arguments to be passed to and from other methods}
}
\value{
\item{Lambda Test statistics}{a \emph{test statistics on the presence of a co-break in the conditional mean}}
\item{Omega Test statistics}{a \emph{test statistics on the presence of a co-break in the conditional variance}}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Aue A., Hormann S., Horvath  L.and Reimherr M. (2009), Break detection in the covariance structure of multivariate time series models. \emph{The Annals of Statistics}. 37: 4046-4087

Bai J., Lumsdaine R. L. and Stock J. H. (1998), Testing For and Dating Common Breaks in Multivariate Time Series. \emph{Review of Economic Studies}. 65: 395-432

Barassi M., Horvath L. and Yuqian Z. (2018), Change-Point Detection in the Conditional Correlation Structure of Multivariate Volatility Models. \emph{Journal of Business \& Economic Statistics}


%% ~put references to the literature/web site here ~
}
\author{
The code was written by Andrea Bucci and Giulio Palomba.
%%  ~~who you are~~
}
\examples{
data(Realized)

testCS <- multiCUMSUM(Realized[,1:10], alpha = 0.95)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~VLSTAR }% use one of  RShowDoc("KEYWORDS")
