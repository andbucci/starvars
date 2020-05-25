\name{summary}
\alias{summary}
\alias{summary.VLSTAR}
\alias{print.summary.VLSTAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary method for objects of class VLSTAR
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\sQuote{\code{summary}} methods for class \sQuote{\code{VLSTAR}}.
}
\usage{
\method{summary}{VLSTAR}(object, \dots)
\method{print.summary}{VLSTAR}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \sQuote{\code{VLSTAR}} obtained through \command{VLSTAR()}.
%%     ~~Describe \code{y1} here~~
}
\item{x}{
A summary object of class \sQuote{\code{VLSTAR}} obtained through \command{summary()}.
%%     ~~Describe \code{y1} here~~
}
\item{\dots}{further arguments to be passed to and from other methods}
}
%\details{

%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A data frame of predictions for each dependent variable and the (1-\eqn{\alpha}) prediction intervals.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{

Terasvirta T. and Yang Y. (2014), Specification, Estimation and Evaluation of Vector Smooth Transition Autoregressive Models with Applications. \emph{CREATES Research Paper 2014-8}
%% ~put references to the literature/web site here ~
}
\author{
The code was written by Andrea Bucci
%%  ~~who you are~~
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{VLSTAR}}
}
\examples{
##
##See 'VLSTAR' examples
##
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~VLSTAR }% use one of  RShowDoc("KEYWORDS")
\concept{VLSTAR}
\concept{Vector Autoregressive Smooth Transition}
