\name{logLik}

\alias{logLik}
\alias{logLik.VLSTAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-Likelihood method
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns the log-Likelihood of a VLSTAR object.
}
\usage{
\method{logLik}{VLSTAR}(object, type = c('Univariate', 'Multivariate'), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \sQuote{\code{VLSTAR}} obtained through \command{VLSTAR()}.
%%     ~~Describe \code{y1} here~~
}
\item{type}{
Type of Log-Likelihood to be showed (univariate or multivariate).
}
\item{\dots}{further arguments to be passed to and from other methods}
}
%\details{

%%  ~~ If necessary, more details than the description above ~~
%}
\details{
  The log-likelihood of a VLSTAR model is defined as:
  \deqn{
\log l(y_t|I_t;\theta)=-\frac{T\tilde{n}}{2}\ln(2\pi)-\frac{T}{2}\ln|\Omega|-\frac{1}{2}\sum_{t=1}^{T}(y_t-\tilde{G}_tB\,z_t)'\Omega^{-1}(y_t-\tilde{G}_tB\,z_t)
  }
}

\value{
  An object with class attribute \code{logLik}.
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
\concept{logLik}

