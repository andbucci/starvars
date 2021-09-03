\name{SSQ}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sum of Squares
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Sum of Squares function to be maximized in the estimation of a VLSTAR model.
}
\usage{
SSQ(param, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
Matrix of initial values for the parameters.
%%     ~~Describe \code{y1} here~~
}
\item{data}{
List of objects useful for the optimization, such as the dependent variables and the estimated covariance matrix.
}
}
%\details{


\value{
  The value of the Sum of Squares.
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
\code{\link{VLSTAR}}, \code{\link{loglike}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~VLSTAR }% use one of  RShowDoc("KEYWORDS")
\concept{VLSTAR}
\concept{Vector Autoregressive Smooth Transition}
\concept{loglike}

