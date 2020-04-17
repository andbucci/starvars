\name{predict.VLSTAR}
\alias{predict.VLSTAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
VLSTAR Prediction
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function allows to obtain one-step ahead forecasts, with the prediction intervals, of a VLSTAR object.
}
\usage{
predict(object, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  An object of class VLSTAR obtained through VLSTARnls or VLSTARml
%%     ~~Describe \code{y1} here~~
}
  \item{alpha}{
  (optional) Confidence level
%%     ~~Describe \code{x1} here~~
}
}
%\details{

%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A data frame of predictions for each dependent variable and the (1-\alpha) prediction intervals.
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{VLSTARml}} and \code{\link{VLSTARnls}} for log-likehood and nonlinear least squares estimation of the VLSTAR model.
}
\examples{
data(Realized)
fit <- VLSTARnls(Realized[-1,1:10], p = 1, n.combi = 10,
n.iter = 50, st = Realized[-nrow(Realized),1])
predict(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~VLSTAR }% use one of  RShowDoc("KEYWORDS")
