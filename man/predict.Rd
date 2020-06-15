\name{predict}
\alias{predict.VLSTAR}
\alias{print.vlstarpred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
VLSTAR Prediction
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
One-step or multi-step ahead forecasts, with interval forecast, of a VLSTAR object.
}
\usage{
\method{predict}{VLSTAR}(object, \dots, n.ahead = 1, conf.lev = 0.95,
st.new = NULL, st.num = NULL, newdata = NULL,
method = c('naive', 'Monte Carlo', 'bootstrap'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  An object of class \sQuote{\code{VLSTAR}} obtained through \command{VLSTAR()}.
}
%%     ~~Describe \code{y1} here~~

\item{\dots}{further arguments to be passed to and from other methods}

\item{n.ahead}{An integer specifying the number of ahead predictions
}
\item{conf.lev}{
Confidence level
}
\item{st.new}{
Vector of new data for the transition variable
}
\item{st.num}{
An integer with the index of dependent variable if \code{st.new} is \code{NULL} and the transition variable is a lag of one of the dependent variables
}
\item{method}{
A character identifying which multi-step ahead method should be used among \code{naive}, \code{Monte Carlo} and \code{bootstrap}
}
\item{newdata}{
\code{data.frame} or \code{matrix} of new data for the exogenous variables
}}
%\details{

%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A \code{list} containing:

\item{forecasts}{\code{data.frame} of predictions for each dependent variable and the (1-\eqn{\alpha}) prediction intervals}
\item{y}{a matrix of values for y}
%% ...
}
\references{

Granger C.W.J. and Terasvirta T. (1993), Modelling Non-Linear Economic Relationships. \emph{Oxford University Press}

Lundbergh S. and Terasvirta T. (2007), Forecasting with Smooth Transition Autoregressive Models. \emph{John Wiley and Sons}

Terasvirta T. and Yang Y. (2014), Specification, Estimation and Evaluation of Vector Smooth Transition Autoregressive Models with Applications. \emph{CREATES Research Paper 2014-8}
%% ~put references to the literature/web site here ~
}
\author{
The code was written by Andrea Bucci and Eduardo Rossi
%%  ~~who you are~~
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{VLSTAR}} for log-likehood and nonlinear least squares estimation of the VLSTAR model.
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
