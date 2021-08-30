\name{startingVLSTAR}
\alias{startingVLSTAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Starting parameters for a VLSTAR model
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function allows the user to obtain the set of starting values of Gamma and C for the convergence algorithm via searching grid.
}
\usage{
startingVLSTAR(y, exo = NULL, p = 1,
m = 2, st = NULL, constant = TRUE,
n.combi = NULL,singlecgamma = FALSE)

}
\value{
An object of class \code{startingVLSTAR}.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  \code{data.frame} or \code{matrix} of dependent variables of dimension \code{(Txn)}
%%     ~~Describe \code{y1} here~~
}
  \item{exo}{
  (optional) \code{data.frame} or \code{matrix} of exogenous variables of dimension \code{(Txk)}
%%     ~~Describe \code{x1} here~~
}
  \item{p}{
  lag order
%%     ~~Describe \code{p} here~~
}
  \item{m}{
  number of regimes
%%     ~~Describe \code{m} here~~
}
  \item{st}{
  single transition variable for all the equation of dimension \code{(Tx1)}
%%     ~~Describe \code{st} here~~
}
  \item{constant}{
  \code{TRUE} or \code{FALSE} to include or not the constant
%%     ~~Describe \code{constant} here~~
}

  \item{n.combi}{
  Number of combination for the searching grid of Gamma and C
%%     ~~Describe \code{n.combi} here~~
}
  \item{singlecgamma}{
  \code{TRUE} or \code{FALSE} to use single gamma and c
}
}

\details{
The searching grid algorithm allows for the optimal choice of the parameters \gamma and c by minimizing the sum of the Squared residuals for each possible combination.

The parameter c is initialized by using the mean of the dependent(s) variable, while \gamma is sampled between 0 and 100.
%%  ~~ If necessary, more details than the description above ~~
}
\references{
Anderson H.M. and Vahid F. (1998), Testing multiple equation systems for common nonlinear components. \emph{Journal of Econometrics}. 84: 1-36

Bacon D.W. and Watts D.G. (1971), Estimating the transition between two intersecting straight lines. \emph{Biometrika}. 58: 525-534

Terasvirta T. and Yang Y. (2014), Specification, Estimation and Evaluation of Vector Smooth Transition Autoregressive Models with Applications. \emph{CREATES Research Paper 2014-8}
%% ~put references to the literature/web site here ~
}
\author{
The code was written by Andrea Bucci
%%  ~~who you are~~
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{VLSTARjoint}} to test the presence of a unique transition variable among equations and \code{\link{predict.VLSTAR}} for details on predictions produced for this model; \code{\link{summary}}, \code{\link{coef}}, \code{\link{plot}}, \code{\link{logLik}}
}
\examples{
\donttest{
data(Realized)
y <- Realized[-1,1:10]
y <- y[-nrow(y),]
st <- Realized[-nrow(Realized),1]
st <- st[-length(st)]
starting <- startingVLSTAR(y, p = 1, n.combi = 3,
singlecgamma = FALSE, st = st)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~VLSTAR }% use one of  RShowDoc("KEYWORDS")
