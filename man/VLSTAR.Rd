\name{VLSTAR}
\alias{VLSTAR}
\alias{print.VLSTAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
VLSTAR- Estimation
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function allows the user to estimate the coefficients of a VLSTAR model with \emph{m} regimes through maximum likelihood or nonlinear least squares. The set of starting values of Gamma and C for the convergence algorithm can be either passed or obtained via searching grid.
}
\usage{
VLSTAR(y, x = NULL, p = NULL, m = 2, st = NULL, constant = TRUE, n.combi = 50,
       n.iter = 500, starting = NULL, epsilon = 10^(-3), exo = FALSE,
       method = c('ML', 'NLS'))
\method{print}{VLSTAR}(object, \dots)

}
\value{
An object of class VLSTAR, with standard methods.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  Data.frame or matrix of dependent variables of dimension Txn
%%     ~~Describe \code{y1} here~~
}
  \item{x}{
  (optional) Data.frame or matrix of exogenous variables of dimension TxK
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
  single transition variable for all the equation of dimension Tx1
%%     ~~Describe \code{st} here~~
}
  \item{constant}{
  TRUE or FALSE to include or not the constant
%%     ~~Describe \code{constant} here~~
}
  \item{n.combi}{
  Number of combination for the searching grid of Gamma and C
%%     ~~Describe \code{n.combi} here~~
}
  \item{n.iter}{
  number of iteration of the algorithm until forced convergence
%%     ~~Describe \code{n.iter} here~~
}
  \item{starting}{
  (optional) set of intial values for Gamma and C
%%     ~~Describe \code{starting} here~~
}
  \item{epsilon}{
  convergence check measure
%%     ~~Describe \code{epsilon} here~~
}
  \item{exo}{
  TRUE or FALSE to specify if exogenous variables enter in the model
%%     ~~Describe \code{exo} here~~
}

  \item{method}{
  Fitting method: maximum likelihood or nonlinear least squares.
%%     ~~Describe \code{exo} here~~
}
\item{object}{An object of class \sQuote{\code{VLSTAR}} obtained through \command{VLSTAR()}.}
\item{\dots}{further arguments to be passed to and from other methods}
}
\details{
The multivariate smooth transition model is an extension of the smooth transition regression model
introduced by Bacon and Watts (1971)  (see also Anderson and Vahid, 1998). The general model is
\deqn{y_{t} = \mu_0+\sum_{j=1}^{p}\Phi_{0,j}\,y_{t-j}+A_0 x_t \cdot G_t(s_t;\gamma,c)[\mu_{1}+\sum_{j=1}^{p}\Phi_{1,j}\,y_{t-j}+A_1x_t]+\varepsilon_t}
where \eqn{\mu_{0}} and \eqn{\mu_{1}} are the \eqn{\tilde{n} \times 1} vectors of intercepts, \eqn{\Phi_{0,j}} and \eqn{\Phi_{1,j}} are square
\eqn{\tilde{n}\times\tilde{n}} matrices of parameters for lags \eqn{j=1,2,\dots,p}, A_0 and A_1 are \eqn{\tilde{n}\times k} matrices of parameters,
x_t is the \eqn{k \times 1} vector of exogenous variables and \eqn{\varepsilon_t} is the innovation. Finally, \eqn{G_t(s_t;\gamma,c)} is a \eqn{\tilde{n}\times \tilde{n}} diagonal matrix of transition function at time \emph{t}, such that
\deqn{G_t(s_t;\gamma,c)=\{G_{1,t}(s_{1,t};\gamma_{1},c_{1}),G_{2,t}(s_{2,t};\gamma_{2},c_{2}),
    \dots,G_{\tilde{n},t}(s_{\tilde{n},t};\gamma_{\tilde{n}},c_{\tilde{n}})\}.}

Each diagonal element \eqn{G_{i,t}^r} is specified as a logistic cumulative density functions, i.e.
\deqn{G_{i,t}^r(s_{i,t}^r; \gamma_i^r, c_i^r) = \left[1 + \exp\big\{-\gamma_i^r(s_{i,t}^r-c_i^r)\big\}\right]^{-1}}
for \eqn{latex}{i = 1,2, \dots, \tilde{n}} and \eqn{r=0,1,\dots,m-1}, so that the first model is a Vector Logistic Smooth Transition AutoRegressive
(VLSTAR) model.
The ML estimator of \eqn{\theta} is obtained by solving the optimization problem
\deqn{\hat{\theta}_{ML} = arg \max_{\theta}log L(\theta)}
where \eqn{log L(\theta)} is the log-likelihood function of VLSTAR model, given by
\deqn{  ll(y_t|I_t;\theta)=-\frac{T\tilde{n}}{2}\ln(2\pi)-\frac{T}{2}\ln|\Omega|-\frac{1}{2}\sum_{t=1}^{T}(y_t-\tilde{G}_tB\,z_t)'\Omega^{-1}(y_t-\tilde{G}_tB\,z_t)}
The NLS estimators of the VLSTAR model are obtained by solving the optimization problem
\deqn{\hat{\theta}_{NLS} = arg \min_{\theta}\sum_{t=1}^{T}(y_t - \Psi_t'B'x_t)'(y_t - \Psi_t'B'x_t).}
Generally, the optimization algorithm may converge to some local minimum. For this reason, providing valid starting values of \eqn{\theta} is crucial. If there is no clear indication on the initial set of parameters, \eqn{\theta}, this can be done by implementing a grid search. Thus, a discrete grid in the parameter space of \eqn{\Gamma} and C is create to obtain the estimates of B conditionally on each point in the grid. The initial pair of \eqn{\Gamma} and C producing the smallest sum of squared residuals is chosen as initial values, then the model is linear in parameters.
The algorithm is the following:
\enumerate{
\item Construction of the grid for \eqn{\Gamma} and C, computing \eqn{\Psi} for each poin in the grid
\item Estimation of \eqn{\hat{B}} in each equation, calculating the residual sum of squares, \eqn{Q_t}
\item Finding the pair of \eqn{\Gamma} and C providing the smallest \eqn{Q_t}
\item Once obtained the starting-values, estimation of parameters, \emph{B}, via nonlinear least squares (NLS)
\item Estimation of \eqn{\Gamma} and C given the parameters found in step 4
\item Repeat step 4 and 5 until convergence.
}
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
fit.VLSTAR <- VLSTAR(y, p = 1, n.combi = 3,
n.iter = 3, st = st, method ='NLS')

# a few methods for VLSTAR
summary(fit.VLSTAR)
plot(fit.VLSTAR)
predict(fit.VLSTAR, newdata = Realized[nrow(Realized),1:10])}
logLik(fit.VLSTAR)
coef(fit.VLSTAR)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~VLSTAR }% use one of  RShowDoc("KEYWORDS")
