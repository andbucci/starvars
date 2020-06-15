\name{VLSTARjoint}
\alias{VLSTARjoint}
\alias{print.VLSTARjoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Joint linearity test
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function allows the user to test linearity against a Vector Smooth Transition Autoregressive Model with a single transition variable.
}
\usage{
VLSTARjoint(y, exo, st, alpha = 0.05)
\method{print}{VLSTARjoint}(x, \dots)
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
  \item{st}{
  single transition variable for all the equation of dimension \code{(Tx1)}
%%     ~~Describe \code{st} here~~
}
  \item{alpha}{
  Confidence level
%%     ~~Describe \code{constant} here~~
}

\item{x}{\sQuote{VLSTARjoint} object }
\item{\dots}{ further arguments to be passed to and from other methods }
}
\value{
An object of class VLSTARjoint.
}
\details{
Given a VLSTAR model with a unique transition variable, \eqn{s_{1t} = s_{2t} = \dots = s_{\widetilde{n}t} = s_t}, a generalization of the
linearity test presented in Luukkonen, Saikkonen and Terasvirta (1988) may be implemented.

Assuming a 2-state VLSTAR model, such that
\deqn{y_t = B_{1}z_t + G_tB_{2}z_t + \varepsilon_t.}
Where the null \eqn{H_{0} : \gamma_{j} = 0}, \eqn{j = 1, \dots, \widetilde{n}}, is such that \eqn{G_t \equiv (1/2)/\widetilde{n}} and the
previous Equation is linear. When the null cannot be rejected, an identification problem of the parameter \eqn{c_{j}} in the transition function
emerges, that can be solved through a first-order Taylor expansion around \eqn{\gamma_{j} = 0}.

The approximation of the logistic function with a first-order Taylor expansion is given by
\deqn{G(s_t; \gamma_{j},c_{j}) = (1/2) + (1/4)\gamma_{j}(s_t-c_{j}) + r_{jt}}
\deqn{= a_{j}s_t + b_{j} + r_{jt}}
where \eqn{a_{j} = \gamma_{j}/4}, \eqn{b_{j} = 1/2 - a_{j}c_{j}} and \eqn{r_{j}} is the error of the approximation. If \eqn{G_t} is specified as follows
\deqn{G_t = diag\big\{a_{1}s_t + b_{1} + r_{1t}, \dots, a_{\widetilde{n}}s_t+b_{\widetilde{n}} + r_{\widetilde{n}t}\big\}}
\deqn{= As_t + B + R_t}
where \eqn{A = diag(a_{1}, \dots, a_{\widetilde{n}})}, \eqn{B = diag(b_{1},\dots, b_{\widetilde{n}})} e \eqn{R_t = diag(r_{1t}, \dots,
r_{\widetilde{n}t})}, \eqn{y_t} can be written as
\deqn{y_t = B_{1}z_t + (As_t + B + R_t)B_{2}z_t+\varepsilon_t}
\deqn{= (B_{1} + BB_{2})z_t+AB_{2}z_ts_t + R_tB_{2}z_t + \varepsilon_t}
\deqn{= \Theta_{0}z_t + \Theta_{1}z_ts_t+\varepsilon_t^{*}}
where \eqn{\Theta_{0} = B_{1} + B_{2}'B}, \eqn{\Theta_{1} = B_{2}'A} and \eqn{\varepsilon_t^{*} = R_tB_{2} + \varepsilon_t}. Under the null,
\eqn{\Theta_{0} = B_{1}} and \eqn{\Theta_{1} = 0}, while the previous model is linear, with \eqn{\varepsilon_t^{*} = \varepsilon_t}. It
follows that the Lagrange multiplier test, under the null, is derived from the score
\deqn{\frac{\partial \log L(\widetilde{\theta})}{\partial \Theta_{1}} = \sum_{t=1}^{T}z_ts_t(y_t - \widetilde{B}_{1}z_t)'\widetilde{\Omega}^{-1} = S(Y - Z\widetilde{B}_{1})\widetilde{\Omega}^{-1},}
where
\deqn{S = z_{1}'s_{1}\\\vdots\\ z_t's_t}
and where \eqn{\widetilde{B}_{1}} and \eqn{\widetilde{\Omega}} are estimated from the model in \eqn{H_{0}}. If \eqn{P_{Z} = Z(Z'Z)^{-1}Z'} is the
projection matrix of Z, the LM test is specified as follows
\deqn{LM = tr\big\{\widetilde{\Omega}^{-1}(Y - Z\widetilde{B}_{1})'S\big[S'(I_t - P_{Z})S\big]^{-1}S'(Y-Z\widetilde{B}_{1})\big\}.}
Under the null, the test statistics is distributed as a \eqn{\chi^{2}} with \eqn{\widetilde{n}(p\cdot\widetilde{n} + k)} degrees of freedom.}
%%  ~~ If necessary, more details than the description above ~~


\references{
Luukkonen R., Saikkonen P. and Terasvirta T. (1988), Testing Linearity Against Smooth Transition Autoregressive Models. \emph{Biometrika}, 75: 491-499

Terasvirta T. and Yang Y. (2015), Linearity and Misspecification Tests for Vector Smooth Transition Regression Models. \emph{CREATES Research Paper 2014-4}
%% ~put references to the literature/web site here ~
}
\author{
The code was written by Andrea Bucci
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{VLSTAR}} for log-likehood and for NLS estimation of the VLSTAR model.
}
\examples{
data(Realized)
VLSTARjoint(Realized[-1,1:10], Realized[-1,11:16], st=Realized[1:(nrow(Realized)-1),1])
}
