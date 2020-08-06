\name{techprices}
\alias{techprices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Daily closing prices of 3 tech stocks.
  %%  ~~function to do ... ~~
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
    This data set contains the series of daily prices of Google, Microsof and Amazon stocks from January 3, 2005 to June 16, 2020, gathered from Yahoo.
}
\usage{
data("techprices")
}
%- maybe also 'usage' for other objects documented here.

\format{
  An \code{xts} object with 3890 daily observations,
  ranging from from January 3, 2005 to June 16, 2020.

    \tabular{rl}{
    \code{Google} \tab daily closing prices of Google (GOOG) stock. \cr
    \code{Microsoft} \tab daily closing prices of Microsoft (MSFT) stock. \cr
    \code{Amazon} \tab daily closing stock prices of Amazon (AMZN) stock. \cr
    }

}

\author{
  Andrea Bucci
  %%  ~~who you are~~
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
