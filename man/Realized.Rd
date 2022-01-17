\name{Realized}
\alias{Realized}
\title{
  Monthly time series used to test VLSTAR models.
}
\description{
    This data set contains the series of realized covariances in 4 stock market indices, i.e. SP-500, Nikkei, DAX, and FTSE, Dividend Yield and Earning Price growth rate, inflation growth rates for U.S., U.K., Japan and Germany, from August 1990 to June 2018.
}
\usage{
  data(Realized)
}

\format{
  A zoo data frame with 334 monthly observations,
  ranging from 1990:M8 until 2018:M6.

    \tabular{rl}{
    \code{SP} \tab Monthly realized variances of S&P 500 index. \cr
    \code{SP-NIKKEI} \tab Monthly realized covariances between S&P 500 and Nikkei. \cr
    \code{SP-FTSE} \tab Monthly realized covariances between S&P 500 and FTSE. \cr
    \code{SP-DAX} \tab Monthly realized covariances between S&P 500 and DAX. \cr
    \code{NIKKEI} \tab Monthly realized variances of Nikkei index. \cr
    \code{NIKKEI-FTSE} \tab Monthly realized covariances between Nikkei and FTSE. \cr
    \code{NIKKEI-DAX} \tab Monthly realized covariances between Nikkei and DAX. \cr
    \code{FTSE} \tab Monthly realized variances of FTSE index. \cr
    \code{FTSE-DAX} \tab Monthly realized covariances between FTSE and DAX. \cr
    \code{DAX} \tab Monthly realized variances of DAX index. \cr
    \code{DP} \tab Monthly Dividends growth rate over the past year relative to current market prices; S&P 500 index. \cr
    \code{EP} \tab Monthly Earnings growth rate over the past year relative to current market prices; S&P500 index. \cr
    \code{Inf_US} \tab US monthly Industrial Production growth. \cr
    \code{Inf_UK} \tab UK monthly Industrial Production growth. \cr
    \code{Inf_JPN} \tab Japan monthly Industrial Production growth. \cr
    \code{Inf_GER} \tab Germany monthly Industrial Production growth. \cr
    }

}

\author{
  Andrea Bucci
}


\seealso{
    \code{\link{rcov}} to build realized covariances from stock prices or returns.
}

\keyword{datasets}
