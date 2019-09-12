multiGW <- function (l1, l2, T, tau, method = c("HAC", "NeweyWest", "Andrews", 
                                      "LumleyHeagerty"), alternative = c("two.sided", "less", 
                                                                         "greater")) 
{
  if (is.matrix(l1) && ncol(l1) > 2) 
    stop("multivariate time series not allowed")
  if (is.matrix(l2) && ncol(l2) > 2) 
    stop("multivariate time series not allowed")
  if (tau < 1) 
    stop("Predictive Horizon must to be a positive integer")
  if (length(l1) != length(l2)) 
    stop("size of l1 and l2 difier")
  alternative <- match.arg(alternative)
  #DNAME <- deparse(substitute(x))
  dif = l1 - l2
  q = length(dif)
  m = T - q
  n = T - tau - m + 1
  delta = mean(dif)
  mod <- lm(dif ~ 0 + rep(1, q))
  if (tau == 1) {
    re = summary(mod)
    STATISTIC = re$coefficients[1, 3]
    if (alternative == "two.sided") 
      PVAL <- 2 * pnorm(-abs(STATISTIC))
    else if (alternative == "less") 
      PVAL <- round(pnorm(STATISTIC), 4)
    else if (alternative == "greater") 
      PVAL <- round(pnorm(STATISTIC, lower.tail = FALSE), 
                    4)
    names(STATISTIC) <- "Normal Standad"
    METHOD <- "Standard Statistic Simple Regression Estimator"
  }
  if (tau > 1) {
    if (method == "HAC") {
      METHOD <- "HAC Covariance matrix Estimation"
      ds = sqrt(vcovHAC(mod)[1, 1])
    }
    if (method == "NeweyWest") {
      METHOD <- "Newey-West HAC Covariance matrix Estimation"
      ds = sqrt(NeweyWest(mod, tau)[1, 1])
    }
    if (method == "LumleyHeagerty") {
      METHOD <- "Lumley HAC Covariance matrix Estimation"
      ds = sqrt(weave(mod)[1, 1])
    }
    if (method == "Andrews") {
      METHOD <- "kernel-based HAC Covariance matrix Estimator"
      ds = sqrt(kernHAC(mod)[1, 1])
    }
    STATISTIC = delta/ds
    if (alternative == "two.sided") 
      PVAL <- 2 * pnorm(-abs(STATISTIC))
    else if (alternative == "less") 
      PVAL <- pnorm(STATISTIC)
    else if (alternative == "greater") 
      PVAL <- pnorm(STATISTIC, lower.tail = FALSE)
    names(STATISTIC) <- "Normal Standar"
  }
  structure(list(statistic = STATISTIC, alternative = alternative, 
                 p.value = PVAL, method = METHOD))
}
