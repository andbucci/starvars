#'@S3method predict VLSTAR

plot.VLSTAR <- function(object, names = NULL, main.fit = NULL, main.acf = NULL, main.pacf = NULL, main.logi = NULL, ylim.fit = NULL, ylim.resid = NULL, lty.fit = NULL, lty.resid = NULL, lty.logi = NULL, lwd.fit = NULL, lwd.resid = NULL, lwd.logi = NULL, lag.acf = NULL, lag.pacf = NULL, col.fit = NULL, col.resid = NULL, col.logi = NULL,  ylab.fit = NULL, ylab.resid = NULL, ylab.acf = NULL, ylab.pacf = NULL, ylab.logi = NULL, xlab.fit = NULL, xlab.resid = NULL, xlab.logi = NULL, mar = par("mar"), oma = par("oma"), adj.mtext = NA, padj.mtext = NA, col.mtext = NA,...){
  op <- par(no.readonly = TRUE)
  resids <- resid(object)
  fitted <- fitted(object)
  y <- object$Data[[1]]
  ynames <- colnames(y)
  colnames(fitted) <- ynames
  colnames(resids) <- ynames
  if (is.null(names)) {
    names <- ynames
  } else {
    names <- as.character(names)
    if (!(all(names %in% ynames))) {
      warning("\nInvalid variable name(s) supplied, using first variable.\n")
      names <- ynames[1]
    }

  }
  nv <- length(names)
  logistic <- matrix(ncol = ncol(object$yoriginal), nrow = length(object$st))
  for(j in 1:ncol(logistic)){
      for (i in 1:nrow(logistic)){
   logistic[i,j] <- (1+exp(-object$Cgamma[j,1]*(object$st[i]-object$Cgamma[j,2])))^-1
  }
  }
  logistic <- apply(logistic, 2, sort, decreasing = F)
  colnames(logistic) <- ynames

  ifelse(is.null(main.fit), main.fit <- paste("Diagram of fit and residuals for", names), main.fit <- rep(main.fit, nv)[1:nv])
  ifelse(is.null(main.acf), main.acf <- rep("ACF Residuals", nv), main.acf <- rep(main.acf, nv)[1:nv])
  ifelse(is.null(main.pacf), main.pacf <- rep("PACF Residuals", nv), main.pacf <- rep(main.pacf, nv)[1:nv])
  ifelse(is.null(main.logi), main.logi <-  paste("Logistic function for", names), main.logi <- rep(main.logi, nv)[1:nv])
  ifelse(is.null(lty.fit), lty.fit <- c(1, 2), lty.fit <- rep(lty.fit, 2)[1:2])
  ifelse(is.null(lty.resid), lty.resid <- c(1, 1), lty.resid <- rep(lty.resid, 2)[1:2])
  ifelse(is.null(lty.logi), lty.logi <- c(1, 1), lty.logi  <- rep(lty.logi, 2)[1:2])
  ifelse(is.null(lwd.fit), lwd.fit <- c(1, 1), lwd.fit <- rep(lwd.fit, 2)[1:2])
  ifelse(is.null(lwd.resid), lwd.resid <- c(1, 1), lwd.resid <- rep(lwd.resid, 2)[1:2])
  ifelse(is.null(lwd.logi), lwd.logi <- c(1, 1), lwd.logi <- rep(lwd.logi, 2)[1:2])
  ifelse(is.null(lag.acf), lag.acf <- 12, lag.acf <- lag.acf)
  ifelse(is.null(lag.pacf), lag.pacf <- 12, lag.pacf <- lag.pacf)
  ifelse(is.null(col.fit), col.fit <- c("black", "blue"), col.fit <- rep(col.fit, 2)[1:2])
  ifelse(is.null(col.resid), col.resid <- c("black", "red"), col.resid <- rep(col.resid, 2)[1:2])
  ifelse(is.null(col.logi), col.logi <- c("black", "red"), col.logi <- rep(col.logi, 2)[1:2])
  ifelse(is.null(ylab.fit), ylab.fit <- rep("", nv), ylab.fit <- rep(ylab.fit, nv)[1:nv])
  ifelse(is.null(ylab.resid), ylab.resid <- rep("", nv), ylab.resid <- rep(ylab.resid, nv)[1:nv])
  ifelse(is.null(ylab.acf), ylab.acf <- rep("", nv), ylab.acf <- rep(ylab.acf, nv)[1:nv])
  ifelse(is.null(ylab.pacf), ylab.pacf <- rep("", nv), ylab.pacf <- rep(ylab.pacf, nv)[1:nv])
  ifelse(is.null(ylab.logi), ylab.logi <- rep("", nv), ylab.logi <- rep(ylab.logi, nv)[1:nv])
  ifelse(is.null(xlab.fit), xlab.fit <- rep("", nv), xlab.fit <- rep(xlab.fit, nv)[1:nv])
  ifelse(is.null(xlab.resid), xlab.resid <- rep("", nv), xlab.resid <- rep(xlab.resid, nv)[1:nv])
  ifelse(is.null(xlab.logi), xlab.logi <- rep("", nv), xlab.logi <- rep(xlab.logi, nv)[1:nv])

  plotest <- function(y, fitted, resids, logistic, main.fit, main.acf, main.pacf, main.logi, ylab.fit, ylab.resid, ylab.logi,
                      ylab.acf, ylab.pacf, xlab.fit, xlab.resid, xlab.logi, adj.mtext, padj.mtext, col.mtext, ...){
    ifelse(is.null(ylim.fit), ylim.fit <- c(min(c(y, fitted)), max(c(y, fitted))), ylim.fit <- ylim.fit)
    ifelse(is.null(ylim.resid), ylim.resid <- c(min(resids), max(resids)), ylim.resid <- ylim.resid)
    layout(matrix(c(1, 1, 2, 2, 3, 4, 5, 5), nrow = 4, ncol = 2, byrow = TRUE))
    par(oma = c(6, 0, 5, 0), mar = c(0, 5.1, 0, 2.1))
    plot.ts(y, main = "", ylim = ylim.fit, ylab = ylab.fit, xlab = xlab.fit, lty = lty.fit[1], lwd = lwd.fit[1], col = col.fit[1], axes = FALSE)
    lines(fitted, col = col.fit[2], lty = lty.fit[2], lwd = lwd.fit[2])
    box()
    axis(2, pretty(c(y, fitted))[-1])
    mtext(main.fit, side = 3, line = 3, adj = adj.mtext, padj = padj.mtext, col = col.mtext, ...)
    plot.ts(resids, main = "", ylim = ylim.resid, ylab = ylab.resid, xlab = xlab.resid, lty = lty.resid[1], lwd = lwd.resid[1], col = col.resid[1])
    abline(h = 0, col = col.resid[2], lty = lty.resid[2], lwd = lwd.resid[2])
    par(mar=c(1, 5.1, 5, 2.1))
    acf(resids, main = main.acf, ylab = ylab.acf, lag.max = lag.acf, ...)
    pacf(resids, main = main.pacf, ylab = ylab.pacf, lag.max = lag.pacf, ...)
    plot(1:length(logistic),logistic, main = main.logi, ylim = c(0,1), ylab = ylab.logi, xlab = xlab.logi, lty = lty.logi[1], lwd = lwd.logi[1], col = col.logi[1], axes = T)
    }
  par(mar = mar, oma = oma)
  for (i in 1:nv) {
    plotest(y = y[, names[i]], fitted = fitted[, names[i]], resids = resids[, names[i]], logistic = logistic[, names[i]], main.fit = main.fit[i], main.acf = main.acf[i], main.pacf = main.pacf[i], main.logi = main.logi[i], ylab.fit = ylab.fit[i], ylab.resid = ylab.resid[i], ylab.acf = ylab.acf[i], ylab.pacf = ylab.pacf[i], ylab.logi = ylab.logi[i], xlab.fit = xlab.fit[i], xlab.resid = xlab.resid[i], xlab.logi = xlab.logi[i], adj.mtext = adj.mtext, padj.mtext = padj.mtext, col.mtext = col.mtext)
    if (nv > 1) par(ask = TRUE)
  }
  on.exit(par(op))
}
