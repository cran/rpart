#SCCS @(#)rpart.poisson.s	1.3 12/13/99
rpart.poisson <- function(y, offset, parms, wt) {
    if (is.matrix(y)) {
	if (ncol(y)!=2) stop("response must be a 2 column matrix or a vector")
	if (!is.null(offset)) y[,1] <- y[,1] + exp(offset)
	}
    else {
	if (is.null(offset)) y <- cbind(1,y)
	else  y <- cbind( exp(offset), y)
	}
    if (any(y[,1] <=0)) stop("Observation time must be >0")
    if (any(y[,2] <0))  stop("Number of events must be >=0")

    if (missing(parms)) parms <- c(shrink=1, method=1)
    else {
	parms <- as.list(parms)
	if (is.null(parms$method)) method <- 1
	else method <- pmatch(parms$method, c("deviance", "sqrt"))
	if (is.na(method)) stop("Invalid error method for Poisson")

	if (is.null(parms$shrink)) shrink <- 2- method
	else shrink <- parms$shrink

	if (!is.numeric(shrink) || shrink <0)
		stop("Invalid shrinkage value")
	parms <- c(shrink=shrink, method=method)
	}

    list(y=y, parms=parms, numresp=2)
    }
