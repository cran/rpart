#SCCS @(#)rpart.poisson.s	1.2 02/16/97
rpart.poisson <- function(y, offset, parms=1) {
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

    list(y=y, parms=parms, numresp=2)
    }
