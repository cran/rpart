#SCCS @(#)rpart.poisson.s	1.6 07/05/01
rpart.poisson <- function(y, offset, parms, wt) {
    if (is.matrix(y)) {
	if (ncol(y)!=2) stop("response must be a 2 column matrix or a vector")
	if (!is.null(offset)) y[,1] <- y[,1] * exp(offset)
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
	if(is.null(names(parms))) stop("You must input a named list for parms")
	parmsNames <- c("method", "shrink")
	indx <- pmatch(names(parms), parmsNames, nomatch= 0)
	if (any(indx==0))
               stop("'parms' component not matched: ", names(parms)[indx==0])
	else names(parms) <- parmsNames[indx]

	if (is.null(parms$method)) method <- 1
	else method <- pmatch(parms$method, c("deviance", "sqrt"))
	if (is.null(method)) stop("Invalid error method for Poisson")

	if (is.null(parms$shrink)) shrink <- 2- method
	else shrink <- parms$shrink

	if (!is.numeric(shrink) || shrink <0)
		stop("Invalid shrinkage value")
	parms <- c(shrink=shrink, method=method)
	}

    list(y=y, parms=parms, numresp=2, numy=2,
	 summary= function(yval, dev, wt, ylevel, digits) {
	     paste("  events=", formatg(yval[,2]),
		",  estimated rate=" , formatg(yval[,1], digits),
		" , mean deviance=",formatg(dev/wt, digits),
		sep = "")
	     },
	 text= function(yval, dev, wt, ylevel, digits, n, use.n) {
	     if(use.n) {paste(formatg(yval[,1],digits),"\n",
				formatg(yval[,2]),"/",n,sep="")} else
		    {paste(formatg(yval[,1],digits))}}
	 )

    }
