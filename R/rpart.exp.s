#SCCS %W% %G%
# rescaled exponential splitting
rpart.exp <- function(y, offset, parms, wt) {

    if (!inherits(y, "Surv")) 
	   stop("Response must be a survival object - use the Surv() function")

    if (ncol(y) !=2) stop("Can't use (start, stop] data")
    n  <- nrow(y)
    if (any(y[,1]<=0)) stop("Observation time must be >0")
    if (sum(y[,2])==0) stop("No deaths in data set")
    #
    # Rescale time so that the event rate in every interval is = 1
    #
    n  <- nrow(y)
    ord <- order(y[,1])
    temp <- .C("rpartexp", as.integer(n),	
		              as.double(y[ord,]),
		              as.double(wt[ord]),
		              newy = double(n),
	                      double(n))
    newy <- double(n)
    newy[ord] <- temp$newy
    if (length(offset)==n)  newy <- newy * exp(offset)

    if (missing(parms)) parms <- c(shrink=1, method=1)
    else {
	parms <- as.list(parms)
	if (is.null(parms$method)) method <- 1
	else method <- pmatch(parms$method, c("deviance", "sqrt"))
	if (is.na(method)) stop("Invalid error method for Poisson")
	
	if (is.null(parms$shrink)) shrink <- 2-method
	else shrink <- parms$shrink
	if (!is.numeric(shrink) || shrink < 0) 
		stop("Invalid shrinkage value")
	parms <- c(shrink=shrink, method=method)
	}
    list(y=cbind(newy, y[,2]), parms=parms, numresp=2, numy=2,
	 summary= function(yval, dev, wt, ylevel, digits) {
	     paste("  events=", formatg(yval[,2]),
		",  estimated rate=" , formatg(yval[,1], digits),
		" , mean deviance=",formatg(dev/wt, digits),
		sep = "")
	     })
    }
