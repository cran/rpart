#SCCS %W% %G%
rpart.anova <- function(y, offset, parms, wt) {
    if (!is.null(offset)) y <- y-offset
    list(y=y, parms=0, numresp=1, numy=1,
	 summary= function(yval, dev, wt, ylevel, digits ) {
	     paste("  mean=", formatg(yval, digits),
		   ", MSE=" , formatg(dev/wt, digits),
		   sep='')
	     })
    }
