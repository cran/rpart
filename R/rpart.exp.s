#SCCS @(#)rpart.exp.s	1.3 01/21/97
# rescaled exponential splitting
rpart.exp <- function(y, offset, parms=1) {
    # if both a time and an offset occur, apply the offset AFTER rescaling
    #
    late.off <-1

    if (!inherits(y, "Surv")) stop("Response must be a survival object - use the Surv() function")


    if (is.matrix(y)) {
	if (ncol(y)!=2) stop("response must be a 2 column matrix or a vector")
	if (!is.null(offset)) late.off <- exp(offset)
	}
    else {
	if (is.null(offset)) stop("No time value given")
	else  y <- cbind( exp(offset), y)
	}

    if (any(y[,1]<=0)) stop("Observation time must be >0")
    if (any(y[,2]!=0 & y[,2]!=1))  stop("Number of events must be 0 or 1")

    stat <- y[,2]
    dtime <- sort(unique(y[stat==1,1]))
    # Next line avoids round off errors that effect pcount.  Pyears2 counts
    #  forward in time, and dtime exactly= death times is not good.
    dtime2 <- dtime + max(dtime)*sqrt(.Machine$double.eps)
    
    n <- length(stat)
    nd<- length(dtime)
    fit <- .C("pyears2", sn=as.integer(n),
			 sny=as.integer(2),
			 sdoevent=as.integer(1),
	                 sy=y,
			 sodim=as.integer(1),
			 ofac=as.integer(0),
			 odim=as.integer(nd),
			 socut=c(0,dtime2),
			 odata=double(n),  #a vector of zeros
			 pyears= double(nd),
			 pn = double(nd),
			 pcount = double(nd),
			 offtable = 0)

    tmp <- y[,1] #no hazard accumulates if no more deaths
    tmp[tmp>max(dtime)] <- max(dtime,na.rm=T)

    y[,1] <- approx(c(0,dtime2), cumsum(c(0, 
                   (fit$pcount/fit$pyears)*diff(c(0,dtime2)))),tmp)$y*late.off 
    list(y=y, parms=parms, numresp=2)
    }
