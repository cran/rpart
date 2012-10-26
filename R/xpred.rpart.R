##
##  Get a set of cross-validated predictions
##  Most of the setup is identical to the rpart routine
##
xpred.rpart <- function(fit, xval=10, cp, return.all=FALSE)
{
    if (!inherits(fit, 'rpart')) stop("Invalid fit object")

    method <- fit$method
    method.int <- pmatch(method, c("anova", "poisson", "class", "user", "exp"))
    if (method.int==5L) method.int <- 2L
    Terms <- fit$terms

    Y <- fit$y
    X <- fit$x
    wt<- fit$wt
    numresp <- fit$numresp

    if (is.null(Y) || is.null(X)) {
	m <- fit$model
	if (is.null(m)) {
	    m <-fit$call[match(c("", 'formula', 'data', 'weights', 'subset',
                                 'na.action'),
                               names(fit$call), nomatch=0)]
	    if (is.null(m$na.action)) m$na.action<- na.rpart
	    m[[1]] <- as.name("model.frame")
	    m <- eval(m, parent.frame())
        }
	if (is.null(X)) X <- rpart.matrix(m)
	if (is.null(wt)) wt <- model.extract(m, "weights")
	if (is.null(Y)) {
	    yflag <- TRUE
	    Y <- model.extract(m, "response")
            offset <- attr(Terms, "offset")
	    if (method != "user") {
		init <- (get(paste("rpart", method, sep='.')))(Y,offset, NULL)
		Y <- init$y
		if (is.matrix(Y)) numy <- ncol(Y) else numy <- 1
            }
        }
	else {
	    yflag <- FALSE
            if (is.matrix(Y)) numy <- ncol(Y) else numy <- 1
        }
    }
    else {
	yflag <- FALSE
        if (is.matrix(Y)) numy <- ncol(Y) else numy <- 1
	offset <- 0
    }

    nobs <- nrow(X)
    nvar <- ncol(X)
    if (length(wt)==0) wt <- rep(1.0, nobs)

    cats <- rep(0, nvar)
    xlevels <- attr(fit, "xlevels")
    if (!is.null(xlevels)){
        cats[match(names(xlevels), dimnames(X)[[2]])] <-
            unlist(lapply(xlevels, length))
    }

    controls <- fit$control
    if (missing(cp)) {
	cp<- fit$cptable[,1L]
	cp <- sqrt(cp * c(10, cp[-length(cp)]))
	cp[1] <- (1+fit$cptable[1,1])/2
    }
    ncp <- length(cp)

    if (length(xval)==1L) {
	## make random groups
	xgroups <- sample(rep(1:xval, length=nobs), nobs, replace=FALSE)
    }
    else if (length(xval) == nrow(X)) {
	xgroups <- xval
	xval <- length(unique(xgroups))
    }
    else {
	## Check to see if observations were removed due to missing
	if (!is.null(fit$na.action)) {
	    ## if na.rpart was used, then na.action will be a vector
	    temp <- as.integer(fit$na.action)
	    xval <- xval[-temp]
	    if (length(xval) == nobs) {
		xgroups <- xval
		xval <- length(unique(xgroups))
            }
	    else stop("Wrong length for 'xval'")
        }
	else stop("Wrong length for 'xval'")
    }

    costs <- fit$call$costs
    if (is.null(costs)) costs <- rep(1.0, nvar)

    parms <- fit$parms
    if (method=="user") {
	mlist <- fit$functions

	## If yflag==TRUE, then y was retrieved from the original data and we
        ##   need to call init to check and possibly transform it.
        ## If false, then we have one of the few differences with the
        ##  rpart setup
        if (yflag) {
            if (length(parms)==0L) init <- mlist$init(Y, offset,,wt)
            else                   init <- mlist$init(Y, offset, parms, wt)
            Y <- init$Y
            numy <- init$numy
            parms <- init$parms
        }
        else {
            if (is.matrix(Y)) numy <- ncol(Y) else numy <- 1
            init <- list(numresp=numresp, numy=numy, parms=parms)
        }
        keep <- rpartcallback(mlist, nobs, init) #assign to a variable to
                                        #stop garbage collection
	method.int <- 4L            #the fourth entry in func_table.h
    }

    ## Finally do the work
    if (is.matrix(Y))  Y<- as.double(t(Y))
    else storage.mode(Y) <- "double"

    storage.mode(X) <- "double"
    storage.mode(wt) <- "double"
    temp <- as.double(unlist(parms))
    if (length(temp)==0) temp <- 0.0    #if NULL, pass a dummy
    pred <- .Call(C_xpred,
                  ncat = as.integer(cats* !fit$ordered),
                  method= as.integer(method.int),
                  as.double(unlist(controls)),
                  temp,
                  as.integer(xval),
                  as.integer(xgroups),
                  Y,
                  X,
                  wt,
                  as.integer(numy),
                  as.double(costs),
                  as.integer(return.all),
                  as.double(cp),
                  as.double(fit$frame[1, "dev"]),
                  as.integer(numresp))

    if (return.all && numresp>1) {
        temp <- array(pred, dim=c(numresp, length(cp), nrow(X)),
                      dimnames=list(NULL, format(cp), dimnames(X)[[1]]))
        aperm(temp)                     # flip the dimensions
    }
    else matrix(pred, nrow=nrow(X), byrow=T,
                dimnames=list(dimnames(X)[[1]], format(cp)))
}
