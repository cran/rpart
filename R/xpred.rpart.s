# SCCS
#
#  Get a set of cross-validated predictions
xpred.rpart <- function(fit, xval=10, cp) {
    if (!inherits(fit, 'rpart')) stop("Invalid fit object")

    method <- fit$method
    method.int <- pmatch(method, c("anova", "poisson", "class", "user", "exp"))
    if (method.int==5) method.int <- 2
    Terms <- fit$terms

    Y <- fit$y
    X <- fit$x
    wt<- fit$wt
    if (is.null(Y) || is.null(X)) {
	m <- fit$model
	if (is.null(m)) {
	    m <-fit$call[match(c("", 'formula', 'data', 'weights', 'subset',
					 'na.action'),
				names(fit$call), nomatch=0)]
	    if (is.null(m$na.action)) m$na.action<- na.rpart
	    m[[1]] <- as.name("model.frame.default")
	    m <- eval(m, parent.frame())
	    }
	if (is.null(X)) X <- rpart.matrix(m)
	if (is.null(wt)) wt <- model.extract(m, "weights")
	if (is.null(Y)) {
	    yflag <- TRUE
	    Y <- model.extract(m, "response")
            offset <- attr(Terms, "offset")
	    if (method != user) {
		init <- (get(paste("rpart", method, sep='.')))(Y,offset, NULL)
		Y <- as.matrix(init$y)
		numy <- ncol(Y)
		}
	    }
	else {
	    yflag <- FALSE
	    Y <- as.matrix(Y)
	    numy <- ncol(Y)
	    offset <- 0
	    }
	}

    nobs <- nrow(X)
    if (length(wt)==0) wt <- rep(1.0, nobs)

    cats <- rep(0, ncol(X))
    xlevels <- attr(fit, "xlevels")
    if (!is.null(xlevels)){
        cats[match(names(xlevels), dimnames(X)[[2]])] <-
		 unlist(lapply(xlevels, length))
        }

    controls <- fit$control
    if (missing(cp)) {
	cp<- fit$cptable[,1]
	cp <- sqrt(cp * c(10, cp[-length(cp)]))
	cp[1] <- (1+fit$cptable[1,1])/2
	}
    ncp <- length(cp)

    if (length(xval)==1) {
	# make random groups
	xgroups <- sample(rep(1:xval, length=nobs), nobs, replace=FALSE)
	}
    else if (length(xval) == nrow(Y)) {
	xgroups <- xval
	xval <- length(unique(xgroups))
	}
    else stop("Invalid value for xval")

    parms <- as.double(fit$parms)
    if (method=='user') {
	mlist <- fit$functions

        if (!is.list(mlist) || length(mlist) <3)
		stop("User written methods must have 3 functions")
	if (is.null(mlist$init) || typeof(mlist$init) != 'closure')
		stop("User written method does not contain an init function")

	# I need to call init to find out numresp, even though I
	#   only intend to use the first element of the response vector
	#  (For the "built-in" routines numresp is hardcoded into C).
	# If yflag==F, the "transformed" Y is already in hand, so we don't
	#   need to replace it with init$y
	if (length(parms)==0) init <- mlist$init(Y, offset,,wt)
	else                  init <- mlist$init(Y, offset, parms, wt)

        keep <- rpartcallback(mlist, nobs, init)

        if(F) {
	numresp <- init$numresp
	numy <-  init$numy
	if (yflag)  Y <- init$y

	if (is.null(mlist$split) || class(mlist$split) != 'function')
		stop("User written method does not contain a split function")
	if (is.null(mlist$eval) || class(mlist$eval) != 'function')
		stop("User written method does not contain an eval function")

	user.eval <- mlist$eval
	user.split <- mlist$split

	numresp <- init$numresp
	numy <-  init$numy
	parms <- init$parms

	#
	# expr2 is an expression that will call the user "evaluation"
	#   function, and check that what comes back is valid
	# expr1 does the same for the user "split" function
	#
	# For speed in the C interface, yback, xback, and wback are
	#  fixed S vectors of a fixed size, and nback tells us how
	#  much of the vector is actually being used on this particular
	#  callback.
	#
	if (numy==1) {
	    expr2 <- Quote({
		temp <- user.eval(yback[1:nback], wback[1:nback], parms)
		if (length(temp$label) != numresp)
			stop("User eval function returned invalid label")
		if (length(temp$deviance) !=1)
			stop("User eval function returned invalid deviance")
		as.numeric(as.vector(c(temp$deviance, temp$label)))
		})
	    expr1 <- Quote({
		if (nback <0) { #categorical variable
		    n2 <- -1*nback
		    temp  <- user.split(yback[1:n2], wback[1:n2],
					xback[1:n2], parms, F)
		    ncat <- length(unique(xback[1:n2]))
		    if (length(temp$goodness) != ncat-1 ||
			length(temp$direction) != ncat)
			    stop("Invalid return from categorical split fcn")
		    }

		else {
		    temp <- user.split(yback[1:nback], wback[1:nback],
				       xback[1:nback], parms, T)
		    if (length(temp$goodness) != (nback-1))
			stop("User split function returned invalid goodness")
		    if (length(temp$direction) != (nback-1))
			stop("User split function returned invalid direction")
		    }
		as.numeric(as.vector(c(temp$goodness, temp$direction)))
		})
	    }
	else {
	    expr2 <- Quote({
		tempy <- matrix(yback[1:(nback*numy)], ncol=numy)
		temp <- user.eval(tempy, wback[1:nback], parms)
		if (length(temp$label) != numresp)
			stop("User eval function returned invalid label")
		if (length(temp$deviance !=1))
			stop("User eval function returned invalid deviance")
		as.numeric(as.vector(c(temp$deviance, temp$label)))
		})
	    expr1 <- Quote({
		if (nback <0) { #categorical variable
		    n2 <- -1*nback
		    tempy <- matrix(yback[1:(n2*numy)], ncol=numy)
		    temp  <- user.split(tempy, wback[1:n2], xback[1:n2],
					parms, F)
		    ncat <- length(unique(xback[1:n2]))
		    if (length(temp$goodness) != ncat-1 ||
			length(temp$direction) != ncat)
			    stop("Invalid return from categorical split fcn")
		    }
		else {
		    tempy <- matrix(yback[1:(nback*numy)], ncol=numy)
		    temp <- user.split(tempy, wback[1:nback],xback[1:nback],
				       parms, T)
		    if (length(temp$goodness) != (nback-1))
			stop("User split function returned invalid goodness")
		    if (length(temp$direction) != (nback-1))
			stop("User split function returned invalid direction")
		    }
		as.numeric(as.vector(c(temp$goodness, temp$direction)))
		})
	    }
	#
	# The vectors nback, wback, xback and yback will have their
	#  contents constantly re-inserted by C code.  It's one way to make
	#  things very fast.  It is dangerous to do this, so they
	#  are tossed into a separate frame to isolate them.  Evaluations of
	#  the above expressions occur in that frame.
	#
	eframe <- new.frame(list(nback = integer(1),
				 wback = double(nobs),
				 xback = double(nobs),
				 yback = double(nobs),
				 user.eval =  user.eval,
				 user.split = user.split,
				 numy = numy,
				 numresp = numresp,
				 parms = parms), protect=T)
	.Call("init_rpcallback", eframe, as.integer(numy),
	                         as.integer(numresp),
	                         expr1, expr2)
	}
    }

    rpfit <- .C("s_xpred",
		    n = as.integer(nobs),
		    nvarx = as.integer(ncol(X)),
		    ncat = as.integer(cats * !fit$ordered),
		    method= as.integer(method.int),
		    as.double(unlist(controls)),
		    parms = as.double(parms),
		    as.integer(xval),
		    as.integer(xgroups),
		    as.double(t(Y)),
		    as.double(X),
		    as.integer(is.na(X)),
		    pred = double(ncp* nobs),
		    as.integer(ncp),
		    as.double(cp * fit$frame[1,"dev"]),
		    error = character(1),
		    wt = as.double(wt),
		    as.integer(numy),
		    NAOK=TRUE )
    if (rpfit$n == -1)  stop(rpfit$error)

    matrix(rpfit$pred, ncol=ncp, byrow=TRUE,
		dimnames=list(dimnames(X)[[1]], format(cp)) )
    }
