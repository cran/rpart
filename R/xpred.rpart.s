# SCCS @(#)xpred.rpart.s	1.14 12/13/99
#
#  Get a set of cross-validated predictions
#

xpred.rpart <- function(fit, xval=10, cp) {
    if (!inherits(fit, 'rpart')) stop("Invalid fit object")

    method <- fit$method
    method.int <- pmatch(method, c("anova", "poisson", "class", "exp"))
    if (method.int==4) method.int <- 2
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
	    m <- eval(m, sys.frame(sys.parent()))
	    }
	if (is.null(X)) X <- rpart.matrix(m)
	if (is.null(wt)) wt <- model.extract(m, "weights")
	if (is.null(Y)) {
	    Y <- model.extract(m, "response")
            offset <- attr(Terms, "offset")
	    init <- (get(paste("rpart", method, sep='.')))(Y,offset, NULL)
	    Y <- init$y
	    }
	}

    Y <- as.matrix(Y)
    nobs <- nrow(Y)
    if (length(wt)==0) wt <- rep(1.0, nobs)

    cats <- rep(0, ncol(X))
    xlevels <- attr(fit, "xlevels")
    if (!is.null(xlevels)){
        cats[as.numeric(names(xlevels))] <- unlist(lapply(xlevels, length))
        }

    controls <- fit$control
    if (missing(cp)) {
	cp<- fit$cptable[,1]
	cp <- sqrt(cp * c(10, cp[-length(cp)]))
	}
    ncp <- length(cp)

    if (length(xval)==1) {
	# make random groups
	xgroups <- sample(rep(1:xval, length=nobs), nobs, replace=F)
	}
    else if (length(xval) == nrow(Y)) {
	xgroups <- xval
	xval <- length(unique(xgroups))
	}
    else stop("Invalid value for xval")

    rpfit <- .C("s_xpred",
		    n = as.integer(nobs),
		    nvarx = as.integer(ncol(X)),
		    ncat = as.integer(cats),
		    method= as.integer(method.int),
		    as.double(unlist(controls))[1:7],
		    parms = as.double(fit$parms),
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
		    NAOK=T )
    if (rpfit$n == -1)  stop(rpfit$error)

    matrix(rpfit$pred, ncol=ncp, byrow=T,
		dimnames=list(dimnames(X)[[1]], format(cp)) )
    }
