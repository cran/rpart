# SCCS  @(#)rpart.s	1.31 04/23/01
#
#  The recursive partitioning function, for S
#
rpart <- function(formula, data=NULL, weights, subset,
		   na.action=na.rpart, method, model=FALSE, x=FALSE, y=TRUE,
		   parms, control, ...) {


    call <- match.call()
    if (is.data.frame(model)) {
	m <- model
	model <- FALSE
	}
    else {
	m <- match.call(expand=FALSE)
	m$model <- m$method <- m$control<- NULL
	m$x <- m$y <- m$parms <- m$... <- NULL
	m$na.action <- na.action
	m[[1]] <- as.name("model.frame.default")
	m <- eval(m, parent.frame())
	}
    Terms <- attr(m, "terms")
    if(any(attr(Terms, "order") > 1))
	stop("Trees cannot handle interaction terms")

    Y <- model.extract(m, "response")
    wt <- model.extract(m, "weights")
    if(length(wt)==0) wt <- rep(1.0, nrow(m))
    offset <- attr(Terms, "offset")
    X <- rpart.matrix(m)
    nobs <- nrow(X)

    if (missing(method)) {
	if (is.factor(Y) || is.character(Y))      method <- 'class'
        else if (is.Surv(Y))   method <- 'exp'
	else if (is.matrix(Y)) method<- 'poisson'
	else                   method<- 'anova'
	}

    if (is.list(method)) {
	# User written split methods
	mlist <- method
	method <- 'user'

	if (missing(parms)) init <- mlist$init(Y, offset,,wt)
	else                init <- mlist$init(Y, offset, parms, wt)

	method.int <- 4      #the fourth entry in func_table.h
if(F) {
	if (!is.list(mlist) || length(mlist) !=3)
		stop("User written methods must have 3 functions")
	if (is.null(mlist$init) || typeof(mlist$init) != 'closure')
		stop("User written method does not contain an init function")
	if (is.null(mlist$split) || typeof(mlist$split) != 'closure')
		stop("User written method does not contain a split function")
	if (is.null(mlist$eval) || typeof(mlist$eval) != 'closure')
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
	    expr2 <- quote({
		temp <- user.eval(yback[1:nback], wback[1:nback], parms)
		if (length(temp$label) != numresp)
			stop("User eval function returned invalid label")
		if (length(temp$deviance) !=1)
			stop("User eval function returned invalid deviance")
		as.numeric(as.vector(c(temp$deviance, temp$label)))
		})
	    expr1 <- quote({
		if (nback <0) { #categorical variable
		    n2 <- -1*nback
		    temp  <- user.split(yback[1:n2], wback[1:n2],
					xback[1:n2], parms, FALSE)
		    ncat <- length(unique(xback[1:n2]))
		    if (length(temp$goodness) != ncat-1 ||
			length(temp$direction) != ncat)
			    stop("Invalid return from categorical split fcn")
		    }

		else {
		    temp <- user.split(yback[1:nback], wback[1:nback],
				       xback[1:nback], parms, TRUE)
		    if (length(temp$goodness) != (nback-1))
			stop("User split function returned invalid goodness")
		    if (length(temp$direction) != (nback-1))
			stop("User split function returned invalid direction")
		    }
		as.numeric(as.vector(c(temp$goodness, temp$direction)))
		})
	    }
	else {
	    expr2 <- quote({
		tempy <- matrix(yback[1:(nback*numy)], ncol=numy)
		temp <- user.eval(tempy, wback[1:nback], parms)
		if (length(temp$label) != numresp)
			stop("User eval function returned invalid label")
		if (length(temp$deviance !=1))
			stop("User eval function returned invalid deviance")
		as.numeric(as.vector(c(temp$deviance, temp$label)))
		})
	    expr1 <- quote({
		if (nback <0) { #categorical variable
		    n2 <- -1*nback
		    tempy <- matrix(yback[1:(n2*numy)], ncol=numy)
		    temp  <- user.split(tempy, wback[1:n2], xback[1:n2],
					parms, FALSE)
		    ncat <- length(unique(xback[1:n2]))
		    if (length(temp$goodness) != ncat-1 ||
			length(temp$direction) != ncat)
			    stop("Invalid return from categorical split fcn")
		    }
		else {
		    tempy <- matrix(yback[1:(nback*numy)], ncol=numy)
		    temp <- user.split(tempy, wback[1:nback],xback[1:nback],
				       parms, TRUE)
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
        rho <- new.env()
        assign("nback", integer(1), envir = rho)
        assign("wback", double(nobs), envir = rho)
        assign("xback", double(nobs), envir = rho)
        assign("yback", double(nobs), envir = rho)
        assign("user.eval", user.eval, envir = rho)
        assign("user.split", user.split, envir = rho)
        assign("numy", numy, envir = rho)
        assign("numresp", numresp, envir = rho)
        assign("parms", parms, envir = rho)
	.Call("init_rpcallback", rho, as.integer(numy),
	                         as.integer(numresp),
	                         expr1, expr2)
    }
        ## assign this to avoid garbage collection
        keep <- rpartcallback(mlist, nobs, init)
    }
    else {
	method.int <- pmatch(method, c("anova", "poisson", "class", "exp"))
	if (is.na(method.int)) stop("Invalid method")
	method <- c("anova", "poisson", "class", "exp")[method.int]
	if (method.int==4) method.int <- 2

	if (missing(parms))
	  init <- (get(paste("rpart", method, sep='.')))(Y,offset, ,wt)
	else
	  init <- (get(paste("rpart", method, sep='.')))(Y,offset, parms, wt)
	}

    Y <- init$y

    xlevels <- attr(X, "column.levels")
    cats <- rep(0,ncol(X))
    if(!is.null(xlevels)) {
	cats[match(names(xlevels), dimnames(X)[[2]])] <-
		  unlist(lapply(xlevels, length))
	}

    controls <- rpart.control(...)
    if (!missing(control)) controls[names(control)] <- control

    xval <- controls$xval
    if (is.null(xval) || (length(xval)==1 && xval==0) || method=='user') {
	xgroups <-0
	xval <- 0
	}
    else if (length(xval)==1) {
	# make random groups
        xgroups <- sample(rep(1:xval, length=nobs), nobs, replace=FALSE)
	}
    else if (length(xval) == nobs) {
	xgroups <- xval
	xval <- length(unique(xgroups))
	}
    else stop("Invalid value for xval")

    #
    # Have s_to_rp consider ordered categories as continuous
    #  A right-hand side variable that is a matrix forms a special case
    # for the code.
    #
    tfun <- function(x) {
	if (is.matrix(x)) rep(is.ordered(x), ncol(x))
	else is.ordered(x)
	}
    isord <- unlist(lapply(m[attr(Terms, 'term.labels')], tfun))
    rpfit <- .C("s_to_rp",
		    n = as.integer(nobs),
		    nvarx = as.integer(ncol(X)),
		    ncat = as.integer(cats* !isord),
		    method= as.integer(method.int),
		    as.double(unlist(controls)),
		    parms = as.double(init$parms),
		    as.integer(xval),
		    as.integer(xgroups),
		    as.double(t(init$y)),
		    as.double(X),
		    as.integer(!is.finite(X)),
		    error = character(1),
		    wt = as.double(wt),
		    as.integer(init$numy),
		    NAOK=TRUE )
    if (rpfit$n == -1)  stop(rpfit$error)

    # rpfit$newX[1:n] contains the final sorted order of the observations
    nodes <- rpfit$n          # total number of nodes
    nsplit<- rpfit$nvarx      # total number of splits, primary and surrogate
    numcp <- rpfit$method     # number of lines in cp table
    ncat  <- rpfit$ncat[1]    #total number of categorical splits
    numresp<- init$numresp    # length of the response vector

    if (nsplit==0) stop("No splits found")
    cpcol <- if (xval>0) 5 else 3
    if (ncat==0) catmat <- 0
    else         catmat <- matrix(integer(1), ncat, max(cats))

    rp    <- .C("s_to_rp2",
		       as.integer(nobs),
		       as.integer(nsplit),
		       as.integer(nodes),
		       as.integer(ncat),
		       as.integer(cats *!isord),
		       as.integer(max(cats)),
		       as.integer(xval),
		       which = integer(nobs),
		       cptable = matrix(double(numcp*cpcol), nrow=cpcol),
		       dsplit =  matrix(double(1),  nsplit,3),
		       isplit =  matrix(integer(1), nsplit,3),
		       csplit =  catmat,
		       dnode  =  matrix(double(1),  nodes, 3+numresp),
		       inode  =  matrix(integer(1), nodes, 6))
    tname <- c("<leaf>", dimnames(X)[[2]])

    if (cpcol==3) temp <- c("CP", "nsplit", "rel error")
    else          temp <- c("CP", "nsplit", "rel error", "xerror", "xstd")
    dimnames(rp$cptable) <- list(temp, 1:numcp)

    splits<- matrix(c(rp$isplit[,2:3], rp$dsplit), ncol=5,
		     dimnames=list(tname[rp$isplit[,1]+1],
			  c("count", "ncat", "improve", "index", "adj")))
    index <- rp$inode[,2]  #points to the first split for each node

    # Now, make ordered categories look like categories again (a printout
    #  choice)
    nadd <- sum(isord[rp$isplit[,1]])
    if (nadd >0) {
	newc <- matrix(integer(1), nadd, max(cats))
	cvar <- rp$isplit[,1]
	indx <- isord[cvar]		     # vector of T/F
	cdir <- splits[indx,2]               # which direction splits went
	ccut <- floor(splits[indx,4])        # cut point
	splits[indx,2] <- cats[cvar[indx]]   #Now, # of categories instead
	splits[indx,4] <- ncat + 1:nadd      # rows to contain the splits

	# Next 4 lines can be done without a loop, but become indecipherable
	for (i in 1:nadd) {
	    newc[i, 1:(cats[(cvar[indx])[i]])] <- -1*as.integer(cdir[i])
	    newc[i, 1:ccut[i]] <- as.integer(cdir[i])
	    }
	if (ncat==0) catmat <- newc
	else         catmat <- rbind(rp$csplit, newc)
	ncat <- ncat + nadd
	}
    else catmat <- rp$csplit

    temp <- ifelse(index==0, 1, index)
    svar <- ifelse(index==0, 0, rp$isplit[temp,1]) #var number for each node
    frame <- data.frame(row.names=rp$inode[,1],
			   var=  factor(svar, 0:ncol(X), tname),
			   n =   rp$inode[,5],
			   wt=   rp$dnode[,3],
			   dev=  rp$dnode[,1],
			   yval= rp$dnode[,4],
			   complexity=rp$dnode[,2],
			   ncompete  = pmax(0, rp$inode[,3]-1),
			   nsurrogate=rp$inode[,4])

    if (method.int ==3 ) {
        numclass <- init$numresp -1
        temp <- rp$dnode[,-(1:4)] %*% diag(init$parms[1:numclass]*
					   sum(init$counts)/init$counts)
        yprob <- matrix(temp /rp$dnode[,3] ,ncol=numclass)
        yval2 <- matrix(rp$dnode[, -(1:3)], ncol=numclass+1)
	frame$yval2 <- cbind(yval2, yprob)
	}
    else if (init$numy >1) frame$yval2 <- rp$dnode[,-(1:3)]

    if (is.null(init$summary))
	    stop("Initialization routine is missing the summary function")
    if (is.null(init$print))
	    functions <- list(summary=init$summary)
    else    functions <- list(summary=init$summary, print=init$print)
    if (!is.null(init$text)) functions <- c(functions, list(text=init$text))
    if (method=='user')	functions <- c(functions, mlist)

    ans <- list(frame = frame,
                where = structure(rp$which, names = row.names(m)),
                call=call, terms=Terms,
    		cptable =  t(rp$cptable),
		splits = splits,
		method = method,
		parms  = init$parms,
		control= controls,
		functions= functions)

    if (ncat>0) ans$csplit <- catmat +2
    if (model) {
	ans$model <- m
	if (missing(y)) y <- FALSE
	}
    if (y) ans$y <- Y
    if (x) {
	ans$x <- X
	ans$wt<- wt
	}
    ans$ordered <- isord
    if (!is.null(xlevels)) attr(ans, 'xlevels') <- xlevels
    if(method=='class') attr(ans, "ylevels") <- init$ylevels
    na.action <- attr(m, "na.action")
    if (length(na.action)) ans$na.action <- na.action
    class(ans) <- c("rpart")
    ans
    }
