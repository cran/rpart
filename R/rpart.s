# SCCS @(#)rpart.s	1.23 02/12/98
#
#  The recursive partitioning function, for S
#
rpart <- function(formula, data=NULL, weights, subset,
		   na.action=na.rpart, method, model=F, x=F, y=T,
		   parms, control, ...) {


    call <- match.call()
    if (is.data.frame(model)) {
	m <- model
	model <- F
	}
    else {
	m <- match.call(expand=F)
	m$model <- m$method <- m$control<- NULL
	m$x <- m$y <- m$parms <- m$... <- NULL
	m$na.action <- na.action
	m[[1]] <- as.name("model.frame.default")
	m <- eval(m, sys.frame(sys.parent()))
	}
    Terms <- attr(m, "terms")
    if(any(attr(Terms, "order") > 1))
	stop("Trees cannot handle interaction terms")

    Y <- model.extract(m, "response")
    if (missing(method)) {
	if (is.factor(Y))      method <- 'class'
        else if (is.Surv(Y))   method <- 'exp'
	else if (is.matrix(Y)) method<- 'poisson'
	else                   method<- 'anova'
	}
    method.int <- pmatch(method, c("anova", "poisson", "class", "exp"))
    if (is.na(method.int)) stop("Invalid method")
    method <- c("anova", "poisson", "class", "exp")[method.int]
    if (method.int==4) method.int <- 2

    w <- model.extract(m, "weights")
    if(length(w)) warning("Weights ignored")
    offset <- attr(Terms, "offset")

    if (missing(parms))
	  init <- (get(paste("rpart", method, sep='.')))(Y,offset)
    else
	  init <- (get(paste("rpart", method, sep='.')))(Y,offset, parms)
    Y <- init$y
    X <- rpart.matrix(m)
    nobs <- nrow(X)

    xlevels <- attr(X, "column.levels")
    cats <- rep(0,ncol(X))	
    if(!is.null(xlevels)) {
	cats[as.numeric(names(xlevels))] <- unlist(lapply(xlevels, length))
	}

    controls <- rpart.control(...)
    if (!missing(control)) controls[names(control)] <- control

    xval <- controls$xval
    if (is.null(xval) || (length(xval)==1 && xval==0)) {
	xgroups <-0
	xval <- 0
	}
    else if (length(xval)==1) {
	# make random groups
        xgroups <- sample(rep(1:xval, length=nobs), nobs, replace=F)
	}
    else if (length(xval) == nobs) {
	xgroups <- xval
	xval <- length(unique(xgroups))
	}
    else stop("Invalid value for xval")

    # 
    # Have s_to_rp consider ordered categories as continuous
    #
    isord <- sapply(m, is.ordered)[-1]
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
		    as.integer(is.na(X)),
		    error = character(1),
		    NAOK=T )
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
		       dsplit =  matrix(double(1),  nsplit,2),
		       isplit =  matrix(integer(1), nsplit,3),
		       csplit =  catmat,
		       dnode  =  matrix(double(1),  nodes, 2+numresp),
		       inode  =  matrix(integer(1), nodes, 6))
    tname <- c("<leaf>", dimnames(X)[[2]])

    if (cpcol==3) temp <- c("CP", "nsplit", "rel error")
    else          temp <- c("CP", "nsplit", "rel error", "xerror", "xstd")
    dimnames(rp$cptable) <- list(temp, 1:numcp)

    splits<- matrix(c(rp$isplit[,2:3], rp$dsplit), ncol=4,
		     dimnames=list(tname[rp$isplit[,1]+1],
			      c("count", "ncat", "improve", "index")))
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

    rplab <- .C("rplabel", as.integer(nodes),
			   as.integer(index),
			   splits[index, c(2,4)],
			   as.integer(ncat),
			   as.integer(catmat),
			   cutleft = character(nodes),
			   cutright= character(nodes))[6:7]

    temp <- ifelse(index==0, 1, index)
    svar <- ifelse(index==0, 0, rp$isplit[temp,1]) #var number for each node
    frame <- data.frame(row.names=rp$inode[,1],
			   var=  factor(svar, 0:ncol(X), tname),
			   n =   rp$inode[,5],
			   dev=  rp$dnode[,1],
			   yval= rp$dnode[,3],
			   complexity=rp$dnode[,2],
			   ncompete  = pmax(0, rp$inode[,3]-1),
			   nsurrogate=rp$inode[,4])

    frame$splits <- matrix(unlist(rplab), ncol=2,
			   dimnames=list(NULL, c("cutleft", "cutright")))
    if (method=='class') {
        numclass <- init$numresp -1
        temp <- rp$dnode[,-(1:3)] %*% diag(init$parms[1:numclass]*nobs /
						 init$counts)
        frame$yprob <- matrix(temp /c(temp %*% rep(1,numclass)) ,
			   ncol=numclass, dimnames=list(NULL, init$ylevels))
        frame$yval2 <- matrix(rp$dnode[, -(1:3)], ncol=numclass,
        		    dimnames=list(NULL, init$ylevels))
	}	
    else if (method=='poisson' | method=='exp') frame$yval2 <- rp$dnode[,4]

    ans <- list(frame = frame, 
                where = structure(rp$which, names = row.names(m)),
                call=call, terms=Terms, 
    		cptable =  t(rp$cptable),
		splits = splits,
		method = method,
		parms  = init$parms,
		control= controls)

    if (ncat>0) ans$csplit <- catmat +2
    if (model) {
	ans$model <- m
	if (missing(y)) y <- F
	}
    if (y) ans$y <- Y
    if (x) ans$x <- X
    ans$control <- controls
    if (!is.null(xlevels)) attr(ans, 'xlevels') <- xlevels
    if(method=='class') attr(ans, "ylevels") <- init$ylevels
    class(ans) <- c("rpart")
    ans
    }
