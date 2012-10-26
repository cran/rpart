#
#  The recursive partitioning function, for R
#
rpart <- function(formula, data, weights, subset,
                  na.action = na.rpart, method,
                  model = FALSE, x = FALSE, y = TRUE,
                  parms, control, cost, ...)
{
    Call <- match.call()
    if (is.data.frame(model)) {
        m <- model
        model <- FALSE
    } else {
        indx <- match(c("formula", "data", "weights", "subset"),
                      names(Call), nomatch = 0L)
        if (indx[1] == 0L) stop("a 'formula' argument is required")
        temp <- Call[c(1L, indx)]        # only keep the arguments we wanted
        temp$na.action <- na.action    # This one has a default
        temp[[1L]] <- as.name('model.frame') # change the function called
        m <- eval(temp, parent.frame()) # or eval.parent
    }

    Terms <- attr(m, "terms")
    if(any(attr(Terms, "order") > 1L))
	stop("Trees cannot handle interaction terms")

    Y <- model.response(m)
    wt <- model.weights(m)
    if(length(wt) == 0L) wt <- rep(1.0, nrow(m))
    offset <- model.offset(m)
    X <- rpart.matrix(m)
    nobs <- nrow(X)
    nvar <- ncol(X)

    if (missing(method)) {
	method <- if (is.factor(Y) || is.character(Y)) "class"
        else if (inherits(Y, "Surv")) "exp"
	else if (is.matrix(Y)) "poisson"
	else "anova"
    }

    if (is.list(method)) {
        ## User-written split methods
	mlist <- method
	method <- "user"

        ## Set up C callback.  Assign the result to a variable to avoid
        ## garbage collection
	if (missing(parms)) init <- mlist$init(Y, offset, wt=wt)
	else                init <- mlist$init(Y, offset, parms, wt)
        keep <- rpartcallback(mlist, nobs, init)

	method.int <- 4L             #the fourth entry in func_table.h
	numresp <- init$numresp
	numy <-  init$numy
	parms <- init$parms
    }
    else {
	method.int <- pmatch(method, c("anova", "poisson", "class", "exp"))
	if (is.na(method.int)) stop("Invalid method")
	method <- c("anova", "poisson", "class", "exp")[method.int]
	if (method.int == 4L) method.int <- 2L

	if (missing(parms))
            init <- (get(paste("rpart", method, sep='.')))(Y, offset, ,wt)
	else
            init <- (get(paste("rpart", method, sep='.')))(Y, offset, parms, wt)
        ns <- asNamespace("rpart")
        if(!is.null(init$print)) environment(init$print) <- ns
        if(!is.null(init$summary)) environment(init$summary) <- ns
        if(!is.null(init$text)) environment(init$text) <- ns
    }

    Y <- init$y

    xlevels <- attr(X, "column.levels")
    cats <- rep(0L, ncol(X))
    if(!is.null(xlevels)) {
	cats[match(names(xlevels), dimnames(X)[[2L]])] <-
            unlist(lapply(xlevels, length))
    }

    ## We want to pass any ... args to rpart.control, but not pass things
    ##  like "dats=mydata" where someone just made a typo.  The use of ...
    ##  is simply to allow things like "cp=.05" with easier typing
    extraArgs <- list(...)
    if (length(extraArgs)) {
	controlargs <- names(formals(rpart.control)) #legal arg names
	indx <- match(names(extraArgs), controlargs, nomatch=0L)
	if (any(indx==0L))
            stop(gettextf("Argument %s not matched", names(extraArgs)[indx==0L]),
                 domain = NA)
    }

    controls <- rpart.control(...)
    if (!missing(control)) controls[names(control)] <- control

    xval <- controls$xval
    if (is.null(xval) || (length(xval) == 1L && xval == 0L) || method=="user") {
	xgroups <- 0L
	xval <- 0L
    }
    else if (length(xval) == 1L) {
                                        # make random groups
        xgroups <- sample(rep(1L:xval, length=nobs), nobs, replace=FALSE)
    }
    else if (length(xval) == nobs) {
	xgroups <- xval
	xval <- length(unique(xgroups))
    }
    else {
        ## Check to see if observations were removed due to missing
	if (!is.null(attr(m, "na.action"))) {
            ## if na.rpart was used, then na.action will be a vector
	    temp <- as.integer(attr(m, "na.action"))
	    xval <- xval[-temp]
	    if (length(xval) == nobs) {
		xgroups <- xval
		xval <- length(unique(xgroups))
            }
	    else stop("Wrong length for 'xval'")
        }
	else stop("Wrong length for 'xval'")
    }

    ##
    ## Incorporate costs
    ##
    if (missing(cost)) cost <- rep(1.0, nvar)
    else {
	if (length(cost) != nvar)
            stop("Cost vector is the wrong length")
	if (any(cost <= 0)) stop("Cost vector must be positive")
    }

    ##
    ## Have C code consider ordered categories as continuous
    ##  A right-hand side variable that is a matrix forms a special case
    ## for the code.
    ##
    tfun <- function(x) {
	if (is.matrix(x)) rep(is.ordered(x), ncol(x))
	else is.ordered(x)
    }
    labs <- sub("^`(.*)`$", "\\1", attr(Terms, 'term.labels')) #beware backticks
    isord <- unlist(lapply(m[labs], tfun))

    storage.mode(X) <- "double"
    storage.mode(wt) <- "double"
    temp <- as.double(unlist(init$parms))
    if (length(temp)==0) temp <- 0.0    #if parms is NULL pass a dummy
    rpfit <- .Call(C_rpart,
                   ncat = as.integer(cats * !isord),
                   method = as.integer(method.int),
                   as.double(unlist(controls)),
                   temp,
                   as.integer(xval),
                   as.integer(xgroups),
                   as.double(t(init$y)),
                   X,
                   wt,
                   as.integer(init$numy),
                   as.double(cost))

    nsplit <- nrow(rpfit$isplit) # total number of splits, primary and surrogate
    ncat <- if (!is.null(rpfit$csplit))
        nrow(rpfit$csplit) # total number of categorical splits
    else 0L
    nodes <- nrow(rpfit$inode)
    if (nsplit == 0L) xval <- 0L # No xvals were done if no splits were found

    numcp <- ncol(rpfit$cptable)
    temp <- if (nrow(rpfit$cptable) == 3L) c("CP", "nsplit", "rel error")
    else  c("CP", "nsplit", "rel error", "xerror", "xstd")
    dimnames(rpfit$cptable) <- list(temp, 1:numcp)

    tname <- c("<leaf>", dimnames(X)[[2]])
    splits<- matrix(c(rpfit$isplit[, 2:3], rpfit$dsplit), ncol = 5L,
                    dimnames=list(tname[rpfit$isplit[, 1L] + 1L],
                    c("count", "ncat", "improve", "index", "adj")))
    index <- rpfit$inode[,2]  #points to the first split for each node

    ## Now, make ordered categories look like categories again (a printout
    ##  choice)
    nadd <- sum(isord[rpfit$isplit[, 1L]])
    if (nadd > 0L) {
	newc <- matrix(1L, nadd, max(cats))
	cvar <- rpfit$isplit[, 1L]
	indx <- isord[cvar]             # vector of TRUE/FALSE
	cdir <- splits[indx, 2L]        # which direction splits went
	ccut <- floor(splits[indx, 4L]) # cut point
	splits[indx, 2L] <- cats[cvar[indx]] # Now, # of categories instead
	splits[indx, 4L] <- ncat + 1L:nadd # rows to contain the splits

        ## Next 4 lines can be done without a loop, but become indecipherable
	for (i in 1L:nadd) {
	    newc[i, 1L:(cats[(cvar[indx])[i]])] <- -1*as.integer(cdir[i])
	    newc[i, 1L:ccut[i]] <- as.integer(cdir[i])
        }
	catmat <- if (ncat == 0L) newc
        else {
            ## newc have more cols than existing categorical splits
            cs <- rpfit$csplit
            ncs <- ncol(cs); ncc <- ncol(newc)
            if (ncs < ncc) cs <- cbind(cs, matrix(1L, nrow(cs), ncc - ncs))
            rbind(cs, newc)
        }
	ncat <- ncat + nadd
    }
    else catmat <- rpfit$csplit

    if (nsplit == 0L) {                    #tree with no splits
	frame <- data.frame(row.names = 1L,
			    var =  "<leaf>",
			    n = rpfit$inode[, 5L],
			    wt = rpfit$dnode[, 3L],
			    dev =  rpfit$dnode[,1L],
			    yval = rpfit$dnode[,4L],
			    complexity =rpfit$dnode[,2L],
			    ncompete  = 0L,
			    nsurrogate = 0L)
    } else {
	temp <- ifelse(index == 0, 1, index)
	svar <- ifelse(index == 0, 0, rpfit$isplit[temp,1L]) #var number
	frame <- data.frame(row.names=rpfit$inode[,1],
			    var =  factor(svar, 0:ncol(X), tname),
			    n =   rpfit$inode[, 5L],
			    wt =   rpfit$dnode[, 3L],
			    dev =  rpfit$dnode[, 1L],
			    yval = rpfit$dnode[, 4L],
			    complexity = rpfit$dnode[, 2L],
			    ncompete = pmax(0L, rpfit$inode[, 3L]-1L),
			    nsurrogate = rpfit$inode[, 4L])
    }
    if (method.int == 3L) {
        ## Create the class probability vector from the class counts, and
        ##   add it to the results
        ## Also scale the P(T) result
        ## The "pmax" 3 lines down is for the case of a factor y which has
        ##   no one at all in one of its classes.  Both the prior and the
        ##   count will be zero, which led to a 0/0.
        numclass <- init$numresp -2L
        nodeprob <- rpfit$dnode[, numclass+5L] / sum(wt) # see ginidev.c
        temp <- pmax(1L,init$counts)    #overall class freq in data
        temp <- rpfit$dnode[, 4L+(1L:numclass)] %*% diag(init$parms$prior/temp)
        yprob <- temp /rowSums(temp)    #necessary with altered priors
        yval2 <- matrix(rpfit$dnode[, 4L+(0L:numclass)], ncol=numclass+1)
	frame$yval2 <- cbind(yval2, yprob, nodeprob)
    }
    else if (init$numresp > 1L)
        frame$yval2 <- rpfit$dnode[,-(1L:3L), drop=FALSE]

    if (is.null(init$summary))
        stop("Initialization routine is missing the 'summary' function")
    functions <- if (is.null(init$print)) list(summary = init$summary)
    else list(summary = init$summary, print = init$print)
    if (!is.null(init$text)) functions <- c(functions, list(text = init$text))
    if (method == "user") functions <- c(functions, mlist)

    where <- rpfit$which
    names(where) <- row.names(m)

    ## FIMXE clean up
    ans <- if (nsplit == 0L) {                  # no 'splits' component
	list(frame = frame,
             where = where,
             call = Call, terms = Terms,
             cptable =  t(rpfit$cptable),
             method = method,
             parms  = init$parms,
             control = controls,
             functions = functions,
             numresp = init$numresp)
    } else {
	list(frame = frame,
             where = where,
             call = Call, terms = Terms,
             cptable = t(rpfit$cptable),
             splits = splits,
             method = method,
             parms  = init$parms,
             control = controls,
             functions = functions,
             numresp = init$numresp)
    }
    if (ncat > 0L) ans$csplit <- catmat + 2L
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
    if(!is.null(attr(m, "na.action"))) ans$na.action <- attr(m, "na.action")
    if (!is.null(xlevels)) attr(ans, "xlevels") <- xlevels
    if (method == "class") attr(ans, "ylevels") <- init$ylevels
    class(ans) <- "rpart"
    ans
}
