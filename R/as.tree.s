#SCCS @(#)as.tree.s	1.2 03/21/97
# Change an rpart object into a tree object

as.tree <- function(x) {
    if (!inherits(x, 'rpart')) stop ("Only applicable to rpart objects")
  
    temp <- match(c('var', 'n', 'dev', 'yval', 'splits', 'yprob'), 
		          dimnames(x$frame)[[2]], nomatch=0)
    x$frame <- x$frame[,temp]
    x$cptable <- NULL
    x$splits <- NULL
    x$parms <- NULL
    x$control <- NULL
    x$csplit <- NULL
    
    # Fix xlevels, the only hard part
    xlevs <- attr(x, 'xlevels')
    nclass <- length(xlevs)
    xvars <- attr(x$terms, 'term.labels')
    nvar <- length(xvars)

    if (length(xlevs)==0) newlev <- vector('list',nvar)  #no factors
    else {
	temp <- match(1:nvar, as.numeric(names(xlevs)))
	newlev <- xlevs[temp]
        }
    names(newlev) <- xvars
    attr(x, 'xlevels') <- newlev

    class(x) <- 'tree'
    x
    }

	
	
