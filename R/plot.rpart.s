plot.rpart <- function(x, uniform=FALSE, branch=1, compress=FALSE,
			     nspace, margin=0, minbranch=.3, ...){
    if(!inherits(x, "rpart"))
	    stop("Not an rpart object")
    if (!is.null(x$frame$splits)) x <- rpconvert(x)  #help for old objects
    if (nrow(x$frame) <= 1)
        stop("fit is not a tree, just a root")

    if (compress & missing(nspace)) nspace <- branch
    if (!compress) nspace <- -1     #means no compression
    if (dev.cur() == 1) get(getOption("device"))()
    assign(paste(".rpart.parms", dev.cur(), sep = "."),
            list(uniform=uniform, branch=branch, nspace=nspace,
		 minbranch=minbranch), envir=.GlobalEnv)

    #define the plot region
    temp <- rpartco(x)
    xx <- temp$x
    yy <- temp$y
    temp1 <- range(xx) + diff(range(xx))*c(-margin, margin)
    temp2 <- range(yy) + diff(range(yy))*c(-margin, margin)
    plot(temp1, temp2, type='n', axes=FALSE, xlab='', ylab='', ...)

    # Draw a series of horseshoes or V's, left son, up, down to right son
    #   NA's in the vector cause lines() to "lift the pen"
    node <- as.numeric(row.names(x$frame))
    temp <- rpart.branch(xx, yy, node, branch)

    if (branch > 0) text(xx[1], yy[1], '|')
    lines(c(temp$x), c(temp$y))
    invisible(list(x=xx, y=yy))
}




