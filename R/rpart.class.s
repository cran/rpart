#SCCS @(#)rpart.class.s	1.3 12/13/99
rpart.class <- function(y, offset, parms, wt) {
    if (!is.null(offset)) stop("No offset allowed in classification models")
    fy <- as.factor(y)
    y <- as.integer(fy)
    numclass <- max(y[!is.na(y)])
    counts <- tapply(wt, factor(y, levels=1:numclass), sum)
    counts <- ifelse(is.na(counts), 0, counts)   #in case of an empty class
    numresp <- 1+numclass
    if (missing(parms) || is.null(parms))
	parms <- c(counts/sum(counts), rep(1,numclass^2)-diag(numclass),1)
    else if (is.list(parms)) {
	if (is.null(parms$prior)) temp <- c(counts/sum(counts))
	else {
	    temp <- parms$prior
	    if (sum(temp) !=1) stop("Priors must sum to 1")
	    if (any(temp<0)) stop("Priors must be >= 0")
	    if (length(temp) != numclass) stop("Wrong length for priors")
	    }

	if (is.null(parms$loss)) temp2<- 1 - diag(numclass)
	else {
	    temp2 <- parms$loss
	    if (length(temp2) != numclass^2)
			    stop("Wrong length for loss matrix")
	    temp2 <- matrix(temp2, ncol=numclass)
	    if (any(diag(temp2)!=0))
			stop("Loss matrix must have zero on diagonals")
	    if (any(temp2<0))
			stop("Loss matrix cannot have negative elements")
	    if (any(apply(temp2,1,sum)==0))
			stop("Loss matrix has a row of zeros")
	    }

	if (is.null(parms$split)) temp3 <- 1
 	    else {
		temp3 <- pmatch(parms$split, c("gini", "information"))
		if (is.null(temp3)) stop("Invalid splitting rule")
		}
	parms <- c(temp, temp2, temp3)
	}
    else stop("Parameter argument must be a list")

    list(y=y, parms=parms, numresp=numclass+1, counts=counts,
		ylevels= levels(fy))
    }

