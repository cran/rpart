## SCCS %W% %G%
predict.rpart <-
function(object, newdata = list(),
	 type = c("vector", "matrix", "class", "probs"))
{
    if(!inherits(object, "rpart"))
	stop("Not legitimate tree")
    type <- match.arg(type)
    if(missing(newdata))
	where <- object$where
    else {
	if(is.null(attr(newdata, "terms"))) {
	    Terms <- delete.response(object$terms)
	    act <- (object$call)$na.action
	    if (is.null(act)) act<- na.rpart
	    newdata <- model.frame(Terms, newdata, na.action = act,
                                      xlev=attr(object, "xlevels"))
        }
	where <- pred.rpart(object, rpart.matrix(newdata))
    }
    frame <- object$frame
    method <- object$method
    ylevels <- attr(object, "ylevels")
    nclass <- length(ylevels)
    if(missing(type) && nclass > 0) type <- "probs"
    if(type == "vector" || (type=="matrix" && is.null(frame$yval2))) {
	pred <- frame$yval[where]
	names(pred) <- names(where)
    }
    else if (type == "matrix") {
	pred <- frame$yval2[where,]
	dimnames(pred) <- list(names(where), NULL)
    }
    else if(type == "class" && nclass > 0) {
	pred <- factor(ylevels[frame$yval[where]], levels=ylevels)
	names(pred) <- names(where)
    }
    else if (type == "probs" && nclass > 0) {
	pred <- frame$yval2[where, 1 + nclass + 1:nclass]
	dimnames(pred) <- list(names(where), ylevels)
    }
    else stop("Invalid prediction for rpart object")

    # Expand out the missing values in the result
    # But only if operating on the original dataset
    if (missing(newdata) && !is.null(object$na.action))
        pred <- naresid(object$na.action, pred)
    pred
}

