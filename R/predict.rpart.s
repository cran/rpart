## SCCS @(#)predict.rpart.s	1.9 02/15/00
predict.rpart <-
function(object, newdata = list(), type = c("vector", "tree", "class"))
{
    if(!inherits(object, "rpart"))
	stop("Not legitimate tree")
    type <- match.arg(type)
    if(missing(newdata) & type == "tree")
	return(object)  #idiot proofing
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
    ylevels <- attr(object,'ylevels')
    if(type == "vector") {
      if(length(ylevels)>0){
	frame <- frame$yprob[where,]
	dimnames(frame)[[1]] <- names(where)
      } else {
	frame <- frame$yval[where]
	names(frame) <- names(where)
      }
      return(frame)
    } else if(type == "class") {
      if(length(ylevels) == 0)
	stop("Type class is only appropriate for classification")
      frame <- factor(ylevels[frame$yval[where]], levels=ylevels)
      names(frame) <- names(where)
      return(frame)
    } else stop("Cannot do rpart objects yet")
}

