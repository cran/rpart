#SCCS  @(#)residuals.rpart.s	1.7 02/11/00


residuals.rpart <- function(object, type)
    {

    if(!inherits(object, "rpart"))
	    stop("Not legitimate rpart object")

    if (object$method=='anova' || object$method=='class')
      { ## code taken directly from residuals.tree
        if(is.null(y <- object$y))
                y <- model.extract(model.frame(object), "response")
        frame <- object$frame

        if(is.null(ylevels <- attr(object, "ylevels"))) {
                tmpr <- y - frame$yval[object$where]    #       y <- unclass(y)
		#Expand out the missing values in the result
                if (!is.null(object$na.action))
               	   tmpr <- naresid(object$na.action, tmpr)
                return(tmpr)
	      }

        if(missing(type))
                type <- "usual"
        else if(is.na(match(type, c("usual", "pearson", "deviance"))))
                stop("Don't know about this type of residual")
        if(type == "usual")
                yhat <- frame$yval[object$where]
        else {
            nclass <- length(ylevels)
            yprob  <- frame$yval2[, 1+nclass + 1:nclass]
            yhat <- yprob[object$where,  ][cbind(seq(y), unclass(y))]
        }
        r <- switch(type,
                usual = as.integer(y != yhat),
                # misclassification
                pearson = (1 - yhat)/yhat,
                # sum((obs-fitted)/fitted)
                deviance = -2 * log(yhat))
        names(r) <- names(y)
       }

    else {
	if(is.null(y <- object$y))
		y <- model.extract(model.frame(object), "response")
	lambdat  <- (object$frame$yval)[object$where] * y[,1]

	events <- y[,2]
	temp <- pmax(events, 1)
	r <- sign(events-lambdat) *
		  sqrt(-2*((events - lambdat) + events*log(lambdat/temp)))
	}

    #Expand out the missing values in the result
    if (!is.null(object$na.action))
	r <- naresid(object$na.action, r)

    r
    }
