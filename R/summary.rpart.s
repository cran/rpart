#SCCS  %W% %G%
summary.rpart <- function(object, cp=0, digits=getOption("digits"), file,  ...)
{
    if(!inherits(object, "rpart")) stop("Not legitimate rpart object")

     if (!missing(file)) {
	  sink(file)
	  on.exit(sink())
	  }

    if(!is.null(object$call)) {
        cat("Call:\n")
        dput(object$call)
        }

     omit <- object$na.action
     n <- object$frame$n
     if (length(omit))
      cat("  n=", n[1], " (", naprint(omit), ")\n\n", sep="")
     else cat("  n=", n[1], "\n\n")

     print(object$cptable, digits=digits)
     ff <- object$frame
     ylevel <- attr(object,'ylevels')
     id <- as.integer(row.names(ff))
     parent.id <- ifelse(id==1,1, floor(id/2))
     parent.cp <- ff$complexity[match(parent.id, id)]
     rows <- (1:length(id))[parent.cp > cp]
     if (length(rows)>0) rows <- rows[order(id[rows])]
     else rows <- 1
     is.leaf <- (ff$var=='<leaf>')
     index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))

     sname <- dimnames(object$splits)[[1]]
     cuts <- vector(mode='character', length=nrow(object$splits))
     temp <- object$splits[ ,2]
     for (i in 1:length(cuts)) {
	 if (temp[i] == -1)
	     cuts[i] <-paste("<", format(signif(object$splits[i,4], digits=digits)))
	 else if (temp[i] ==1)
	     cuts[i] <-paste("<", format(signif(object$splits[i,4], digits=digits)))
	 else cuts[i]<- paste("splits as ",
	     paste(c("L", "-", "R")[object$csplit[object$splits[i,4], 1:temp[i]]],
                   collapse='', sep=''), collapse='')
	 }
    # S-PLUS 4.0 can't handle null vectors here
    if(any(temp<2)) cuts[temp<2 ] <- format(cuts[temp<2],justify="left")
    cuts <- paste(cuts, ifelse(temp >=2, ",",
			 ifelse(temp==1, " to the right,", " to the left, ")),
			 sep = '')

    if (is.null(ff$yval2))
	    tprint <- object$functions$summary(ff$yval[rows], ff$dev[rows],
					  ff$wt[rows], ylevel, digits)
    else
	    tprint <- object$functions$summary(ff$yval2[rows,], ff$dev[rows],
					  ff$wt[rows], ylevel, digits)

    for (ii in 1:length(rows)) {
	i <- rows[ii]
	nn <- ff$n[i]
	twt <- ff$wt[i]
	cat("\nNode number ", id[i], ": ", nn, " observations", sep='')
	if (ff$complexity[i] < cp || is.leaf[i]) cat("\n")
	else cat(",    complexity param=",
		       format(signif(ff$complexity[i], digits)), "\n", sep="")

	cat(tprint[ii], "\n")
	if (ff$complexity[i] > cp && !is.leaf[i] ){
	    sons <- 2*id[i] + c(0,1)
	    sons.n <- ff$n[match(sons, id)]
	    cat("  left son=", sons[1], " (", sons.n[1], " obs)",
		" right son=", sons[2], " (", sons.n[2], " obs)", sep='')
	    j <- nn - (sons.n[1] + sons.n[2])
	    if (j>1) cat(", ", j, " observations remain\n", sep='')
	    else if (j==1) cat(", 1 observation remains\n")
	    else     cat("\n")
	    cat("  Primary splits:\n")
	    j <- seq(index[i], length=1+ff$ncompete[i])
	    if (all(nchar(cuts[j]) < 25))
		  temp <- format(cuts[j], justify="left")
	    else  temp <- cuts[j]
	    cat(paste("      ", format(sname[j], justify="left"), " ", temp,
		      " improve=", format(signif(object$splits[j,3], digits)),
		      ", (", nn - object$splits[j,1], " missing)", sep=''),
		      sep="\n")
	    if (ff$nsurrogate[i] >0) {
		cat("  Surrogate splits:\n")
		j <- seq(1 +index[i] + ff$ncompete[i], length=ff$nsurrogate[i])
		agree <- object$splits[j,3]
		adj   <- object$splits[j,5]
		if (all(nchar(cuts[j]) < 25))
		      temp <- format(cuts[j], justify="left")
		else  temp <- cuts[j]
		cat(paste("      ", format(sname[j], justify="left"), " ",
		      temp,
		      " agree=", format(round(agree, 3)),
                      ", adj=" , format(round(adj, 3)),
		      ", (", object$splits[j,1], " split)", sep=''),
		      sep="\n")
		}
	    }
	}
    cat("\n")
    invisible(object)
    }
