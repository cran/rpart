#SCCS  @(#)summary.rpart.s	1.9 02/12/98
summary.rpart <- function(x, cp=0, digits=.Options$digits-3,
			file,  ...) {
    if(!inherits(x, "rpart")) stop("Not legitimate rpart object")

     if (!missing(file)) {
	  sink(file)
	  on.exit(sink())
	  }

    if(!is.null(x$call)) {
        cat("Call:\n")
        dput(x$call)
        cat('\n')
        }

     print(x$cptable, digits=digits)
     ff <- x$frame
     ylevel <- attr(x,'ylevels')
     id <- as.integer(row.names(ff))
     parent.id <- ifelse(id==1,1, floor(id/2))
     parent.cp <- ff$complexity[match(parent.id, id)]
     rows <- (1:length(id))[parent.cp > cp]
     if (length(rows)>0) rows <- rows[order(id[rows])]
     else rows <- 1
     is.leaf <- (ff$var=='<leaf>')
     index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))

     sname <- dimnames(x$splits)[[1]]
     cuts <- vector(mode='character', length=nrow(x$splits))
     temp <- x$splits[ ,2]
     for (i in 1:length(cuts)) {
	 if (temp[i] == -1)
	     cuts[i] <-paste("<", format(signif(x$splits[i,4], digits=digits)))
	 else if (temp[i] ==1)
	     cuts[i] <-paste("<", format(signif(x$splits[i,4], digits=digits)))
	 else cuts[i]<- paste("splits as ",
	     paste(c("L", "-", "R")[x$csplit[x$splits[i,4], 1:temp[i]]],
                   collapse='', sep=''), collapse='')
	 }
    # S-PLUS 4.0 can't handle null vectors here
    if(any(temp<2)) cuts[temp<2 ] <- format(cuts[temp<2])
     cuts <- paste(cuts, ifelse(temp >=2, ",",
			 ifelse(temp==1, " to the right,", " to the left, ")),
			 sep = '')
     for (i in rows) {
	nn <- ff$n[i]
	cat("\nNode number ", id[i], ": ", nn, " observations", sep='')
	if (ff$complexity[i] < cp || is.leaf[i]) cat("\n")
	else cat(",    complexity param=",
		       format(signif(ff$complexity[i], digits)), "\n", sep="")

	if (x$method=='anova')
	    cat("  mean=", format(signif(ff$yval[i], digits)),
		" , SS/n=" , format(signif(ff$dev[i]/nn, digits)),"\n",
		sep = "")
	else if (x$method=='class') {
            if(!is.null(ylevel)) 
	       yval <- ylevel[ff$yval]
	    else
	       yval <- ff$yval
	    cat("  predicted class=", format(yval[i]),
		" expected loss=", format(signif(ff$dev[i]/nn, digits)),"\n",
		"    class counts: ", format(ff$yval2[i,]),"\n",
		"   probabilities: ", format(round(ff$yprob[i,], digits)),"\n")
	  }
	else if (x$method=='poisson'|x$method=='exp')
	    cat("  events=", format(ff$yval2[i]),
		",  estimated rate=" , format(signif(ff$yval[i], digits)),
		" , deviance/n=" , format(signif(ff$dev[i]/nn, digits)),"\n",
		sep = "")
	if (ff$complexity[i] > cp && !is.leaf[i] ){
	    sons <- 2*id[i] + c(0,1)
	    sons.n <- ff$n[match(sons, id)]
	    cat("  left son=", sons[1], " (", sons.n[1], " obs)",
		" right son=", sons[2], " (", sons.n[2], " obs)", sep='')
	    j <- nn - (sons.n[1] + sons.n[2])
	    if (j>0) cat(", ", j, " observations remain\n", sep='')
	    else     cat("\n")
	    cat("  Primary splits:\n")
	    j <- seq(index[i], length=1+ff$ncompete[i])
	    cat(paste("      ", format(sname[j]), " ", cuts[j],
		      " improve=", format(signif(x$splits[j,3], digits)),
		      ", (", nn - x$splits[j,1], " missing)", sep=''),
		      sep="\n")
	    if (ff$nsurrogate[i] >0) {
		cat("  Surrogate splits:\n")
		j <- seq(1 +index[i] + ff$ncompete[i], length=ff$nsurrogate[i])
		agree <- x$splits[j,3]
# I had to remove the "adjusted": to be correct the temp variable must be
#   based on "node->lastsurrogate", which is not retained in the S object
		temp  <- max(sons.n)/ nn
		cat(paste("      ", format(sname[j]), " ", cuts[j],
		      " agree=", format(signif(agree, digits)),
#                     ", adj=" , format(signif((agree-temp)/(1-temp), digits)),
		      ", (", x$splits[j,1], " split)", sep=''),
		      sep="\n")
		}
	    }
	}
    cat("\n")
    invisible(x)
    }
