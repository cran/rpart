#SCCS  @(#)print.rpart.s	1.11 09/03/97
print.rpart <- function(x, pretty=0, spaces=2, cp, 
               digits=.Options$digits-3, ...) {
    if(!inherits(x, "rpart")) stop("Not legitimate rpart object")

    #This is stolen, unabashedly, from print.tree
    if (x$method=='class')
         cat("node), split, n, loss, yval, (yprob)\n")
    else cat("node), split, n, deviance, yval\n")
    cat("      * denotes terminal node\n\n")
  
    if (!missing(cp)) x <- prune.rpart(x, cp=cp)
    frame <- x$frame
    ylevel <- attr(x,'ylevels')
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    indent <- paste(rep(" ", spaces * 32), collapse = "")   
    #32 is the maximal depth
    if(length(node) > 1) {
        indent <- substring(indent, 1, spaces * seq(depth))
        indent <- paste(c("", indent[depth]), format(node), ")", sep = "")
        }
    else indent <- paste(format(node), ")", sep = "")
    if (x$method=='class') {
        if(!is.null(ylevel)) 
           yval <- paste(as.character(ylevel[frame$yval]),
                                  " (", sep = "")
        else
           yval <- paste(as.character(frame$yval),
                                     " (", sep = "")
        yprob <- format(frame$yprob,digits=digits)
        for(i in 1:ncol(yprob))
            yval <- paste(yval, yprob[, i])
        yval <- paste(yval, ")")
        }
    else yval <- format(signif(frame$yval, digits = digits))
    term <- rep(" ", length(depth))
    term[frame$var == "<leaf>"] <- "*"
    z <- labels(x, pretty = pretty)
    n <- frame$n
    z <- paste(indent, z, n, format(signif(frame$dev, digits = digits)), 
            yval, term)
    cat(z, sep = "\n")
    return(invisible(x))
    #end of the theft
    }
