# SCCS 03/03/98 @(#)post.rpart.s	1.11
#
post.rpart <- function(tree, title., 
		       filename=paste(deparse(substitute(tree)),".ps",sep=""),
		       digits=.Options$digits - 3, pretty=T, 
		       use.n=T,  horizontal=T, ...)
{
  if(filename !=""){
    postscript(file = filename, horizontal=horizontal, ...)
    par(mar=c(2,2,4,2)+.1)
    on.exit(dev.off())
  } else {
    oldpar <- par(mar=c(2,2,4,2)+.1)
    on.exit(invisible(par(oldpar)))
  }

  plot(tree, uniform=T, branch=.2, compress=T, margin=.1)
  text(tree, all=T, use.n=use.n, fancy=T, digits=digits, pretty=pretty)
  method <- tree$method

  if(missing(title.)) {
    temp  <- attr(tree$terms,'variables')[2]	      
    title(paste("Endpoint =",temp),cex=.8)
  } else if (title. !="") title(title.,cex=.8)
}

