# SCCS %W% %G%
#
#  Interactively snip off part of a tree
#

snip.rpart.mouse <- function(tree, 
		      parms=paste(".rpart.parms", dev.cur(), sep = ".")) {
    xy <- rpartco(tree)
    toss <- NULL
    ff <- tree$frame
    if (exists(parms)) {
        parms <- get(parms)
	branch <- parms$branch
	}
    else branch <- 1

    node <- as.numeric(row.names(tree$frame))
    draw <- rpart.branch(xy$x,xy$y, node, branch)

    lastchoice <- 0
    while (length(choose <- identify.rpart(tree)) >0 ) {
	if (ff$var[choose] == '<leaf>') {
		cat("Terminal node -- try again\n")
		next
		}

	if (choose != lastchoice) {
	    # print out some info on the click
	    cat("node number:", node[choose], " n=", ff$n[choose], "\n")
	    cat("    response=", format(ff$yval[choose]))
	    if (is.null(ff$yval2)) cat ("\n")
	    else if (is.matrix(ff$yval2)) 
		  cat(" (", format(ff$yval2[choose,]), ")\n")
	    else  cat(" (", format(ff$yval2[choose]), ")\n")
	    cat("    Error (dev) = ", format(ff$dev[choose]), "\n")
	    lastchoice <- choose
	    }
	else {
	    # second click-- erase all of the descendants
	    #   (stolen from snip.tree)
	    id  <- node[choose]
	    id2 <- node
	    while (any(id2>1)) {
		id2 <- floor(id2/2)
		temp  <- (match(id2, id, nomatch=0) >0)	
  	        id <- c(id, node[temp])
		id2[temp] <- 0
		}
	    temp <- match(id, node[ff$var != '<leaf>'], nomatch=0)
	    lines(c(draw$x[,temp]), c(draw$y[,temp]), col=0)
	    toss <- c(toss, node[choose])
	    }
	}
    toss
    }

identify.rpart <- function(x)
{
  xy <- treeco(x)
#  identify(xy$x, xy$y, n=1, plot=F)
  uin <- 1/xyinch()
  repeat{
    XY <- locator(1)
    if(is.null(XY)) return(NULL)
    else {
      ux <- xy$x
      uy <- xy$y
      xp <- XY$x
      yp <- XY$y
      d2 <- ((xp - ux) * uin[1])^2 + ((yp - uy) * uin[2])^2
      dist <- min(d2)
      indx <- if(dist > 0.25) 0 else seq(along=ux)[d2 == dist][1]
      if(!indx) {
        cat("No node close to point, try again\n")
      } else return(indx)
    }
  }
}
