#SCCS @(#)prune.rpart.s	1.8 02/27/98
prune.rpart <- function(tree, cp, ...)
{
    ff <- tree$frame
    id <- as.integer(row.names(ff))
    toss <- id[ff$complexity <= cp &  ff$var!='<leaf>']#not a leaf
    if (length(toss)==0) return(tree)   #all the tree is retained

    newx <- snip.rpart(tree, toss)

    ## Now cut down the CP table
    temp <- pmax(tree$cptable[,1], cp)
    keep <- match(unique(temp), temp)
    newx$cptable <- tree$cptable[keep,]
    newx$cptable[max(keep),1] <- cp

    newx
}
