# SCCS @(#)labels.rpart.s	1.2 01/21/97
# Differs from labels.tree only in the processing of xlevels
#
labels.rpart <- function(object, pretty = TRUE, collapse = TRUE)
{
    if(!inherits(object, "rpart"))
	stop("Not legitimate rpart object")
    frame <- object$frame
    xlevels <- attr(object, "xlevels")
    var <- frame$var
    splits <- frame$splits
    if(!is.null(pretty)) {
      if(pretty) xlevels <- lapply(xlevels, abbreviate, pretty)
      var2 <- as.character(as.numeric(var)-1)
      fix <- grep("^:", splits[, 1])
      for(i in fix)
        for(j in 1:2) {
          # split :cde into c(3,4,5) and look up levels.
          sh <- splits[i, j]
          nc <- nchar(sh)
          sh <- substring(sh, 2:nc, 2:nc)
          xl <- xlevels[[ var2[i] ]][ match(sh, letters) ]
          splits[i, j] <- paste(":", paste(as.vector(xl), collapse=","), sep="")
        }
    }
    var <- as.character(var)
    if(!collapse)
	return(array(paste(var, splits, sep = ""), dim(splits)))
    node <- as.numeric(row.names(frame))
    parent <- match((node %/% 2), node)
    odd <- as.logical(node %% 2)
    node[odd] <- paste(var[parent[odd]], splits[parent[odd], 2], sep = "")
    node[!odd] <- paste(var[parent[!odd]], splits[parent[!odd], 1], sep =
		"")
    node[1] <- "root"
    node
}
