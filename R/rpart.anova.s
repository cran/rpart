#SCCS @(#)rpart.anova.s	1.1 09/19/95
rpart.anova <- function(y, offset, parms) {
    if (!is.null(offset)) y <- y-offset
    list(y=y, parms=0, numresp=1)
    }
