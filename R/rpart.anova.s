#SCCS @(#)rpart.anova.s	1.2 12/13/99
rpart.anova <- function(y, offset, parms, wt) {
    if (!is.null(offset)) y <- y-offset
    list(y=y, parms=0, numresp=1)
    }
