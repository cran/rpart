na.rpart <- function(x)
{
    y <- x[[1]]
    if (is.matrix(y)){
        keep <- apply(y, 1, function(x) !all(is.na(x)))
    } else keep <- !is.na(y)
    xmiss <- is.na(x[,seq(along=names(x))[-1], drop=FALSE])
    keep <- keep & ((xmiss %*% rep(1, ncol(xmiss))) < ncol(xmiss))
    if (all(keep)) x
    else {
        temp <- seq(keep)[!keep]
        names(temp) <- row.names(x)[!keep]
#the methods for this group are all the same as for na.omit
        attr(temp, "class") <- c('na.rpart', 'omit')
        structure(x[keep,], na.action=temp)
    }
}
