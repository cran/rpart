##
## Do Rpart predictions given a tree and a matrix of predictors
pred.rpart <- function(fit, x)
{
    frame <- fit$frame
    if (nrow(frame)==1L)                # root only tree
        return(structure(rep(1, nrow(x), names = rownames(x))))

    nc <- frame[, c('ncompete', 'nsurrogate')]
    frame$index <- 1 + c(0, cumsum((frame$var != "<leaf>") +
                                   nc[[1L]] + nc[[2L]]))[-(nrow(frame)+1L)]
    frame$index[frame$var == "<leaf>"] <- 0
    vnum <- match(dimnames(fit$splits)[[1L]], dimnames(x)[[2L]])
    if (any(is.na(vnum)))
        stop("Tree has variables not found in new data")
    temp <- .C(C_pred_rpart,
               as.integer(dim(x)),
               as.integer(dim(frame)[1L]),
               as.integer(dim(fit$splits)),
               as.integer(if(is.null(fit$csplit)) rep(0L, 2)
                          else dim(fit$csplit)),
               as.integer(row.names(frame)),
               as.integer(unlist(frame[, c('n', 'ncompete', 'nsurrogate', 'index')])),
               as.integer(vnum),
               as.double(fit$splits),
               as.integer(fit$csplit -2L),
               as.integer((fit$control)$usesurrogate),
               as.double(x),
               as.integer(is.na(x)),
               where = integer(dim(x)[1L]),
               NAOK =TRUE)
    temp <- temp$where
    names(temp) <- rownames(x)
    temp
}
