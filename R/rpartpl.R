rpartpl <- function(x, y, node, erase = FALSE, ...)
{

    if (!is.null(x$x) && !is.null(x$y)) {
        node <- y
        y <- x$y
        x <- x$x
    }
    parent <- match((node %/% 2L), node)
    sibling <- match(ifelse(node %% 2L, node - 1L, node + 1L), node)
    xx <- rbind(x, x, x[sibling], x[sibling], NA)
    yy <- rbind(y, y[parent], y[parent], y[sibling], NA)
    if (any(erase)) {
        ## erase denotes the set of nodes to be erased
        lines(c(xx[, erase]), c(yy[, erase]), col = 0L)
        return(x = x[!erase], y = y[!erase])
    }
    plot(range(x), range(y), type = "n", axes = FALSE, xlab = "", ylab = "")
    text(x[1L], y[1L], "|", ...)
    lines(c(xx[, -1L]), c(yy[, -1L]), ...)
    list(x = x, y = y)
}

