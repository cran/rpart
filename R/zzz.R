.First.lib <- function(lib, pkg) library.dynam("rpart", pkg, lib)

is.Surv <- function(x) inherits(x, "Surv")

tree.depth <- function (nodes)
{
    depth <- floor(log(nodes, base = 2) + 1e-7)
    as.vector(depth - min(depth))
}

string.bounding.box <- function(s)
{
    s2 <- strsplit(s, "\n")
    rows <- sapply(s2, length)
    columns <- sapply(s2, function(x) max(nchar(x)))
    list(columns=columns, rows=rows)
}
