.First.lib <- function(lib, pkg) library.dynam("rpart", pkg, lib)

if(version$major == "0" && version$minor < "0.63")
  labels <- function(object, ...) UseMethod("labels")

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
