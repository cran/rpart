\name{rpart-internal}
\alias{is.Surv}
\alias{string.bounding.box}
\alias{tree.depth}
\alias{descendants}
\alias{node.match}
\title{Internal Functions Used by the `rpart' Package}
\description{
  Internal functions used in the \code{rpart} package.
}
\usage{
is.Surv(x)
string.bounding.box(s)
tree.depth(nodes)
descendants(nodes, include=TRUE)
node.match(nodes, nodelist, leaves, print.it = TRUE)
}
\note{
  These are not intended for use by end users.
  \code{is.Surv} duplicates \code{survival}, \code{tree.depth},
  \code{descendants} and \code{node.match} duplicate \code{tree}.
}
\keyword{internal}
