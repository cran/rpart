\name{rpart-internal}
\alias{is.Surv}
\alias{string.bounding.box}
\alias{tree.depth}
\title{Internal Functions Used by the `rpart' Package}
\description{
  Internal functions used in the \code{rpart} package.
}
\usage{
is.Surv(x)
string.bounding.box(s)
tree.depth(nodes)
}
\note{
  These are not intended for use by end users.
  \code{is.Surv} duplicated \code{survival}, \code{tree.depth}
  duplicates \code{tree}.
}
\keyword{internal}
