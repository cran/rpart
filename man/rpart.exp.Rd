\name{rpart.exp}
\alias{rpart.exp}
\title{
  Set up for `exp' Method for an Rpart Model
}
\description{
  Set up for the \code{exp} method for an \code{rpart} model.
}
\usage{
rpart.exp(y, offset, parms, wt)
}
\arguments{
  \item{y}{the responses, a vector or a two-column matrix.}
  \item{offset}{an offset, or \code{NULL}.}
  \item{parms}{a list of parameters.}
  \item{wt}{case weights.}
}
\value{
  Returns a list consisting of
  \item{y}{a matrix giving the observation times and number of events,}
  \item{parms}{a named list of parameters, with components \code{shrink}
    and \code{method} (1 for \code{"deviance"}, 2 for \code{"sqrt"}),}
  \item{numresp}{the number of responses, here \code{2},}
  \item{numy}{the number of columns in \code{y}, here \code{2},}
  \item{summary}{a function to be invoked by \code{\link{summary.rpart}},}
  \item{text}{a function to be invoked by \code{\link{text.rpart}}.}
}
\details{
  \code{rpart.exp} is to be used only by the function \code{rpart}.  It
  transforms the data to look like a unit-rate Poisson process so that
  the \code{rpart.exp} method can be used.
}
\keyword{internal}
