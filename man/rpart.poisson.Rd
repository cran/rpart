\name{rpart.poisson}
\alias{rpart.poisson}
\title{
  Set up for `poisson' Method for an Rpart Model 
}
\description{
  Set up for \code{poisson} method for an \code{rpart} model. 
}
\usage{
rpart.poisson(y, offset, parms, wt)
}
\arguments{
  \item{y}{the responses, a vector or a two-column matrix.  If a matrix
    the first column is the observation times and the second is the counts.}
  \item{offset}{an offset, or \code{NULL}.  The offset is interpreted as
  a log observation time.}
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
  \code{rpart.poisson} is used only by the function \code{rpart}.
}
\keyword{tree}
% docclass is function
