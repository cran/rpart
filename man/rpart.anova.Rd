\name{rpart.anova}
\alias{rpart.anova}
\title{
  Set up for `anova' Method for an Rpart Model
}
\description{
  Set up for \code{anova} method for an \code{rpart} model.
}
\usage{
rpart.anova(y, offset, parms, wt)
}
\arguments{
  \item{y}{the responses.}
  \item{offset}{an offset, or \code{NULL}.}
  \item{parms}{a list of parameters, usually empty.}
  \item{wt}{case weights.}
}
\value{
  Returns a list consisting of
  \item{y}{(subtracting \code{offset} if necessary),}
  \item{parms}{as input,}
  \item{numresp}{the number of responses, here \code{1},}
  \item{summary}{a function to be invoked by \code{\link{summary.rpart}},}
  \item{text}{a function to be invoked by \code{\link{text.rpart}}.}
}
\details{
  \code{rpart.anova} is to be used only by the function \code{rpart}.
}
\keyword{internal}
