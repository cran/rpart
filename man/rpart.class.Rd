\name{rpart.class}
\alias{rpart.class}
\title{
  Set up for `class' Method for an Rpart Model
}
\description{
  Set up for the \code{class} method for an \code{rpart} model.
}
\usage{
rpart.class(y, offset, parms, wt)
}
\arguments{
  \item{y}{the responses.}
  \item{offset}{\code{NULL}: anything else is an error.}
  \item{parms}{a named list of parameters.}
  \item{wt}{case weights.}
}
\value{
  Returns a list consisting of
  \item{y}{\code{y} as factor codes,}
  \item{parms}{a named list of parameters.  This will have components
    \code{prior}, \code{loss} and \code{split} (and more if supplied).}
  \item{numresp}{the number of responses (one plus the number of classes).}
  \item{counts}{a vector of counts for each class.}
  \item{ylevels}{the response levels corresponding to the codes given in
    \code{y}.}
  \item{summary}{a function to be invoked by \code{\link{summary.rpart}},}
  \item{text}{a function to be invoked by \code{\link{text.rpart}}.}
}
\details{
  \code{rpart.class} is to be used only by the function \code{rpart}. It
  validates the input parameter list: if that is missing or \code{NULL}
  default values are supplied.
}
\seealso{
\code{\link{rpart}} 
}
\keyword{internal}
