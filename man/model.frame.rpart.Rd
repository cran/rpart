\name{model.frame.rpart}
\alias{model.frame.rpart}
\title{
Create a Model Frame for fitting an Rpart Object
}
\usage{
model.frame.rpart(formula, \dots)
}
\description{
  \code{model.frame} method for \code{rpart} objects.
}
\arguments{
\item{formula}{
  the formula or other object defining what terms should be included in
  the model frame.
}}
\value{
  An object of class \code{model.frame} suitable for fitting an
  \code{rpart}-based model.  Called by \code{rpart} and other
  functions for manipulating \code{rpart} objects. 
}
\details{
  This function is a method for the generic function \code{model.frame} for
  class \code{rpart}.  It can be invoked by calling \code{model.frame}
  for an object of the appropriate class, or directly by calling
  \code{model.frame.rpart} regardless of the class of the object.
}
\seealso{
  \code{\link{model.frame}}
}
\keyword{rpart}





