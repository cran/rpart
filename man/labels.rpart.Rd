\name{labels.rpart}
\alias{labels.rpart}
\title{
  Create Split Labels For an Rpart Object
}
\usage{
labels.rpart(object, pretty = TRUE, collapse = TRUE)
}
\description{
  Create split labels for an \code{rpart} object.
}
\arguments{
\item{object}{
  fitted model object of class \code{rpart}.  This is assumed to be the result
  of some function that produces an object with the same named
  components as that returned by the \code{rpart} function.
}
\item{pretty}{
  an integer denoting the extent to which factor levels will be
  abbreviated in the character strings defining the splits;
  (0) signifies no abbreviation of levels.  A \code{NULL} signifies using
  elements of letters to represent the different factor levels.  The
  default (TRUE) indicates the maximum possible abbreviation.
}
\item{collapse}{
  logical.  If \code{TRUE} (default), the returned value is a vector of labels
  for each node, otherwise a matrix of left and right split labels is
  returned. 
}}
\value{
  Vector of split labels (\code{collapse=RUET}) or matrix of left and
  right splits (\code{collapse=FALSE}) for the supplied \code{rpart}
  object.  This function is called by printing methods for \code{rpart}
  and is not intended to be called directly by the users.
}
\seealso{
\code{\link{abbreviate}}
}
\keyword{rpart}




