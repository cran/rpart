\name{snip.rpart.mouse}
\alias{snip.rpart.mouse}
\title{
  Internal Functions called from `snip.rpart' 
}
\description{
  Internal function called from \code{snip.rpart}. 
}
\usage{
snip.rpart.mouse(tree, parms=paste(".rpart.parms", dev.cur(), sep = "."))
}
\arguments{
\item{tree}{
  fitted model object of class \code{rpart}.  This is assumed to be the result
  of some function that produces an object with the same named
  components as that returned by the \code{rpart} function.
}
\item{parms}{
  Parameters carried over from the function \code{plot.rpart}
  (\code{uniform}, \code{branch}, \code{compress}, \code{nspace},
  \code{margin}, \dots{}).  See help page for \code{plot.rpart}
  for more details. 
}}
\value{
  Called by \code{snip.rpart}
}
\details{
  This is not intended to be called directly by users
}
\seealso{
  \code{\link{snip.rpart}}, \code{\link{plot.rpart}}
}
\keyword{rpart}

