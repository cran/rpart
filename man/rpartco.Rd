\name{rpartco}
\alias{rpartco}
\title{
Calculate Plotting Coordinates for Nodes of an Rpart Object
}
\usage{
rpartco(tree, parms=paste(".rpart.parms", dev.cur(), sep = "."))
}
\arguments{
\item{tree}{
fitted model object of class \code{rpart}.  This is assumed to be the result
of some function that produces an object with the same named
components as that returned by the \code{rpart} function.
}
\item{parms}{
Parameters carried over from the function \code{plot.rpart} (\code{uniform}, \code{branch},
\code{compress}, \code{nspace}, \code{margin}, \dots{}).  See help page for \code{plot.rpart} for more
details. 
}}
\value{
Called by \code{plot.rpart}.
}
\details{
This is not intended to be called directly by users.
}
\seealso{
\code{\link{plot.rpart}}.
}
\keyword{rpart,}
% Converted by Sd2Rd version 0.3-1.
