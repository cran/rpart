\name{print.rpart}
\alias{print.rpart}
\title{
Print an Rpart Object
}
\description{
This function prints an \code{rpart} object.  It is a method for the generic
function \code{print} of class \code{rpart}.
}
\usage{
print.rpart(x, pretty=0, spaces=2, cp=0, digits= .Options$digits - 3, \dots)
}
\arguments{
\item{x}{
fitted model object of class \code{rpart}.  This is assumed to be the result
of some function that produces an object with the same named
components as that returned by the \code{rpart} function.
}
\item{pretty}{
a logical or \code{NULL} denoting the extent to which factor levels in split
labels will be abbreviated.  A value of \code{FALSE} signifies no
abbreviation.  A \code{NULL} signifies using elements of letters to represent
the different factor levels.
}
\item{spaces}{
the number of spaces to indent nodes of increasing depth.
}
\item{digits}{
the number of digits of numbers to print.
}
\item{cp}{
prune all nodes with a complexity less than \code{cp} from the printout.
}}
\section{Side Effects}{
A semi-graphical layout of the contents of \code{x$frame} is
printed.  Indentation is used to convey the tree topology.
Information for each node includes the node number, split, size,
deviance, and fitted value.  For the \code{"class"} method, \code{yprob} is also
printed. 
}
\details{
This function is a method for the generic function \code{print} for class
\code{"rpart"}.  It can be invoked by calling print for an object of the
appropriate class, or directly by calling \code{print.rpart} regardless of
the class of the object.
}
\seealso{
\code{\link{print}}, \code{\link{rpart.object}}, \code{\link{summary.rpart}}, \code{\link{printcp}}
}
\examples{
data(car.test.frame)
z.auto <- rpart(Mileage ~ Weight, car.test.frame)
z.auto
\dontrun{node), split, n, deviance, yval
      * denotes terminal node


1) root 60 1355.00 24.58  
  2) Weight>2567.5 45  361.20 22.47  
    4) Weight>3087.5 22   61.32 20.41 *
    5) Weight<3087.5 23  117.70 24.43 *
  3) Weight<2567.5 15  186.90 30.93 *
}
}
\keyword{rpart}
% Converted by Sd2Rd version 0.3-1.
