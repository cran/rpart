\name{predict.rpart}
\alias{predict.rpart}
\alias{pred.rpart}
\title{
Predictions from a Fitted Rpart Object
}
\description{
Returns a vector of predicted responses from a fitted \code{rpart} object.
}
\usage{
predict.rpart(object, newdata=list(), type=c("vector", "tree"))
}
\arguments{
\item{object}{
fitted model object of class \code{rpart}.  This is assumed to be the result
of some function that produces an object with the same named
components as that returned by the \code{rpart} function.
}
\item{newdata}{
data frame containing the values at which predictions are required.
The predictors referred to in the right side of \code{formula(object)} must
be present by name in \code{newdata}.  If missing, the fitted values are returned.
}
\item{type}{
character string denoting whether the predictions are returned as a
vector (default) or as an \code{rpart} object.
}}
\value{
if \code{type="vector"}:
vector of predicted responses or, if the response is a factor, matrix
of predicted class probabilities.  This new object is obtained by
dropping newdata down object.  For factor predictors, if an
observation contains a level not used to grow the tree, it is left at
the deepest possible node and \code{frame$yval} at the node is the
prediction.
if \code{type="tree"}
an object of class \code{rpart} is returned with new values for \code{frame$n} and
\code{frame$yval} (and \code{frame$yprob} if it exists).
}
\details{
This function is a method for the generic function predict for class
\code{rpart}. It can be invoked by calling predict for an object of the
appropriate class, or directly by calling \code{predict.rpart} regardless of
the class of the object.
}
\seealso{
\code{\link{predict}}, \code{\link{rpart.object}}
}
\examples{
data(car.test.frame)
z.auto <- rpart(Mileage ~ Weight, car.test.frame)
predict(z.auto)
}
\keyword{rpart,}
% Converted by Sd2Rd version 0.3-1.






