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
\method{predict}{rpart}(object, newdata=list(),
       type=c("vector", "matrix", "class", "probs"))
}
\arguments{
  \item{object}{
    fitted model object of class \code{rpart}.  This is assumed to be
    the result of some function that produces an object with the same named 
    components as that returned by the \code{rpart} function. 
  }
  \item{newdata}{
    data frame containing the values at which predictions are required. 
    The predictors referred to in the right side of
    \code{formula(object)} must be present by name in \code{newdata}.
    If missing, the fitted values are returned. 
  }
  \item{type}{
    character string denoting whether the predictions are returned as a 
    vector (default except for classifcation trees) or as an
    \code{rpart} object, or as a matrix or (for
    classification trees) a factor of classes or a matrix of class
    probabilities (default for classification trees).
  }
}
\value{
  A new object is obtained by
  dropping \code{newdata} down the object.  For factor predictors, if an
  observation contains a level not used to grow the tree, it is left at
  the deepest possible node and \code{frame$yval} at the node is the
  prediction.

  If \code{type="vector"}:\cr
  vector of predicted responses.
  
%    If \code{type="tree"}:\cr
%    an object of class \code{rpart} is returned with new values for
%    \code{frame$n} and \code{frame$yval} (and \code{frame$yval2} if it
%    exists).

  If \code{type="matrix"}:\cr
  a matrix of the full responses (\code{frame$yval2} if this exists,
  otherwise \code{frame$yval}).   For a Poisson tree this is
  (event rate, \#events), and for a classification tree it is
  (predicted class, class counts, class probabilities).

  If \code{type="class"}:\cr
  (for a classification tree) a factor of classifications based on the
  responses.

  If \code{type="probs"}:\cr
  (for a classification tree) a matrix of class probabilities.
}
\details{
  This function is a method for the generic function predict for class
  \code{rpart}. It can be invoked by calling \code{predict} for an object
  of the appropriate class, or directly by calling \code{predict.rpart}
  regardless of the class of the object.
}
\seealso{
\code{\link{predict}},  \code{\link{rpart.object}}   
}
\examples{
data(car.test.frame)
z.auto <- rpart(Mileage ~ Weight, car.test.frame)
predict(z.auto)

data(kyphosis)
fit <- rpart(Kyphosis ~ Age + Number + Start, data=kyphosis)
predict(fit, type="vector") # level numbers
predict(fit, type="class")  # factor
predict(fit, type="probs")  # class probabilities (default)
predict(fit, type="matrix") # level number, class frequencies, probabilities

data(iris)
sub <- c(sample(1:50, 25), sample(51:100, 25), sample(101:150, 25))
fit <- rpart(Species ~ ., data=iris, subset=sub)
fit
table(predict(fit, iris[-sub,], type="class"), iris[-sub, "Species"])
}
\keyword{tree}
% docclass is function
% Converted by Sd2Rd version 1.20.
