\name{residuals.rpart}
\alias{residuals.rpart}
\title{
  Residuals From a Fitted Rpart Object 
}
\usage{
residuals.rpart(object, type)
}
\description{
  Method for \code{residuals} for an \code{rpart} object.
}
\arguments{
\item{object}{
fitted model object of class \code{rpart}.  This is assumed to be the result 
of some function that produces an object with the same named 
components as that returned by the \code{rpart} function. 
}
\item{type}{
  Indicates method of residual for a classification tree (no distinction
  is made for methods \code{"anova"}, \code{"poisson"}, or
  \code{"exp"}).  Allows the
  following types of residuals: \code{"usual"} (misclassification errors),
  \code{"pearson"} (sum of the observed minus fitted squared over fitted), and
  \code{"deviance"} (minus twice logarithm of fitted).
}}
\value{
vector of residuals of type \code{type} from a fitted \code{rpart} object. 
}
\examples{
data(solder)
fit <- rpart(skips ~ Opening + Solder + Mask + PadType + Panel,
       data=solder, method='anova')
summary(residuals(fit))
plot(predict(fit),residuals(fit))
}
\keyword{tree}
% docclass is function
% Converted by Sd2Rd version 1.20.
