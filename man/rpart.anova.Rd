\name{rpart.anova}
\alias{rpart.anova}
\title{
  Set up for `anova' Method for an Rpart Model
}
\description{
  Set up for \code{anova} method for an \code{rpart} model.
\usage{
rpart.anova(y, offset, parms, wt)
}
\value{
  Returns a list consisting of \code{y} (substracting \code{offset}
  if necessary), \code{parms}, and \code{numresp=1}).
  \code{rpart.anova} is used only by the function \code{rpart}.
}
\keyword{tree}
% docclass is function
% Converted by Sd2Rd version 1.20.
