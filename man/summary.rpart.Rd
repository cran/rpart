\name{summary.rpart}
\alias{summary.rpart}
\title{
Summarize a Fitted Rpart Object
}
\description{
Returns a detailed listing of a fitted \code{rpart} object.
}
\usage{
summary.rpart(object, cp=0, file, digits)
}
\arguments{
\item{object}{
fitted model object of class \code{rpart}.  This is assumed to be the result
of some funciton that produces an object with the same named
components as that returned by the \code{rpart} function.
}
\item{cp}{
trim nodes with a complexity of less than \code{cp} from the listing.
}
\item{file}{
write the output to a given file name. (Full listings of a tree are
often quite long).
}}
\details{
This function is a method for the generic function summary for class
\code{"rpart"}.  It can be invoked by calling \code{summary} for an object of the
appropriate class, or directly by calling \code{summary.rpart} regardless of
the class of the object.
}
\seealso{
\code{\link{summary}}, \code{\link{rpart.object}}, \code{\link{printcp}}.
}
\examples{
data(car.test.frame)
z.auto <- rpart(Mileage ~ Weight, car.test.frame)
summary(z.auto)
\dontrun{
      CP nsplit rel error xerror    xstd 
1 0.5953      0    1.0000 1.0078 0.17510
2 0.1345      1    0.4047 0.4983 0.07396
3 0.0100      2    0.2701 0.3946 0.06354

Node number 1: 60 observations,    complexity param=0.5953
  mean=24.58 , SS/n=22.58
  left son=2 (45 obs) right son=3 (15 obs)
  Primary splits:
      Weight < 2568 to the right, improve=0.5953, (0 missing)

Node number 2: 45 observations,    complexity param=0.1345
  mean=22.47 , SS/n=8.027
  left son=4 (22 obs) right son=5 (23 obs)
  Primary splits:
      Weight < 3088 to the right, improve=0.5045, (0 missing)

Node number 3: 15 observations
  mean=30.93 , SS/n=12.46

Node number 4: 22 observations
  mean=20.41 , SS/n=2.787

Node number 5: 23 observations
  mean=24.43 , SS/n=5.115
}
}
\keyword{rpart}
% Converted by Sd2Rd version 0.3-1.








