\name{summary.rpart}
\alias{summary.rpart}
\title{
  Summarize a Fitted Rpart Object
}
\description{
  Returns a detailed listing of a fitted \code{rpart} object.
}
\usage{
summary.rpart(object, cp=0, digits, file, \dots)
}
\arguments{
\item{object}{
  fitted model object of class \code{rpart}.  This is assumed to be the result
  of some function that produces an object with the same named
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
  \code{"rpart"}.  It can be invoked by calling \code{summary}
  for an object of the appropriate class, or directly by calling
  \code{summary.rpart} regardless of the class of the object.
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
1 0.59535      0    1.0000 1.0197 0.17680
2 0.13453      1    0.4047 0.6141 0.10905
3 0.01283      2    0.2701 0.4711 0.08668
4 0.01000      3    0.2573 0.4710 0.08638


Node number 1: 60 observations,    complexity param=0.5953
  mean=24.58, MSE=22.58
  left son=2 (45 obs) right son=3 (15 obs)
  Primary splits:
      Weight < 2568 to the right, improve=0.5953, (0 missing)


Node number 2: 45 observations,    complexity param=0.1345
  mean=22.47, MSE=8.027
  left son=4 (22 obs) right son=5 (23 obs)
  Primary splits:
      Weight < 3088 to the right, improve=0.5045, (0 missing)


Node number 3: 15 observations
  mean=30.93, MSE=12.46


Node number 4: 22 observations
  mean=20.41, MSE=2.787


Node number 5: 23 observations,    complexity param=0.01283
  mean=24.43, MSE=5.115
  left son=10 (15 obs) right son=11 (8 obs)
  Primary splits:
      Weight < 2748 to the right, improve=0.1477, (0 missing)


Node number 10: 15 observations
  mean=23.8, MSE=4.027


Node number 11: 8 observations
  mean=25.62, MSE=4.984
}}
\keyword{rpart}
