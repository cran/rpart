\name{car.test.frame}
\alias{car.test.frame}
\non_function{}
\title{Automobile Data from `Consumer Reports' 1990}
\description{
  The \code{car.test.frame} data frame has 60 rows and 8 columns,
  giving data on makes of cars taken from the April, 1990 issue of
  \emph{Consumer Reports}. This is part of a larger dataset, some
  columns of which are given in \code{\link{cu.summary}}.
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{Price}{
a numeric vector giving the list price in US dollars of a standard model
}
\item{Country}{
of origin, a factor with levels
\code{France} 
\code{Germany} 
\code{Japan} 
\code{Japan/USA} 
\code{Korea} 
\code{Mexico} 
\code{Sweden} 
\code{USA} 
}
\item{Reliability}{
a numeric vector coded \code{1} to \code{5}.
}
\item{Mileage}{
fuel consumption miles per US gallon, as tested. 
}
\item{Type}{
a factor with levels
\code{Compact} 
\code{Large} 
\code{Medium} 
\code{Small} 
\code{Sporty} 
\code{Van} 
}
\item{Weight}{
kerb weight in pounds.
}
\item{Disp.}{
the engine capacity (displacement) in litres.
}
\item{HP}{
the net horsepower of the vehicle.
}
}
\source{
  \emph{Consumer Reports}, April, 1990, pp. 235--288 quoted in

  John M. Chambers and Trevor J. Hastie eds. (1992)
  \emph{Statistical Models in S}, Wadsworth and Brooks/Cole,
  Pacific Grove, CA 1992, pp. 46--47.
}
\seealso{\code{\link{cu.summary}}}
\examples{
data(car.test.frame)
z.auto <- rpart(Mileage ~ Weight, car.test.frame)
summary(z.auto)
}
\keyword(datasets)
