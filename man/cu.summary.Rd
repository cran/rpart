\name{cu.summary}
\alias{cu.summary}
\non_function{}
\title{Automobile Data from `Consumer Reports' 1990}
\description{
  The \code{cu.summary} data frame has 117 rows and 5 columns,
  giving data on makes of cars taken from the April, 1990 issue of
  \emph{Consumer Reports}.
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
\code{Brazil} 
\code{England} 
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
an ordered factor with levels
\code{Much worse} < \code{worse} < \code{average} < \code{better} < \code{Much better}
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
}
\source{
  \emph{Consumer Reports}, April, 1990, pp. 235--288 quoted in

  John M. Chambers and Trevor J. Hastie eds. (1992)
  \emph{Statistical Models in S}, Wadsworth and Brooks/Cole,
  Pacific Grove, CA 1992, pp. 46--47.
}
\seealso{\code{\link{car.test.frame}}}
\examples{
data(cu.summary)
fit <- rpart(Price ~ Mileage + Type + Country, cu.summary)
plot(fit, compress=T)
text(fit, use.n=T)
}
\keyword(datasets)
