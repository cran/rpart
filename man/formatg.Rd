\name{formatg}
\alias{formatg}
\title{Print out numbers in C format}
\description{
  Print out numbers in C format.
}
\usage{
formatg(x, digits=getOption("digits"),
        format=paste("%.", digits, "g", sep = ""))
}
\arguments{
  \item{x}{numeric vector}
  \item{digits}{number of significant digits}
  \item{format}{\code{printf}-style numeric format.}
}
\value{
  A character vector.
}
\seealso{\code{\link{formatC}}}

\keyword{tree}
