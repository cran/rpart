\name{post.rpart}
\alias{post.rpart}
\title{
PostScript Presentation Plot of an Rpart Object
}
\description{
Generates a PostScript presentation plot of an \code{rpart} object.
}
\usage{
post.rpart(tree, title, 
           filename=paste(deparse(substitute(tree)), ".ps", sep = ""), 
           digits= .Options$digits -3, pretty=T,
	   use.n=T, horizontal=T, \dots)
}
\arguments{
\item{tree}{
fitted model object of class \code{rpart}.  This is assumed to be the result
of some function that produces an object with the same named
components as that returned by the \code{rpart} function.
}
\item{title}{
a title which appears at the top of the plot.  By default, the
name of the \code{rpart} endpoint is printed out.
}
\item{filename}{
ASCII file to contain the output.  By default, the name of the file is
the name of the object given by \code{rpart} (with the suffix \code{.ps} added). If
\code{filename = ""}, the plot appears on the current graphical device.
}
\item{digits}{
number of significant digits to include in numerical data.
}
\item{pretty}{
an integer denoting the extent to which factor levels will be
abbreviated in the character strings defining the splits;
(0) signifies no abbreviation of levels.  A \code{NULL} signifies using
elements of letters to represent the different factor levels.  The
default (\code{T}) indicates the maximum possible abbreviation.
}
\item{use.n}{
Logical.  If \code{TRUE} (default), adds to label 
(#events level1/ #events level2/etc. for method \code{class}, \code{n} for 
method \code{anova}, and # events/n for methods \code{poisson} and \code{exp}).
}
\item{horizontal}{
Logical.  If \code{TRUE} (default), plot is horizontal. If \code{FALSE}, plot
appears as landscape.
}
\item{...}{
other arguments to the \code{postscript} function.
}}
\section{Side Effects}{
a plot of \code{rpart} is created using the \code{postscript} driver.
}
\details{
The plot created uses the functions \code{plot.rpart} and \code{text.rpart} (with
the \code{fancy} option).  The settings were chosen because they looked
good to us, but other options may be better, depending on the \code{rpart}
object.  
Users are encouraged to write their own function containing favorite
options.
}
\seealso{
\code{\link{plot.rpart}}, \code{\link{rpart}}, \code{\link{text.rpart}}, \code{\link{abbreviate}}
}
\examples{
data(car.test.frame)
z.auto <- rpart(Mileage ~ Weight, car.test.frame)
post(z.auto, file="")   # display tree on active device
   # now construct postscript version on file "pretty"
   # with no title
post(z.auto, file="pretty", title=" ")
z.hp <- rpart(Mileage ~ Weight + HP, car.test.frame)
post(z.hp)
}
\keyword{rpart}
% Converted by Sd2Rd version 0.3-1.
