\name{text.rpart}
\alias{text.rpart}
\title{
Place Text on a Dendrogram
}
\description{
Labels the current plot of the tree dendrogram with text.
}
\usage{
text.rpart(x, splits=T, label="yval", FUN=text, all=F, pretty=NULL, 
digits=.Options$digits - 3, use.n=F, fancy=F, fwidth=0.8, fheight=0.8, \dots)
}
\arguments{
\item{x}{
fitted model object of class \code{rpart}.  This is assumed to be the result
of some function that produces an object with the same named
components as that returned by the \code{rpart} function.
}
\item{splits}{
logical flag.  If \code{TRUE} (default), then the splits in the tree are
labeled with the criterion for the split.
}
\item{label}{
a column name of \code{x$frame};  values of this will label the nodes.
}
\item{FUN}{
the name of a labeling function, e.g. \code{text}.  
}
\item{all}{
Logical.  If \code{TRUE}, all nodes are labeled, otherwise just terminal nodes.
}
\item{pretty}{
an integer denoting the extent to which factor levels in split labels
will be abbreviated.  A value of (0) signifies no abbreviation.  A
\code{NULL}, the default, signifies using elements of letters to represent
the different factor levels.
}
\item{digits}{
number of significant digits to include in numerical labels.
}
\item{use.n}{
Logical.  If \code{TRUE} (default), adds to label 
(#events level1/ #events level2/etc. for \code{class}, \code{n} for \code{anova}, and
# events/n for \code{poisson} and \code{exp}).
}
\item{fancy}{
Logical.  If \code{TRUE}, nodes are represented by ellipses (interior nodes)
and rectangles (leaves) and labeled by \code{yval}.  The edges connecting the
nodes are labeled by left and right splits.
}
\item{fwidth}{
Relates to option fancy and the width of the ellipses and rectangles.
If \code{fwidth < 1} then it is a scaling factor (default = .8).  If 
\code{fwidth > 1} then it representes the number of character widths 
(for current graphical device) to use.
}
\item{fheight}{
Relates to option fancy and the height of the ellipses and rectangles.
If \code{fheight <1} then it is a scaling factor (default = .8).  If
\code{fheight > 1} then it represents the number of character heights 
(for current graphical device) to use.
}
\item{...}{
Graphical parameters may also be supplied as arguments to this
function (see \code{par}).
}}
\section{Side Effects}{
the current plot of a tree dendrogram is labeled.
}
\seealso{
  \code{\link{text}}, \code{\link{plot.rpart}}, \code{\link{rpart}},
  \code{\link{post.rpart}}, \code{\link{abbreviate}} 
}
\examples{
data(freeny)
freen.tr <- rpart(y ~ ., freeny)
plot(freen.tr)
text(freen.tr, use.n=T, all=T)
}
\keyword{rpart}
% Converted by Sd2Rd version 0.3-1.
