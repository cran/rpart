

\name{rpart.control}
\alias{rpart.control}
\title{
Control for Rpart Models.
}
\description{
Various parameters that control aspects of the \code{rpart} fit.
}
\usage{
rpart.control(minsplit=20, minbucket=round(minsplit/3), cp=0.01,
              maxcompete=4, maxsurrogate=5, usesurrogate=2, xval=10)
}
\arguments{
\item{minsplit}{
the minimum number of observations that must exist in a node, in order for
a split to be attempted.
}
\item{minbucket}{
the minimum number of observations in any terminal \code{<leaf>} node.
If only one of \code{minbucket} or \code{minsplit} is specified, the code either sets
\code{minsplit} to \code{minbucket*3} or \code{minbucket} to \code{minsplit/3}, as appropriate.
}
\item{cp}{
complexity parameter.  Any split that does not decrease the overall lack of
fit by a factor of \code{cp} is not attempted.  For instance, with \code{anova} splitting,
this means that the overall Rsquare must increase by \code{cp} at each step.
The main role of this parameter is to save computing time by pruning off
splits that are obviously not worthwhile.
Essentially,the user informs the program that any split which does not improve
the fit by \code{cp} will likely be pruned off by cross-validation, and that hence
the program need not pursue it.
}
\item{maxcompete}{
the number of competitor splits retained in the output.  It is useful to
know not just which split was chosen, but which variable came in second,
third, etc.
}
\item{maxsurrogate}{
the number of surrogate splits retained in the output.  If this is set to
zero the compute time will be shortened, since approximately half of the
computational time (other than setup) is used in the search for surrogate
splits.
}
\item{usesurrogate}{
how to use surrogates in the splitting process.  0= display only; an
observation with a missing value for the primary split rule is not sent
further down the tree.  1= use surrogates, in order, to split subjects
missing the primary variable; if all surrogates are missing the observation
is not split.  2= if all surrogates are missing, then send the observation
in the majority direction.  A value of 0 corresponds to the action of \code{tree},
and 2 to the recommendations of Breiman, et.al.
}
\item{xval}{
number of cross-validations
}}
\value{
a list containing the options
}
\seealso{
\code{\link{rpart}}
}
\keyword{rpart}
% Converted by Sd2Rd version 0.3-1.
