\name{rpart.branch}
\alias{rpart.branch}
\title{
  Draws horseshoes - left child, up over, down to right child
}
\description{
  Draws a branch of a tree.
}
\usage{
rpart.branch(x, y, node, branch)
}
\arguments{
\item{x}{
x coordinates of nodes 
}
\item{y}{
y coordinates of nodes 
}
\item{node}{
description of node number 
}
\item{branch}{
controls shape of horseshoe (1 is default, 0 has no shoulder) 
}
}
\value{
  Internal function, called by \code{plot.rpart}.  Not usually called by user.
}
\keyword{tree}
% docclass is function
% Converted by Sd2Rd version 1.20.
