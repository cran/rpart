# SCCS @(#)plotcp.s	1.1 02/08/98
# Contributed by B.D. Ripley 97/07/17
#
plotcp <- function(x, minline=TRUE, lty=3, col=1,
		   upper=c("size", "splits", "none"), ...)
{
  if(!inherits(x, "rpart")) stop("Not legitimate rpart object")
  upper <- match.arg(upper)
  p.rpart <- x$cptable
  if(ncol(p.rpart) < 5)
    stop("cptable does not contain cross-validation results")
  xstd <- p.rpart[, 5]
  xerror <- p.rpart[, 4]
  nsplit <- p.rpart[, 2]
  ns <- seq(along=nsplit)
  cp0 <- p.rpart[ ,1]
  cp <- sqrt(cp0 * c(Inf, cp0[-length(cp0)]))
  ylim <- c(min(xerror - xstd) - 0.1, max(xerror + xstd) + 0.1)
  plot(ns, xerror, axes = FALSE, xlab = "cp", ylab =
       "X-val Relative Error", ylim = ylim, type = "o", ...)
  box()
  axis(2, ...)
  segments(ns, xerror - xstd, ns, xerror + xstd)
  axis(1, at = ns, lab = as.character(signif(cp, 2)), ...)
  switch(upper,
	 size = {
           axis(3, at = ns, lab = as.character(nsplit+1), ...)
           mtext("size of tree", side=3, line=3)
	 },
	 splits = {
           axis(3, at = ns, lab = as.character(nsplit), ...)
           mtext("number of splits", side=3, line=3)
	 },)
  minpos <- min(seq(along=xerror)[xerror==min(xerror)])
  if(minline) abline(h=(xerror+xstd)[minpos], lty=lty, col=col)
  invisible()
}