#
# Simplest weight test: treble the weights
#
#  By using the unshrunken estimates the weights will nearly cancel
#   out:  frame$wt, frame$dev, frame$yval2, and improvement will all
#   be threefold larger, other things will be the same.
# The improvement is the splits matrix, column 3, rows with n>0.  Other
#   rows are surrogate splits.

tempc <- rpart.control(maxsurrogate=0, cp=0, xval=0)
fit1 <- rpart(Surv(pgtime, pgstat) ~ age + eet + g2+grade+gleason +ploidy,
                stagec, control=tempc,
                method='poisson', parms=list(shrink=0))
wts <- rep(3, nrow(stagec))
fit1b <- rpart(Surv(pgtime, pgstat) ~ age + eet + g2+grade+gleason +ploidy,
                stagec, control= tempc, parms=list(shrink=0),
                method='poisson', weights=wts)
fit1b$frame$wt   <- fit1b$frame$wt/3
fit1b$frame$dev  <- fit1b$frame$dev/3
fit1b$frame$yval2<- fit1b$frame$yval2/3
fit1b$splits[,3] <- fit1b$splits[,3]/3
all.equal(fit1[-3], fit1b[-3])   #all but the "call"

#
# Compare a pair of multiply weighted fits
#  In this one, the lengths of where and y won't match
# I have to set minsplit to the smallest possible, because otherwise
#  the replicated data set will sometimes have enough "n" to split, but
#  the weighted one won't.  Use of CP keeps the degenerate splits
#  (n=2, several covariates with exactly the same improvement) at bay.
# For larger trees, the weighted split will sometimes have fewer
#  surrogates, because of the "at least two obs" rule.
#
.Random.seed <- c(53, 33, 62, 6, 45, 2, 43, 33, 16, 10, 39, 3)
wts <- sample(1:5, nrow(stagec), replace=T)
temp <- rep(1:nrow(stagec), wts)             #row replicates
xgrp <- rep(1:10, length=146)[order(runif(146))]
xgrp2<- rep(xgrp, wts)
#  Direct: replicate rows in the data set, and use unweighted
fit2 <- rpart(Surv(pgtime, pgstat) ~ age + eet + g2+grade+gleason +ploidy,
	      control=rpart.control(minsplit=2, xval=xgrp2, cp=.025),
	      data=stagec[temp,], method='poisson')
#  Weighted
fit2b<- rpart(Surv(pgtime, pgstat) ~ age + eet + g2+grade+gleason +ploidy,
	      control=rpart.control(minsplit=2, xval=xgrp, cp=.025),
	      data=stagec, method='poisson', weight=wts)

all.equal(fit2$frame[-2],  fit2b$frame[-2])  # the "n" component won't match
all.equal(fit2$cptable,    fit2b$cptable)
all.equal(fit2$splits[,-1],fit2b$splits[,-1]) #fails
all.equal(fit2$splits[-24,-1],fit2b$splits[-24,-1]) #ok
all.equal(fit2$csplit,    fit2b$csplit)
# Line 24 is a surrogate split in a group whose 2 smallest ages are
#  47 and 48.  The weighted fit won't split there because it wants to
#  send at least 2 obs to the left; the replicate fit thinks that there
#  are several 47's.





