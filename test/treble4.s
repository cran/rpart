#
# Treble test for class trees with 2 outcomes
#
fit1 <- rpart(Kyphosis ~ Age + Number + Start, data=kyphosis, 
               control=rpart.control(maxsurrogate=0, cp=0, xval=0),
               parms=list(prior=c(.7,.3), 
                          loss=matrix(c(0,1,2,0),nrow=2,ncol=2)))
wts <- rep(3, nrow(kyphosis))
fit1b <- rpart(Kyphosis ~ Age + Number + Start, data=kyphosis, 
                control=rpart.control(maxsurrogate=0, cp=0, xval=0),
	       weights=wts,
               parms=list(prior=c(.7,.3), 
                          loss=matrix(c(0,1,2,0),nrow=2,ncol=2)))
fit1b$frame$wt   <- fit1b$frame$wt/3
fit1b$frame$dev  <- fit1b$frame$dev/3
fit1b$frame$yval2<- fit1b$frame$yval2/3
fit1b$splits[,3] <- fit1b$splits[,3]/3
all.equal(fit1[-3], fit1b[-3])   #all but the "call"


# Now for a set of non-equal weights
nn <- nrow(kyphosis)
wts <- sample(1:5, nn, replace=T)
temp <- rep(1:nn, wts)             #row replicates
xgrp <- rep(1:10, length=nn)[order(runif(nn))]
xgrp2<- rep(xgrp, wts)
tempc <- rpart.control(minsplit=2, xval=xgrp2, maxsurrogate=0)
#  Direct: replicate rows in the data set, and use unweighted
fit2 <- rpart(Kyphosis ~ Age + Number + Start, data=kyphosis[temp,], 
               control=tempc, 
               parms=list(prior=c(.7,.3), 
                          loss=matrix(c(0,1,2,0),nrow=2,ncol=2)))
#  Weighted
tempc <- rpart.control(minsplit=2, xval=xgrp, maxsurrogate=0)
fit2b <- rpart(Kyphosis ~ Age + Number + Start, data=kyphosis, 
               control=tempc, weights=wts,
               parms=list(prior=c(.7,.3), 
                          loss=matrix(c(0,1,2,0),nrow=2,ncol=2)))

all.equal(fit2$frame[-2],  fit2b$frame[-2])  # the "n" component won't match
all.equal(fit2$cptable,    fit2b$cptable)
all.equal(fit2$splits[,-1],fit2b$splits[,-1]) #fails
all.equal(fit2$csplit,    fit2b$csplit)
