#
# Test weights in a regression problem
#

xgrp <- rep(1:10,5)
fit4 <- rpart(income ~ population + region + illiteracy +life + murder +
                        hs.grad + frost , mystate,
                   control=rpart.control(minsplit=10, xval=xgrp))
wts <- rep(3, nrow(mystate))
fit4b <-  rpart(income ~ population + region + illiteracy +life + murder +
                        hs.grad + frost , mystate,
                   control=rpart.control(minsplit=10, xval=xgrp), weights=wts)
fit4b$frame$wt   <- fit4b$frame$wt/3
fit4b$frame$dev  <- fit4b$frame$dev/3
fit4b$cptable[,5] <- fit4b$cptable[,5] * sqrt(3)
temp <- c('frame', 'where', 'splits', 'csplit', 'cptable')
all.equal(fit4[temp], fit4b[temp])  


# Next is a very simple case, but worth keeping
dummy <- data.frame(y=1:10, x1=c(10:4, 1:3), x2=c(1,3,5,7,9,2,4,6,8,0))

xx1 <- rpart(y ~ x1 + x2, dummy, minsplit=4, xval=0)
xx2 <- rpart(y ~ x1 + x2, dummy, weights=rep(2,10), minsplit=4, xval=0)

all.equal(xx1$frame$dev, c(82.5, 10, 2, .5, 10, .5, 2))
all.equal(xx2$frame$dev, c(82.5, 10, 2, .5, 10, .5, 2)*2)
summary(xx2)


# Now for a set of non-equal weights
nn <- nrow(mystate)
wts <- sample(1:5, nn, replace=T)
temp <- rep(1:nn, wts)             #row replicates
xgrp <- rep(1:10, length=nn)[order(runif(nn))]
xgrp2<- rep(xgrp, wts)
tempc <- rpart.control(minsplit=2, xval=xgrp2, maxsurrogate=0)
#  Direct: replicate rows in the data set, and use unweighted
fit5 <-  rpart(income ~ population + region + illiteracy +life + murder +
                        hs.grad + frost , data=mystate[temp,], control=tempc)
#  Weighted
tempc <- rpart.control(minsplit=2, xval=xgrp, maxsurrogate=0)
fit5b <-  rpart(income ~ population + region + illiteracy +life + murder +
                        hs.grad + frost , data=mystate, control=tempc,
                        weights=wts)
all.equal(fit5$frame[-2],  fit5b$frame[-2])  # the "n" component won't match
all.equal(fit5$cptable,    fit5b$cptable)
all.equal(fit5$splits[,-1],fit5b$splits[,-1]) #fails
all.equal(fit5$csplit,    fit5b$csplit)
