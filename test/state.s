#
# In order to compare the x-vals estimates of the mainline and S versions,
#  it is necessary that we use stratified xval sets (like the mainline
#  does).  

mystate <- data.frame(state.x77, region=factor(state.region))
names(mystate) <- c("population","income" , "illiteracy","life" ,
       "murder", "hs.grad", "frost",     "area",      "region")

xvals <- 1:nrow(mystate)
xvals[order(mystate$income)] <- rep(1:10, length=nrow(mystate))

mystate <- data.frame(state.x77, region=factor(state.region))
names(mystate) <- c("population","income" , "illiteracy","life" ,
       "murder", "hs.grad", "frost",     "area",      "region")

fit4 <- rpart(income ~ population + region + illiteracy +life + murder +
			hs.grad + frost , mystate,
		   control=rpart.control(minsplit=10, xval=xvals))

summary(fit4)


#
# Check out xpred.rpart
#
meany <- mean(mystate$income)
xpr <- xpred.rpart(fit4, xval=xvals)
xpr2 <- (xpr - mystate$income)^2
risk0 <- mean((mystate$income - meany)^2)
xpmean <- as.vector(apply(xpr2, 2, mean))   #kill the names
all.equal(xpmean/risk0, as.vector(fit4$cptable[,'xerror']))

xpstd <- as.vector(apply((sweep(xpr2, 2, xpmean))^2, 2, sum))
xpstd <- sqrt(xpstd)/(50*risk0)
all.equal(xpstd, as.vector(fit4$cptable[,'xstd']))

#
# recreate subset #3 of the xval
#
tfit4 <- rpart(income ~ population + region + illiteracy +life + murder +
			hs.grad + frost , mystate,  subset=(xvals!=3),
		   control=rpart.control(minsplit=10, xval=0))
tpred <- predict(tfit4, mystate[xvals==3,])
all.equal(tpred, xpr[xvals==3,ncol(xpr)])

# How much does this differ from the "real" formula, more complex,
#   found on page 309 of Breiman et al. ?
#xtemp <- (xpr2/outer(rep(1,50),xpmean)) -  ((mystate$income - meany)^2)/risk0
#real.se<- xpmean* sqrt(apply(xtemp^2,2,sum))/(risk0*50)

