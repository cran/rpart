#
# Test out some Poisson xpreds
#
xx <- xpred.rpart(fit1)

xevent <- xx * stagec$pgtime   # expected number of events for the subject

temp0 <- apply( stagec$pgstat - xevent, 2, mean)
temp1 <- apply((stagec$pgstat - sqrt(xevent))^2 ,2, mean)
temp2 <- apply((stagec$pgstat - xevent)^2/xevent, 2, mean)

temp4 <- cbind(fit1$cptable, temp0, temp1/temp1[1], temp2/temp2[1])
dimnames(temp4) <- list(NULL, c(dimnames(fit1$cptable)[[2]], 'bias',
	'scaled', 'Pearson'))

