fit5 <- rpart(factor(pgstat) ~  age + eet + g2+grade+gleason +ploidy,
	  stagec)

fit5

fit6 <- rpart(factor(pgstat) ~  age + eet + g2+grade+gleason +ploidy,
		stagec, parm=list(prior=c(.5,.5)))
summary(fit6)
