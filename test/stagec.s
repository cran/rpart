#
# Read the data
#
#   Time to progression in years
#   status  1=progressed, 0= censored
#   age
#   early endocrine therapy   1=no 2=yes
#   % of cells in g2 phase, from flow cytometry
#   tumor grade (Farrow) 1,2,3,4
#   Gleason score (competing grading system)
#   ploidy

stagec <- read.table('data.stagec',  col.names=c("pgtime", "pgstat", "age",
			"eet", "g2", "grade", "gleason", "ploidy"))
stagec$ploidy <- factor(stagec$ploidy, levels=1:3,
				labels=c("diploid", "tetraploid", "aneuploid"))

cox0 <- coxph(Surv(pgtime, pgstat) ~ 1, stagec)
cox1 <- coxph(Surv(pgtime, pgstat) ~ age + eet + g2 + grade + ploidy, stagec)
cox1

fit1 <- rpart(Surv(pgtime, pgstat) ~ age + eet + g2+grade+gleason +ploidy,
		stagec, control=rpart.control(usesurrogate=0, cp=0),
		method='poisson')
fit1
summary(fit1)

fit2 <- rpart(cox0$residual ~ age + eet + g2+grade+gleason +ploidy,
		stagec)
fit2
summary(fit2)



fit3 <- rpart(Surv(pgtime, pgstat) ~ age + eet + g2+grade+gleason +ploidy,
		stagec, control=rpart.control(usesurrogate=1, cp=.001))

summary(fit3)

