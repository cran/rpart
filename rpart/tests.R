#  version for use with R.

### This file tests many of the rpart/tree functions.  Plots are
### commented out so that the program can be run via BATCH.
###
### The assumption is that you are running this in the rpart 
### directory that you created.  If you are running this from another 
### directory, the first line of the program will need to be modified.
###

library(survival4)

kyphosis <- read.table("kyphosis.dat", header=T)
match <- 
function (x, table, nomatch = NA) 
.Internal(match(as.character(x), table, nomatch))

cu.summary <- read.table("cu.dat", header=T, sep=",", row.names=1)
cu.summary$Reliability <-
  as.ordered(factor(as.character(cu.summary$Reliability),
                    levels=c("Much worse","worse","average","better","Much better")))
cu.summary$Type <- factor(cu.summary$Type,
                          levels = c("Compact","Large","Medium","Small","Sporty","Van"))
solder.balance <- read.table("solder.dat", header=T)
solder.balance$Opening <- factor(solder.balance$Opening,
                                  levels = c("S", "M", "L"), ordered=T)
solder.balance$Solder <- factor(solder.balance$Solder,
                                 levels=c("Thin", "Thick"), ordered=T)
solder.balance$PadType <- factor(solder.balance$PadType, levels=
  c("W4","D4","L4","D6","L6","D7","L7","L8","W9","L9"))
solder.balance$Panel <- factor(solder.balance$Panel)

lung <- read.table("lung.dat", header=T)

library(rpart)

# Read in test data stagec
#
#   Time to progression in years
#   status  1=progressed, 0= censored
#   age
#   early endocrine therapy   1=no 2=yes
#   % of cells in g2 phase, from flow cytometry
#   tumor grade (Farrow) 1,2,3,4
#   Gleason score (competing grading system)
#   ploidy

stagec <- read.table("data.stagec",  col.names=c("pgtime", "pgstat", "age",
                     "eet", "g2", "grade", "gleason", "ploidy"))

stagec$ploidy <- factor(stagec$ploidy, levels=1:3,
                        labels=c("diploid", "tetraploid", "aneuploid"))

######################################################################


#.Random.seed <- c(17, 30, 53, 49, 12,  2, 22, 34, 23, 61, 29,  2)

# fit classification tree to data in kyphosis data frame
# method = "class"
 
c.fit <- rpart(Kyphosis ~ Age + Number + Start, data=kyphosis, 
               control=rpart.control(minsplit=10,cp=0.001,xval=5),
               parms=list(prior=c(.7,.3), 
                          loss=matrix(c(0,1,2,0),nrow=2,ncol=2)))
  
printcp(c.fit)
print(c.fit)
summary(c.fit,cp=.01)
summary(residuals(c.fit))
apply(xpred.rpart(c.fit),2,table)


######################################################################

#.Random.seed <- c(17, 30, 53, 49, 12,  2, 22, 34, 23, 61, 29,  2)

# fit classification tree to data in stagec data frame
# method = "class"

c.fit2 <- rpart(factor(pgstat) ~  age + eet + g2+grade+gleason +ploidy,
                data=stagec)

printcp(c.fit2)
print(c.fit2)
summary(c.fit2)
summary(xpred.rpart(c.fit2))

c.fit3 <- rpart(factor(pgstat) ~  age + eet + g2+grade+gleason +ploidy,
                stagec, parm=list(prior=c(.5,.5)))

printcp(c.fit3)
print(c.fit3)
print(c.fit3,cp=.05)
summary(c.fit3)
summary(c.fit3,cp=.05)

######################################################################

#.Random.seed <- c(17, 30, 53, 49, 12,  2, 22, 34, 23, 61, 29,  2)

# fit classification tree to data in cu.summary data frame
# method = "class" (multiple levels, endpoint is a factor)

## NOTE: If you have a large dataset, turning the crossvalidation and surrogate
##       splits off will speed things up.

c.fit4 <- rpart(Reliability ~ Price + Country + Mileage + Type, 
                data=cu.summary, control=rpart.control(xval=0,maxsurrogate=0))
  
printcp(c.fit4)
print(c.fit4)
print(c.fit4,cp=.025)
summary(c.fit4)
summary(c.fit4,cp=.025)
summary(residuals(c.fit4))
table(xpred.rpart(c.fit4))

######################################################################

# fit regression tree to all variables
# method = "anova"

#.Random.seed <- c(17, 30, 53, 49, 12,  2, 22, 34, 23, 61, 29,  2)
a.fit <- rpart(skips ~ Opening + Solder + Mask + PadType + Panel, 
               data=solder.balance, 
               control=rpart.control(xval=10))  

printcp(a.fit)
print(a.fit)
print(a.fit,cp=.10)
summary(a.fit)
summary(a.fit,cp=.10)
summary(residuals(a.fit))
  
xmat <- xpred.rpart(a.fit)
xerr <- (xmat - solder.balance$skips)^2
xsum <- apply(xerr,2,sum)/var(solder.balance$skips)  
xsum/xsum[1]  ## compare to rel error from printcp(a.fit)


## Does pruning work?
prune(a.fit, cp=.25)

## Check predict function
summary(predict(a.fit))
summary(solder.balance$skips)

summary(predict(a.fit,newdata=solder.balance[1,]))

######################################################################


# fit poisson tree to all variables
# method = "poisson"

#.Random.seed <- c(17, 30, 53, 49, 12,  2, 22, 34, 23, 61, 29,  2)
p.fit <- rpart(cbind(time,status) ~ inst + age + sex + ph.ecog + ph.karno +
               pat.karno + meal.cal + wt.loss, method="poisson", data=lung)

  
printcp(p.fit)
print(p.fit)
summary(p.fit)
summary(p.fit,cp=.025)
summary(residuals(p.fit))
apply(xpred.rpart(p.fit),2,summary)

######################################################################


# fit poisson tree to all variables rescaling time 
# (more like survival results)
# method = "exp"

#.Random.seed <- c(17, 30, 53, 49, 12,  2, 22, 34, 23, 61, 29,  2)

e.fit <- rpart(Surv(time,status) ~ inst + age + sex + ph.ecog + ph.karno +
               pat.karno + meal.cal + wt.loss, method="exp", data=lung)
  
printcp(e.fit)
print(e.fit)
print(e.fit,cp=.025)
summary(e.fit)
summary(e.fit,cp=.025)
summary(residuals(e.fit))
apply(xpred.rpart(e.fit),2,mean)

######################################################################

# fit poisson tree to all variables rescaling time 
# (more like survival results) Example #2
# method = "exp" and method="poisson"


fit1 <- rpart(Surv(pgtime, pgstat) ~ age + eet + g2+grade+gleason +ploidy,
                data=stagec, control=rpart.control(usesurrogate=0, cp=0),
                method="poisson")

printcp(fit1)
print(fit1)
summary(fit1,cp=.05)


fit2 <- rpart(Surv(pgtime, pgstat) ~ age + eet + g2+grade+gleason +ploidy,
                data=stagec, control=rpart.control(usesurrogate=1, cp=.001))

printcp(fit2)
print(fit2)
print(fit2,cp=.01)
summary(fit2,cp=.025)
 


###########################################################################


if(F) {

### the following code shows various ways to plot your results

plot(c.fit)
text(c.fit)

plot(c.fit2)
text(c.fit2,use.n=T)

plot(c.fit4,uniform=T)
text(c.fit4,use.n=T)

plot(a.fit,uniform=T)
text(a.fit,use.n=T,all=T)

plot(p.fit,branch=.5,compress=T,uniform=T)
text(p.fit,digits=2,use.n=T)

plot(e.fit,compress=T)
text(e.fit,digits=2)

par(mfrow=c(2,1))
rsq.rpart(a.fit)
par(mfrow=c(1,1))

## post.rpart creates a postscript file with a nicer looking tree
 
post.rpart(p.fit,file="",use.n=F) 
post.rpart(a.fit,file="",use.n=T) 
post.rpart(c.fit,title=" ",file="",use.n=T)
post(c.fit4,file='')
post.rpart(e.fit,file="",use.n=T)

## create a trimmed tree using snip.rpart 
##  (double click on left mouse button to remove branch, 
##   click on middle button to end trimming)

plot(p.fit)
text(p.fit)
p2.fit <- snip.rpart(p.fit)
summary(p2.fit)
plot(p2.fit)
text(p2.fit,use.n=T)

## look at some other plotting functions ...

meanvar(a.fit)
plot(predict(a.fit),residuals(a.fit))
abline(h=0,lty=2)

}

rm(a.fit,p.fit,e.fit,c.fit,c.fit2,c.fit3,c.fit4,fit1,fit2,xmat,xerr,xsum)

