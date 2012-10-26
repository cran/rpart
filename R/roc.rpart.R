roc.rpart <- function(object,plot.ok=TRUE,x.orient=1)
{

    if(!inherits(object, 'rpart') & !(object$method=='class') &
       !(length(attr(object,"ylevels")) ==2))
        stop('Not legitimate \"rpart\" tree and endpoint not a 2 level-factor')

    ss.compare <- function(a, b)	a >= b

    endnodes <- object$frame$splits[,1]==""
    truth <- object$frame$yval2[endnodes,]
    cutoffs <- sort(unique(c(0,1,object$frame$yprob[endnodes,2])))
    pred.np <- outer(cutoffs,object$frame$yprob[endnodes,2], ss.compare)

    last.r <- dim(pred.np)[1]
    last.c <- dim(pred.np)[2]
    if(sum(pred.np[1,])>0) {
        pred.np <- rbind(matrix(FALSE,nrow=1,ncol=last.c),pred.np)
	cutoffs <- c(NA,cutoffs)
    }
    last.r <- dim(pred.np)[1]
    last.c <- dim(pred.np)[2]
    if(sum(pred.np[last.r,])<last.c) {
        pred.np <- rbind(pred.np,matrix(TRUE,nrow=1,ncol=last.c))
	cutoffs <- c(cutoffs,NA)
    }

    cutoff.n <- length(cutoffs)
    ## set up some empty matrices ##
    sensitivity <- matrix(0,nrow=cutoff.n,ncol=1)
    specificity <- matrix(0,nrow=cutoff.n,ncol=1)
    negpred <- matrix(0,nrow=cutoff.n,ncol=1)
    pospred <- matrix(0,nrow=cutoff.n,ncol=1)
    tpcp<- matrix(0,nrow=cutoff.n,ncol=1)
    tncp<- matrix(0,nrow=cutoff.n,ncol=1)
    tpcn<- matrix(0,nrow=cutoff.n,ncol=1)
    tncn<- matrix(0,nrow=cutoff.n,ncol=1)
    ss.table <- array(0,c(cutoff.n,2,2))

    for (i in 1:cutoff.n) {

        ss.table <- matrix(0,nrow=2,ncol=2)

        ss.table[1,1] <- sum(truth[pred.np[i,],1])
        ss.table[2,1] <- sum(truth[!pred.np[i,],1])
        ss.table[1,2] <- sum(truth[pred.np[i,],2])
        ss.table[2,2] <- sum(truth[!pred.np[i,],2])

        sensitivity[i] <- ss.table[2,2]/(ss.table[2,2] + ss.table[1,2])
        specificity[i] <- ss.table[1,1]/(ss.table[1,1] + ss.table[2,1])
        negpred[i] <- ss.table[1,1]/(ss.table[1,1] + ss.table[1,2])
        pospred[i] <- ss.table[2,2]/(ss.table[2,2] + ss.table[2,1])
        tpcp[i]<-ss.table[2,2]
        tncp[i]<-ss.table[2,1]
        tpcn[i]<-ss.table[1,2]
        tncn[i]<-ss.table[1,1]
    }

    if(plot.ok) {

        o.par <- par(pty='s')
        on.exit(par(o.par))
        if(x.orient==1){
            plot(1-specificity,sensitivity,type='o', xlim=c(0,1),
                 ylim=c(0,1), ylab='Sensitivity', xlab='1-Specificity')
        }
        if(x.orient==2){
            plot(specificity,sensitivity,type='o', xlim=c(0,1),
                 ylim=c(0,1), ylab='Sensitivity', xlab='Specificity')
        }
    }

    data.frame(cutoffs=format(round(cutoffs,3)), tpcp,tncp,tpcn,tncn,
               sensitivity=format(round(sensitivity,2)),
               specificity=format(round(specificity,2)),
               pospred=format(round(pospred,2)),
               negpred=format(round(negpred,2)))


}

