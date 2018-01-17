# ===================================================================
# === Ancilary indicators ===========================================
# ===================================================================

slp<-function(x,mat,ind){

  lm(y~x1,data.frame(x1=1:length(ind),y=log(mat[x,ind])))$coef[2]

}

AAV<-function(x,mat,ind){
  ni<-length(ind)
  mean(abs((mat[x,ind[2:ni]]-mat[x,ind[1:(ni-1)]])/mat[x,ind[1:(ni-1)]]))
}

mu<-function(x,mat,ind){
  mean(mat[x,ind])
}

getinds<-function(PPD,styr,res){

  nsim<-dim(PPD@Cat)[1]
  proyears<-dim(PPD@Cat)[2]-styr+1

  if(res>proyears)message(paste0("The temporal resolution for posterior predictive data calculation (",res,") is higher than the number of projected years (",proyears,"). Only one time step of indicators are calculated for ",proyears, " projected years."))
  np<-floor(proyears/res)

  tsd<- c("Cat","Cat","Cat", "Ind","ML")     #, "Lc")
  stat<-c("slp","AAV","mu",  "slp",  "slp")  #,"slp")
  ntsd<-length(tsd)
  inds<-array(NA,c(ntsd,np,nsim))

  for(i in 1:ntsd){
    for(pp in 1:np){
      ind<-styr+((pp-1)*res)+1:res
      inds[i,pp,]<-sapply(1:nsim,get(stat[i]),mat=slot(PPD,tsd[i]),ind=ind)
    }
  }

  inds

}

CC<-function(indPPD,indData,pp=1,nams=c("YS", "AAVY","YM","IS","MLS"),res){

  namst<-paste(rep(nams,pp),rep((1:pp)*res,each=length(nams)))
  cols<-c("#ff000050","#0000ff50")
  ntsd<-dim(indPPD)[1]
  ni<-pp*ntsd
  ind2PPD<-matrix(indPPD[,1:pp,],nrow=ni)
  ind2Data<-matrix(indData[,1:pp,],nrow=ni)
  par(mfrow=c(ni-1,ni-1),mai=rep(0,4),omi=c(0.5,0.75,0,0.05))

  for(i in 2:ni){

    for(j in 1:(ni-1)){

      if(j==i|j>i){

        plot(1,1,col='white',axes=F)

      }else{

        #coly=cols[ceiling(posmean(cor(mcmc@rawdat[1:maxn,keep1[i]],mcmc@rawdat[1:maxn,keep2[j]]))*ncols)]
        xlim<-quantile(c(ind2PPD[i,],ind2Data[i,]),c(0.02,0.98))
        ylim<-quantile(c(ind2PPD[j,],ind2Data[j,]),c(0.02,0.98))
        plot(ind2PPD[i,],ind2PPD[j,],pch=19,xlim=xlim,ylim=ylim,cex=0.8,col=cols[1],axes=F)
        points(ind2Data[i,],ind2Data[j,],pch=19,cex=0.8,col=cols[2])

      }

      if(j==1)mtext(namst[i],2,line=2,cex=0.6,las=2)
      if(i==ni)mtext(namst[j],1,line=1,cex=0.6,las=2)
      #if(j==1)mtext(i,2,line=2,cex=0.5,las=2)
      #if(i==nplotted)mtext(j,1,line=1,cex=0.5,las=2)

    }

  }

}

Probs<-function(indPPD,indData,alpha=0.05){

  np<-dim(indPPD)[2]
  nsim<-dim(indPPD)[3]
  #PRB<-array(NA,c(4,np)) # True Negative, False Negative, False Positive, True Positive
  PRB<-array(NA,c(2,np)) # True Negative, False Negative, False Positive, True Positive
  mah<-array(NA,c(2,np,nsim))

  for(pp in 1:np){

    ntsd<-dim(indPPD)[1]
    ni<-pp*ntsd

    ind3PPD<-t(matrix(indPPD[,1:pp,],nrow=ni))
    ind3Data<-t(matrix(indData[,1:pp,],nrow=ni))

    # NULL = TRUE  (true negatives, false negatives)
    covr <- cov.mcd(ind3PPD)
    mahN <- mahalanobis(ind3PPD, center = covr$center, cov = covr$cov, tol = 1e-25)
    mahA <- mahalanobis(ind3Data, center = covr$center, cov = covr$cov, tol = 1e-25)

    mah[1,pp,]<-mahN
    mah[2,pp,]<-mahA

    Thres<-quantile(mahN,1-alpha)
    PRB[1,pp]<-sum(mahN>Thres)/length(mahN)   # False positive
    PRB[2,pp]<-sum(mahA>Thres)/length(mahA)   # True positive

  }

  list(mah=mah,PRB=PRB)

}

#' Exceptional Circumstances
#'
#' @param MSE An object of class MSE
#' @param hzn Time horizon for posterior data
#' @param alpha Probability of incorrectly rejecting the null hypothesis of normal data when it is true
#' @importFrom MVN mardiaTest roystonTest hzTest
#' @importFrom robustbase covMcd
#' @return Indicators of MSE misspecification
#' @author T. Carruthers
#' @references Carruthers et al. 2018
#' @export
EC=function(PPD,Data,res=6,styr=NA,LB=0.9,UB=0.95,plot=T){

  indPPD<-getinds(PPD,styr=styr,res=res)
  indData<-getinds(Data,styr=styr,res=res)

  if(plot)CC(indPPD,indData,pp=2,res=res)

  out<-Probs(indPPD,indData)
  np<-dim(indPPD)[2]
  par(mfrow=c(np,1))

  for(pp in 1:np){

    xlim=quantile(out$mah[,pp,],c(0.01,0.99))
    densN<-density(out$mah[1,pp,],adj=1,from=0)
    densN$y<-densN$y/max(densN$y)
    densA<-density(out$mah[2,pp,],adj=1,from=0)
    densA$y<-densA$y/max(densA$y)

    plot(densN,xlim=xlim,col='blue',main="")
    lines(densA,xlim=xlim,col='red')

  }


}




# ====================================================================
