# ===================================================================
# === Ancilary indicators ===========================================
# ===================================================================

slp<-function(x,mat,ind){

  lm(y~x1,data.frame(x1=1:length(ind),y=log(mat[x,ind])))$coef[2]

}


slp2<-function(x,mat,ind){

  x1<-1:length(ind)
  y=log(mat[x,ind])
  mux<-mean(x1)
  muy<-mean(y)
  SS<-sum((x1-mux)^2)
  (1/SS)*sum((x1-mux)*(y-muy))

}

AAV<-function(x,mat,ind){
  ni<-length(ind)
  mean(abs((mat[x,ind[2:ni]]-mat[x,ind[1:(ni-1)]])/mat[x,ind[1:(ni-1)]]))
}

mu<-function(x,mat,ind){
  log(mean(mat[x,ind]))
}

getinds<-function(PPD,styr,res, tsd= c("Cat","Cat","Cat","Ind","ML"),stat=c("slp","AAV","mu","slp", "slp")){

  nsim<-dim(PPD@Cat)[1]
  proyears<-dim(PPD@Cat)[2]-styr+1

  if(res>proyears)message(paste0("The temporal resolution for posterior predictive data calculation (",res,") is higher than the number of projected years (",proyears,"). Only one time step of indicators are calculated for ",proyears, " projected years."))
  np<-floor(proyears/res)

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

CC<-function(indPPD,indData,pp=1,dnam=c("CS","CV","CM","IS","IM","MLS","MLM"),res){

  if(pp>1)namst<-paste(rep(dnam,pp),rep((1:pp)*res,each=length(dnam)))
  if(pp==1)namst=dnam
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
        xlim<-quantile(c(ind2PPD[j,],ind2Data[j,]),c(0.02,0.98))
        ylim<-quantile(c(ind2PPD[i,],ind2Data[i,]),c(0.02,0.98))
        plot(ind2PPD[j,],ind2PPD[i,],pch=19,xlim=xlim,ylim=ylim,cex=0.8,col=cols[1],axes=F)
        points(ind2Data[j,],ind2Data[i,],pch=19,cex=0.8,col=cols[2])

      }
      if(i==2&j==(ni-1)){
        legend('center',legend=c("Null - stable M", "Alternative - inc M"),text.col=c("blue","red"),bty='n')

      }

      if(j==1)mtext(namst[i],2,line=2,cex=0.6,las=2)
      if(i==ni)mtext(namst[j],1,line=1,cex=0.6,las=2)
      #if(j==1)mtext(i,2,line=2,cex=0.5,las=2)
      #if(i==nplotted)mtext(j,1,line=1,cex=0.5,las=2)

    }

  }

}

Probs<-function(indPPD,indData,alpha=0.05){

  ntsd<-dim(indPPD)[1]
  np<-dim(indPPD)[2]
  nsim<-dim(indPPD)[3]
  #PRB<-array(NA,c(4,np)) # False Positive, True Positive
  PRB<-array(NA,c(2,np))  # False Positive, True Positive
  mah<-array(NA,c(2,np,nsim))

  for(pp in 1:np){

    keep<-array(TRUE,c(ntsd,pp))
    for(i in 1:ntsd){
     for(j in 1:pp){
       if(dip(indPPD[i,j,])>0.065)keep[i,j]=FALSE
      }
    }

    ni<-sum(keep)
    keepind<-as.matrix(expand.grid(1:ntsd,1:pp,1:nsim))[rep(as.vector(keep),nsim),]

    ind3PPD<-t(matrix(indPPD[keepind],nrow=ni))
    ind3Data<-t(matrix(indData[keepind],nrow=ni))

    # NULL = TRUE  (true negatives, false negatives)
    #covr <- cov.mcd(ind3PPD)
    covr <- cov(ind3PPD)
    #test<-svd(covr)
    mu<-apply(ind3PPD,2,median)
    #mahN <- mahalanobis(ind3PPD, center = covr$center, cov = covr$cov, tol = 1e-25)
    #mahA <- mahalanobis(ind3Data, center = covr$center, cov = covr$cov, tol = 1e-25)
    #mahN1 <- mahalanobis(ind3PPD, center = mu, cov = covr, tol = 1e-25)
    #mahA1 <- mahalanobis(ind3Data, center = mu, cov = covr, tol = 1e-25)

    mahN <- mahalanobis_robust(ind3PPD, center = mu, cov = covr)
    mahA <- mahalanobis_robust(ind3Data, center = mu, cov = covr)
    mahN<-extreme.outlier(mahN)
    mahA<-extreme.outlier(mahA)

    mah[1,pp,]<-mahN
    mah[2,pp,]<-mahA

    Thres<-quantile(mahN,1-alpha,na.rm=T)
    #plot(density(mah[1,pp,],na.rm=T),type="l",col="blue")
    #lines(density(mah[2,pp,],na.rm=T),col="red")
    #abline(v=Thres)


    PRB[1,pp]<-mean(mah[1,pp,]>Thres,na.rm=T)   # False positive
    PRB[2,pp]<-mean(mah[2,pp,]>Thres,na.rm=T)   # True positive

  }

  return(list(mah=mah,PRB=PRB,keep=keep))

}

#' Exceptional Circumstances
#'
#' @param MSE An object of class MSE
#' @param hzn Time horizon for posterior data
#' @param alpha Probability of incorrectly rejecting the null hypothesis of normal data when it is true
#' @importFrom MASS cov.mcd
#' @importFrom corpcor pseudoinverse
#' @return Indicators of MSE misspecification
#' @author T. Carruthers
#' @references Carruthers et al. 2018
#' @export
PRBcalc=function(MSE_null,MSE_alt,
                 tsd= c("Cat","Cat","Cat","Ind","ML"),
                 stat=c("slp","AAV","mu","slp", "slp"),
                 dnam=c("CS","CV","CM","IS","MLS"),
                 res=6,alpha=0.05,plotCC=F){

  styr<-MSE_null@nyears+1

  outlist<-new('list')

  for(mm in 1:MSE_null@nMPs){

    PPD<-MSE_null@Misc[[mm]]
    Data<-MSE_alt@Misc[[mm]]

    indPPD<-getinds(PPD,styr=styr,res=res,tsd=tsd,stat=stat)
    indData<-getinds(Data,styr=styr,res=res,tsd=tsd,stat=stat)

    if(plotCC)CC(indPPD,indData,pp=2,res=res)

    out<-Probs(indPPD,indData,alpha=alpha)

    outlist[[mm]]<-out

    message(paste(mm, "of",MSE_null@nMPs,"MPs processed"))

  }

  return(outlist)

}

#' @importFrom corpcor pseudoinverse
mahalanobis_robust<-function (x, center, cov, inverted = FALSE) {

  x <- if (is.vector(x))
    matrix(x, ncol = length(x))
  else as.matrix(x)
  if (!identical(center, FALSE))
    x <- sweep(x, 2L, center)

  invcov <- pseudoinverse(cov)
  setNames(rowSums(x %*% invcov * x), rownames(x))

}



mahplot<-function(outlist,MPs,res=6){

  nMP<-length(outlist)
  ncol<-floor(nMP^(1/3))
  nrow<-ceiling(nMP/ncol)
  par(mai=c(0.2,0.2,0.2,0.01),omi=c(0.2,0.6,0.01,0.01))
  layout(matrix(1:(nrow*ncol*2),nrow=nrow*2),heights=rep(c(2,1),nrow))
  pmin<-min(c(0.95,sapply(outlist,function(x)min(x$PRB[2,]))))
  plabs<-matrix(paste0("(",letters[1:(nMP*2)],")"),nrow=2)

  for(mm in 1:nMP){
    xaxis=F
    yaxis=F
    if(mm<(nMP/ncol)|mm==(nMP/ncol))yaxis=T
    if(mm%in%((1:ncol)*(nMP/ncol)))xaxis=T

    mahdensplot(outlist[[mm]],xaxis=F,yaxis=yaxis,res=res)
    if(yaxis)mtext("Distance D",2,cex=0.9,line=3.5)
    mtext(MPs[mm],3,font=2,line=0.3,cex=0.95)
    mtext(plabs[1,mm],3,line=0.07,adj=0.01,cex=0.88)

    PRBplot(outlist[[mm]],xaxis=xaxis,yaxis=yaxis,res=res,ylim=c(pmin,1))
    if(yaxis)mtext("Power",2,cex=0.9,line=3.5)
    mtext(plabs[2,mm],3,line=0.07,adj=0.01,cex=0.88)

  }

}

mahdensplot<-function(out,adj=0.9,alpha=0.05,xaxis=F,yaxis=F,res=6){

  heightadj<-0.7
  np<-dim(out$PRB)[2]

  densN<-new('list')
  densA<-new('list')
  ylimy=array(0,c(np,2))
  thresh<-rep(0,np)
  tfac<-4
  xref<-rep(0,np)

  xmaxN<-xmaxA<-rep(0,np)


  for(pp in 1:np){
    thresh[pp]<-quantile(out$mah[1,pp,!is.na(out$mah[1,pp,])],1-alpha)
    xref[pp]<-thresh[pp]*tfac
    ylimy[pp,2]<-xref[pp]#quantile(out$mah[1,pp,],0.95,na.rm=T)*1.2
    densN[[pp]]<-density(out$mah[1,pp,],adj=adj,from=0,to=xref[pp],na.rm=T)
    densA[[pp]]<-density(out$mah[2,pp,],adj=adj,from=0,to=xref[pp],na.rm=T)
    ymaxN<-max(densN[[pp]]$y)
    ymaxA<-max(densA[[pp]]$y)
    xmaxN[pp]=ymaxN
    xmaxA[pp]=ymaxA
  }

  for(pp in 1:np){
    #densN[[pp]]$y<-densN[[pp]]$y/max(densN[[pp]]$y)*heightadj
    #densA[[pp]]$y<-densA[[pp]]$y/max(densA[[pp]]$y)*heightadj
    densN[[pp]]$y<-densN[[pp]]$y/xmaxN[pp]*heightadj
    densA[[pp]]$y<-densA[[pp]]$y/xmaxN[pp]*heightadj
  }

  plot(c(0.5,np+0.1),c(0,tfac),col='white',axes=F)

  if(xaxis){
    labs<-paste("Yrs",(1:np)*res-(res-1),"-",(1:np)*res)
    axis(1,1:np,labs)

  }
  if(yaxis){
    axis(2,c(0,1,10E10),c("0","V",""))
  }

  lines(c(-100,100),rep(1,2),lty=2,lwd=1)
  for(pp in 1:np){

    trs<-0.5+pp-1

    lines(densN[[pp]]$y+trs,densN[[pp]]$x/thresh[pp],col='blue')
    lines(densA[[pp]]$y+trs,densA[[pp]]$x/thresh[pp],col='red')

    polyN<-getsegment(densN[[pp]],thresh[pp],lower=T)
    polyA<-getsegment(densA[[pp]],thresh[pp],lower=F)
    polygon(polyN$x+trs,polyN$y/thresh[pp],border=NA,col="#0000ff50")

    if(length(polyA$x)>10)polygon(polyA$x+trs,polyA$y/thresh[pp],border=NA,col="#ff000050")

  }
}


PRBplot<-function(out,res,xaxis=F,yaxis=F,ylim=c(0,1)){
  np<-dim(out$PRB)[2]
  plot(c(0.5,np+0.5),range(out$PRB),col='white',axes=F,xlab="",ylab="",ylim=ylim)
  #lines((1:np)-0.2,1-out$PRB[1,],col="#0000ff50",lwd=2)
  abline(h=0.8,lty=2)
  lines((1:np)-0.2,out$PRB[2,],col="#ff000050",lwd=2)
  axis(1,at=c(-1000,10000))
  if(xaxis){
    labs<-paste0("Yrs ",(1:np)*res-(res-1),"-",(1:np)*res)
    at = axTicks(1)
    mtext(side = 1, text = labs, at = at, line = 1,cex=0.9)

    #axis(1,1:np,labs)
  }
  axis(2)
  #if(yaxis)axis(2)

}


getsegment<-function(densobj,thresh,lower=T){

  if(lower){
    cond<-densobj$x<thresh
  }else{
    cond<-densobj$x>thresh
  }

  xs<-c(0,densobj$y[cond],0)
  ys<-densobj$x[cond]
  ys<-c(ys[1],ys,ys[length(ys)])

  list(x=xs,y=ys)

}

extreme.outlier<-function(x){
  sd<-mean(abs(x-median(x)))
  cond<-(x<(median(x)-3*sd(x)))|(x>(median(x)+3*sd(x)))
  x[cond]<-NA
  x
}

# ====================================================================
