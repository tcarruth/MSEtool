# ===================================================================
# === Ancilary indicators ===========================================
# ===================================================================

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
EC=function(MSE,res=6,MP=1,LB=0.9,UB=0.95){

  PPD<-MSE@Misc[[MP]]
  nsim<-MSE@nsim
  nyears<-MSE@nyears
  proyears<-dim(PPD@Cat)[2]-nyears

  if(res>proyears)message(paste0("The temporal resolution for posterior predictive data calculation (",res,") is higher than the number of projected years (",proyears,"). Only one time step of indicators are calculated for ",proyears, " projected years."))
  np<-floor(proyears/res)

  tsd<-c("Cat","Ind","ML","Lc")
  ntsd<-length(tsd)
  inds<-array(NA,c(ntsd,np,nsim))

  for(i in 1:ntsd){
    for(pp in 1:np){
      ind<-nyears+((pp-1)*res)+1:res
      inds[i,pp,]<-sapply(1:nsim,getgrad,mat=slot(PPD,tsd[i]),ind=ind)
    }
  }

  pp<-1
  tol = 1e-25
  dat<-t(inds[,pp,])
  covr <- covMcd(dat)
  mah <- mahalanobis(dat, center = covr$center, cov = covr$cov, tol = tol)


  ncols <- 200  # nrow(Stat) * ncol(Stat)
  #refcol <- colorRampPalette(c("green", "red"))(ncols)
  refcol <- rainbow(ncols, start = 0, end = 0.36)[ncols:1]
  cols<-rep(refcol[1],nsim)
  LB<-quantile(mah,0.8)
  UB<-quantile(mah,0.95)
  cols[mah>UB]<-refcol[ncols]
  cond<-mah>LB&mah<UB
  cols[cond]<-refcol[ceiling((mah[cond]-LB)/(UB-LB)*ncols)]
  pchs<-rep(19,nsim)
  pchs[mah>UB]<-4
  lwds<-rep(1,nsim)
  lwds[mah>UB]<-2
  #mardiaTest(dat,qqplot=T)
  #roystonTest(dat)
  #hzTest(dat)

  pp<-1
  par(mfrow=c(3,3),mar=rep(4,4))
  plot(inds[1,pp,],inds[2,pp,],xlab="Cat grad",ylab="Ind grad",col=cols,pch=pchs,lwd=lwds)
  plot(inds[1,pp,],inds[3,pp,],xlab="Cat grad",ylab="ML grad",col=cols,pch=pchs,lwd=lwds)
  plot(inds[1,pp,],inds[4,pp,],xlab="Cat grad",ylab="Lc grad",col=cols,pch=pchs,lwd=lwds)
  plot(1,1,col="white",axes=F)
  plot(inds[2,pp,],inds[3,pp,],xlab="Ind grad",ylab="ML grad",col=cols,pch=pchs,lwd=lwds)
  plot(inds[2,pp,],inds[4,pp,],xlab="Ind grad",ylab="Lc grad",col=cols,pch=pchs,lwd=lwds)
  plot(1,1,col="white",axes=F)
  plot(1,1,col="white",axes=F)
  plot(inds[3,pp,],inds[4,pp,],xlab="ML grad",ylab="Lc grad",col=cols,pch=pchs,lwd=lwds)

  (1:nsim)[mah>UB]
}



getgrad<-function(x,mat,ind){

  lm(y~x1,data.frame(x1=1:length(ind),y=log(mat[x,ind])))$coef[2]

}




# ====================================================================
