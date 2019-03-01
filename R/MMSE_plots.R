
plotquant<-function(x,p=c(0.05,0.25,0.75,0.95),yrs,qcol,lcol,addline=T,ablines=NA){
  ny<-length(yrs)
  qs<-apply(x,2,quantile,p=p[c(1,4)])
  qsi<-apply(x,2,quantile,p=p[2:3])
  polygon(c(yrs,yrs[ny:1]),c(qs[1,],qs[2,ny:1]),border=NA,col='#b3ecff')

  polygon(c(yrs,yrs[ny:1]),c(qsi[1,],qsi[2,ny:1]),border=NA,col=qcol)
  if(!is.na(ablines[1]))abline(h=ablines,col='#99999980')

  if(addline)for(i in 1:2)lines(yrs,x[i,],col=lcol,lty=i)
  lines(yrs,apply(x,2,quantile,p=0.5),lwd=2,col="white")
}

Pplot3<-function(MSEobj,maxcol=6,qcol=rgb(0.4,0.8,0.95), lcol= "dodgerblue4",curyr=2018,quants=c(0.05,0.25,0.75,0.95),MPcols='black',addline=F,maxrow=NA){

  if(is.na(maxcol))maxcol=ceiling(length(MSEobj@MPs)/0.5) # defaults to portrait 1:2
  MPs<-MSEobj@MPs
  nMPs<-length(MPs)
  yrs<-curyr+(1:MSEobj@proyears)

  MPcols[MPcols=="green"]<-"darkgreen"

  plots<-split(1:nMPs, ceiling(seq_along(1:nMPs)/maxcol))

  nr<-length(plots)*2
  if(!is.na(maxrow))nr=max(nr,maxrow)
  nc<-maxcol

  mat<-array(0,c(nc,nr*1.5))
  ind<-floor(0.5+(1:nr)*1.5)
  mat[,ind]<-1:(nr*nc)
  mat<-t(mat)
  ht<-rep(0.2,nr*1.5)
  ht[ind]<-1
  layout(mat,heights=ht)
  par(mai=c(0.3,0.3,0.01,0.01),omi=c(0.5,0.5,0.05,0.05))

  B_BMSY<-MSEobj@B_BMSY
  Yd<-MSEobj@C/ array(rep(MSEobj@C[,,1],MSEobj@proyears),dim(MSEobj@C))#MSEobj@OM$RefY
  Yd[is.na(Yd)]<-0

  Blims <- c(0,quantile(B_BMSY,0.95))
  Ylims<- c(0,quantile(Yd,0.95))

  for(pp in 1:length(plots)){

    toplot<-unlist(plots[pp])
    nt<-length(toplot)

    for(i in toplot){
      plot(range(yrs),Blims,col="white")
      plotquant(B_BMSY[,i,],p=quants,yrs,qcol,lcol,ablines=c(0.5,1),addline=addline)
      mtext(MSEobj@MPs[i],3,line=0.2,font=2,col=MPcols[i])
      if(i==toplot[1])mtext("B/BMSY",2,line=2.3)
    }
    if(nt<maxcol)for(i in 1:(maxcol-nt))plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="y label", xlab="x lablel",axes=F)

    for(i in toplot){
      plot(range(yrs),Ylims,col="white")
      plotquant(Yd[,i,],p=quants,yrs,qcol,lcol,ablines=1,addline=addline)
      if(i==toplot[1])mtext("Yield relative to today",2,line=2.3)
    }
    if(nt<maxcol)for(i in 1:(maxcol-nt))plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="y label", xlab="x lablel",axes=F)

  }

  mtext("Projection Year",1,line=0.7,outer=T)

}


BY<-function(MMSEobj,maxcol=6,qcol=rgb(0.4,0.8,0.95), lcol= "dodgerblue4",curyr=2018,quants=c(0.05,0.25,0.75,0.95),MPcols='black',addline=F,maxrow=NA){

  if(is.na(maxcol))maxcol=ceiling(length(MSEobj@MPs)/0.5) # defaults to portrait 1:2
  MPs<-MSEobj@MPs
  nMPs<-length(MPs)
  yrs<-curyr+(1:MSEobj@proyears)

  MPcols[MPcols=="green"]<-"darkgreen"

  plots<-split(1:nMPs, ceiling(seq_along(1:nMPs)/maxcol))

  nr<-length(plots)*2
  if(!is.na(maxrow))nr=max(nr,maxrow)
  nc<-maxcol

  mat<-array(0,c(nc,nr*1.5))
  ind<-floor(0.5+(1:nr)*1.5)
  mat[,ind]<-1:(nr*nc)
  mat<-t(mat)
  ht<-rep(0.2,nr*1.5)
  ht[ind]<-1
  layout(mat,heights=ht)
  par(mai=c(0.3,0.3,0.01,0.01),omi=c(0.5,0.5,0.05,0.05))

  B_BMSY<-MSEobj@B_BMSY
  Yd<-MSEobj@C/ array(rep(MSEobj@C[,,1],MSEobj@proyears),dim(MSEobj@C))#MSEobj@OM$RefY
  Yd[is.na(Yd)]<-0

  Blims <- c(0,quantile(B_BMSY,0.95))
  Ylims<- c(0,quantile(Yd,0.95))

  for(pp in 1:length(plots)){

    toplot<-unlist(plots[pp])
    nt<-length(toplot)

    for(i in toplot){
      plot(range(yrs),Blims,col="white")
      plotquant(B_BMSY[,i,],p=quants,yrs,qcol,lcol,ablines=c(0.5,1),addline=addline)
      mtext(MSEobj@MPs[i],3,line=0.2,font=2,col=MPcols[i])
      if(i==toplot[1])mtext("B/BMSY",2,line=2.3)
    }
    if(nt<maxcol)for(i in 1:(maxcol-nt))plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="y label", xlab="x lablel",axes=F)

    for(i in toplot){
      plot(range(yrs),Ylims,col="white")
      plotquant(Yd[,i,],p=quants,yrs,qcol,lcol,ablines=1,addline=addline)
      if(i==toplot[1])mtext("Yield relative to today",2,line=2.3)
    }
    if(nt<maxcol)for(i in 1:(maxcol-nt))plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="y label", xlab="x lablel",axes=F)

  }

  mtext("Projection Year",1,line=0.7,outer=T)

}



multidebug<-function(MSEsingle,MSEmulti,p=1,f=1,MPno=1,maxsims=4){

  MPno=1
  p=1
  f=1
  par(mfrow=c(2,2),mar=rep(2,4))

  SSBhs<-apply(slot(MSEsingle,"SSB_hist"),c(1,3),sum)
  SSBhm<-apply(slot(MSEmulti,"SSB_hist")[,p,,,],c(1,3),sum)
  matplot(t(SSBhs),type='l',main="SSBhist: runMSE")
  matplot(t(SSBhm),type='l',main="SSBhist: multiMSE")

  SSBs<-slot(MSEsingle,"SSB")
  SSBm<-slot(MSEmulti,"SSB")
  matplot(t(SSBs[,MPno,]),type='l',main="SSB: runMSE")
  matplot(t(SSBm[,p,MPno,]),type='l',main="SSB: multiMSE")

}




