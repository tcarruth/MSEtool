# multiMSE funcs

# Outdated function for producing a stock-specific historical reconstruction
HistMulti<-function(x,FleetP,StockP,maxage,allyears,nyears,nf,MPAc,qfrac, qs,maxF){

  Vuln=retA=array(NA,c(nf,maxage,nyears))
  Effind<-t(matrix(unlist(lapply(FleetP,function(dat,x)dat['Find'][[1]][x,],x=x)),ncol=nf))
  Spat_targc<-unlist(lapply(FleetP,function(dat,x)dat['Spat_targ'][[1]][x],x=x))

  for(ff in 1:nf){ # kind of lazy but not THAT slow
    Vuln[ff,,]<-FleetP[[ff]]$V[x,,1:nyears]
    retA[ff,,]<-FleetP[[ff]]$retA[x,,1:nyears]
  }

  Effdist<-qfrac[x,]*Effind
  Efftot<-apply(Effdist,2,sum)

  MPAind<-TEG(c(nf,nyears,nareas))
  MPAtemp<-array(1/nf,dim(MPAc))
  MPAtemp[MPAind]=(MPAc[MPAind]*Effdist[MPAind[,1:2]])/Efftot[MPAind[,2]] # weighted by effort and fleet exposure
  MPAf<-apply(MPAtemp,2:3,sum)

  Vulnf<-Retc<-array(NA,c(maxage,nyears))
  Vind<-TEG(c(nf,maxage,nyears))
  Vtemp<-Vuln[,,1:nyears]
  Vtemp[Vind]<-(Vuln[Vind]*Effdist[Vind[,c(1,3)]])/Efftot[Vind[,3]]
  Vulnf[,1:nyears]<-apply(Vtemp,2:3,sum)
  #Vulnf[,(nyears+1):allyears]<-Vulnf[,nyears]  # Future vulnerability is the same
  Vulnf<-nlz(Vulnf,2,"max")       # normalize to max 1

  Rtemp<-retA[,,1:nyears]
  Rtemp[Vind]<-(retA[Vind]*Effdist[Vind[,c(1,3)]])/Efftot[Vind[,3]]
  Retc[,1:nyears]<-apply(Rtemp,2:3,sum)
  #Retc[,(nyears+1):allyears]<-Retc[,nyears]   # Future retention is the same

  Spat_targf<-sum(apply(Effdist,1,sum)*Spat_targc)/sum(Effdist) # Approximation according to historical F by fleet

  popdynCPP(nareas, maxage, Ncurr=N[x,p,,1,], nyears,
            M_age=StockP$M_ageArray[x,,], Asize_c=StockP$Asize[x,], MatAge=StockP$Mat_age[x,,],
            WtAge=StockP$Wt_age[x,,],Vuln=Vulnf, Retc=Retc, Prec=StockP$Perr_y[x,], movc=StockP$mov[x,,,], SRrelc=StockP$SRrel[x],
            Effind=Efftot,  Spat_targc=Spat_targf, hc=StockPars[[p]]$hs[x], R0c=R0a[x,],
            SSBpRc=StockP$SSBpR[x,], aRc=StockP$aR[x,], bRc=StockP$bR[x,], Qc=qs[x], Fapic=0, MPA=MPAf, maxF=maxF,
            control=1, SSB0c=StockP$SSB0[x])

}

# New function for reconstructing historical MICE model:

#' Reconstruct historical dynamics
#'
#' Function that reconstructs historical stock trends from fitted qs and all other parameters including MICE components
#'
#' @param x Integer, the simulation number
#' @param StockPars A list of sampled stock parameters, one list element per stock
#' @param FleetPars A hierarcical list of sampled fleet parameters, first list level is stock, second is fleet
#' @param np The number of stocks
#' @param nf The number of fleets
#' @param nareas The number of areas
#' @param maxage The maximum number of modelled ages
#' @param nyears The number of historical 'spool-up' years (from unfished to now)
#' @param N An array of stock numbers [nsim,np,maxage,nyears,nareas] - only the values from the first year are used
#' @param VF An array of vulnerability [nsim,np,nf,maxage,nyears+proyears]
#' @param FretA An array of retention [nsim,np,nf,maxage,nyears+proyears]
#' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
#' @param MPA An of spatial closures by year [np,nf,nyears+proyears,nareas]
#' @param CatchFrac A list of stock-specific fleet fractions of current catch list[[stock]][nsim, nf]
#' @param bounds Bounds for total q estimation
#' @param Rel A list of inter-stock relationships see slot Rel of MOM object class
#' @author T.Carruthers
#' @keywords internal
HistMICE<-function(x,StockPars, FleetPars, np,nf, nareas, maxage, nyears, N, VF, FretA, maxF=0.9, MPA,Rel,qs,qfrac){

  Nx<-array(N[x,,,,],dim(N)[2:5])
  VFx<-array(VF[x,,,,],dim(VF)[2:5])
  FretAx<-array(FretA[x,,,,],dim(VF)[2:5])
  #NIL(StockPars,"K")

  Kx<-matrix(unlist(lapply(StockPars,function(dat)dat['K'])),ncol=np)[x,]
  Linfx<-matrix(unlist(lapply(StockPars,function(dat)dat['Linf'])),ncol=np)[x,]
  t0x<-matrix(unlist(lapply(StockPars,function(dat)dat['t0'])),ncol=np)[x,]
  Mx<-matrix(unlist(lapply(StockPars,function(dat)dat['M'])),ncol=np)[x,]
  R0x<-matrix(unlist(lapply(StockPars,function(dat)dat['R0'])),ncol=np)[x,]

  hsx<-matrix(unlist(lapply(StockPars,function(dat)dat['hs'])),ncol=np)[x,]
  ax<-matrix(unlist(lapply(StockPars,function(dat)dat['a'])),ncol=np)[1,]
  bx<-matrix(unlist(lapply(StockPars,function(dat)dat['b'])),ncol=np)[1,]
  SRrelx<-matrix(unlist(lapply(StockPars,function(dat)dat['SRrel'])),ncol=np)[x,]

  distx<-SSBpRx<-R0ax<-aRx<-bRx<-array(NA,c(np,nareas))
  Perrx<-array(NA,c(np,nyears+maxage))
  movx<-array(NA,c(np,maxage,nareas,nareas))

  for(p in 1:np){
    distx[p,]<-StockPars[[p]]$R0a[x,]/sum(StockPars[[p]]$R0a[x,])
    Perrx[p,]<-StockPars[[p]]$Perr_y[x,1:(nyears+maxage)]
    movx[p,,,]<-StockPars[[p]]$mov[x,,,]
    SSBpRx[p,]<-StockPars[[p]]$SSBpR[x,]
    R0ax[p,]<-StockPars[[p]]$R0a[x,]
    aRx[p,]<-StockPars[[p]]$aR[x,]
    bRx[p,]<-StockPars[[p]]$bR[x,]
  }

  M_ageArrayx<-Mat_agex<-array(NA,c(np,maxage,nyears))
  Effind<-array(NA,c(np,nf,nyears))
  Spat_targ<-array(NA,c(np,nf))

  for(p in 1:np){
    Mat_agex[p,,]<-StockPars[[p]]$Mat_age[x,,1:nyears]
    M_ageArrayx[p,,]<-StockPars[[p]]$M_ageArray[x,,1:nyears]
    Effind[p,,]<-t(matrix(unlist(lapply(FleetPars[[p]],function(dat,x)dat['Find'][[1]][x,],x=x)),ncol=nf))
    Spat_targ[p,]<-unlist(lapply(FleetPars[[p]],function(dat,x)dat['Spat_targ'][[1]][x],x=x))
  }

  popdynMICE(qs=qs[x,],qfrac=qfrac[x,,],np,nf,nyears,nareas,maxage,Nx,VFx,FretAx,Effind,movx,Spat_targ,M_ageArrayx,Mat_agex,Kx,Linfx,t0x,Mx,R0x,R0ax,SSBpRx,hsx,aRx, bRx,ax,bx,Perrx,SRrelx,Rel)

}




#' Slot in list: get the slot values from a list of objects
#'
#' Create of vector of values that correspond with a slot in a list of objects
#'
#' @param listy A list of objects
#' @param sloty A character vector representing the slot name
#' @author T. Carruthers
#' @export
SIL<-function(listy,sloty) {
  if(class(listy[[1]])=="list"){
    out<-NULL
    for(i in 1:length(listy))out<-c(out,unlist(lapply(listy[[i]],function(x)slot(x,sloty))))
  }else{
    out<-unlist(lapply(listy,function(x)slot(x,sloty)))
  }
  out
}


#' Item in list: get the list values from a list of lists
#'
#' Create of vector of values that correspond with a slot in a list of objects
#'
#' @param listy A list of objects
#' @param namey A character vector representing the list item's name
#' @author T. Carruthers
#' @export
NIL<-function(listy,namey){

 unlist(lapply(listy,function(x)x[namey]))

}

#' Sum over list: get the list values from a list of lists
#'
#' Create of vector of values that correspond with a named position in a list of objects
#'
#' @param listy A list of objects
#' @param namey A character vector representing the list item's name
#' @author T. Carruthers
#' @export
SOL<-function(listy,namey){
  out<-listy[[1]][namey][[1]]
  for(i in 2:length(listy))out<-out+listy[[i]][namey][[1]]
  out
}

#' Toms expand grid
#'
#' Create an indexing grid from just a vector of maximum dimension sizes
#'
#' @param vec A vector of maximum array sizes
#' @author T. Carruthers
#' @export
TEG<-function(vec){ # make index for list calculation
  dim<-new('list')
  ndims<-length(vec)
  for(i in 1:ndims)dim[[i]]<-1:vec[i]
  as.matrix(expand.grid(dim))
}



#' Normalize
#'
#' Normalize an array over certain dimensions
#'
#' @param listy A list of objects
#' @param namey A character vector representing the list item's name
#' @author T. Carruthers
#' @export
nlz<-function(arr,dims=NULL,func="max"){
  arrdim<-dim(arr)
  ndim<-length(arrdim)
  if(is.null(dims))dims=1:(ndim-1)
  agg<-apply(arr,dims,func)
  out<-arr
  ind<-TEG(arrdim)
  out[ind]<-out[ind]/agg[ind[,dims]]
  out
}






