# multiMSE funcs


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

 out<-unlist(lapply(listy,function(x)x[namey]))

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






