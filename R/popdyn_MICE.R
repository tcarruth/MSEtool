#' Population dynamics for a MICE model (multiyear)
#'
#' Calls popdynOneMICE iteratively to reconstruct a time series given MICE model inputs
#'
#' @param qs Total catchability
#' @param qfrac Vector [fleet], the fraction of total qs by fleet
#' @param np Integer, the number of stocks
#' @param nf Integer, number of fleets
#' @param nyears Integer, number of historical years (unfished til today)
#' @param nareas Integer, the number of spatial areas
#' @param maxage Integer, maximum modelled age
#' @param Nx Array [stock, age, year, area] of stock numbers
#' @param VFx Array [fleet, age, year, area] of the vulnerability curve
#' @param FretAx Array [fleet, age, year, area] of the retention curve
#' @param Effind Array [fleet, year] of effort
#' @param movx Array [stock,age,area,area] of movement transitions
#' @param Spat_targ Matrix [stock, fleet] of spatial targetting parameter (0 evenly spatial distributed, 1 proportional to vulnerable biomass)
#' @param M_ageArrayx Array [stock, age,year] of Natural mortality rate at age
#' @param Mat_agex Array [stock, age, year] of maturity (spawning fraction) age age
#' @param Kx Vector [stock] of von B growth parameter K
#' @param Linf Vector [stock] of von B asymptotic length parameter Linf
#' @param t0 Vector [stock] of von B theoretical age at zero length (t0)
#' @param Mx Vector [stock] mature natural mortality rate
#' @param R0x Vector [stock] unfished recruitment
#' @param R0ax Matrix [stock, area] unfished recruitment by area
#' @param SSBpRx Matrix [stock, area] spawning biomass per recruit by area
#' @param hsx Vector [stock] steepness of the stock recruitment curve
#' @param aRx Vector [stock] stock recruitment parameter alpha (for Ricker curve)
#' @param bRx Vector [stock] stock recruitment parameter beta (for Ricker curve)
#' @param ax Vector [stock] weight-length parameter a W=aL^b
#' @param bx Vector [stock] weight-length parameter b W=aL^b
#' @param Perrx Matrix [stock, year] process error - the lognormal factor for recruitment strength
#' @param SRrelx Integer vector [stock] the form of the stock recruitment relationship (1 = Beverton-Holt, 2= Ricker)
#' @param Rel A list of inter-stock relationships see slot Rel of MOM object class
#' @author T.Carruthers
#' @keywords internal
#' @export
popdynMICE<-function(qs,qfrac,np,nf,nyears,nareas,maxage,Nx,VFx,FretAx,Effind,movx,Spat_targ,M_ageArrayx,Mat_agex,Kx,Linfx,t0x,Mx,R0x,R0ax,SSBpRx,hsx,aRx, bRx,ax,bx,Perrx,SRrelx,Rel){

  Bx<-SSBx<-array(NA,dim(Nx))
  Fy<-array(NA,c(np,nf,nyears))
  FMy<-array(NA,c(np,nf,maxage,nyears,nareas))
  Wt_agey<-array(NA,c(np,maxage,nyears))
  Ky<-Linfy<-t0y<-My<-hsy<-ay <-by<-array(NA,c(np,nyears))
  Ky[,1]<-Kx; Linfy[,1]<-Linfx; t0y[,1]<-t0x; My[,1]<-Mx; hsy[,1]<-hsx; ay[,1]<-ax; by[,1]<-bx

  Len_age<-matrix(Linfx*(1-exp(-(rep(1:maxage,each=np)-t0x)*(Kx))),nrow=np)
  Wt_agey[,,1]<-ax*Len_age^bx

  for(y in 2:nyears){

    # y<-y+1
    Fy[,,y-1]<-Effind[,,y-1]*qs*qfrac
    # y<-2; M_agecur=M_ageArrayx[,,y-1];Mat_agecur=Mat_agex[,,y-1];    PerrYrp=Perrx[,y+maxage-2]
    Vcur=array(VFx[,,,y-1],dim(VFx)[1:3])
    Retcur=array(FretAx[,,,y-1],dim(FretAx)[1:3])
    Fcur=array(Fy[,,y-1],dim(Fy)[1:2])
    Ncur=array(Nx[,,y-1,],dim(Nx)[c(1:2,4)])
    M_agecur=array(M_ageArrayx[,,y-1],dim(M_ageArrayx)[1:2])
    Mat_agecur<-array(Mat_agex[,,y-1],dim(Mat_agex)[1:2])

    out<-popdynOneMICE(np,nf,nareas, maxage, Ncur=Ncur, Vcur=Vcur, Retcur=Retcur, Fcur=Fcur, PerrYrp=Perrx[,y+maxage-2], hsx=hsy[,y-1], aRx=aRx, bRx=bRx,
                       movx=movx, Spat_targ=Spat_targ, SRrelx=SRrelx, M_agecur=M_agecur, Mat_agecur=Mat_agecur,
                       Kx=Ky[,y-1], Linfx=Linfy[,y-1], t0x=t0y[,y-1], Mx=My[,y-1], R0x=R0x,R0ax=R0ax,SSBpRx=SSBpRx,ax=ay[,y-1],
                       bx=by[,y-1],Rel=Rel)

    Nx[,,y,]<-out$Nnext
    Wt_agey[,,y]<-out$Wt_age
    M_ageArrayx[,,y]<-out$M_agecurx
    Ky[,y]<-out$Kx; Linfy[,y]<-out$Linfx; t0y[,y]<-out$t0x; My[,y]<-out$Mx; hsy[,y]<-out$hsx; ay[,y]<-out$ax; by[,y]<-out$bx
    FMy[,,,y,]<-out$FMx

  }

  Nind<-TEG(dim(Nx))
  Bx[Nind]<-Nx[Nind]*Wt_agey[Nind[,1:3]]
  SSBx[Nind]<-Bx[Nind]*Mat_agex[Nind[,1:3]]
  # matplot(t(apply(SSBx,c(1,3),sum)))
  # matplot(t(apply(Nx,c(1,3),sum)))

  list(Nx=Nx,Bx=Bx,SSBx=SSBx,Fy=Fy,FMy=FMy,Ky=Ky,Linfy=Linfy,t0y=t0y,My=My,hsy=hsy,ay=ay,by=by)
}



#' Population dynamics for a MICE model (single year)
#'
#' Completes a single iteration of recruitment, mortality, fishing and movement given MICE model inputs
#'
#' @param np Integer, the number of stocks
#' @param nf Integer, number of fleets
#' @param nyears Integer, number of historical years (unfished til today)
#' @param maxage Integer, maximum modelled age
#' @param Ncur Array [stock, age, area] of stock numbers
#' @param Vcur Array [fleet, age, area] of the vulnerability curve
#' @param Retcur Array [fleet, age, area] of the retention curve
#' @param Fcur Array [stock, fleet] fishing mortality rate
#' @param PerrYrp Vector [stock] process error - the lognormal factor for recruitment strength
#' @param hsx Vector [stock] steepness of the stock recruitment curve
#' @param aRx Vector [stock] stock recruitment parameter alpha (for Ricker curve)
#' @param bRx Vector [stock] stock recruitment parameter beta (for Ricker curve)
#' @param movx Array [stock,age,area,area] of movement transitions
#' @param Spat_targ Matrix [stock, fleet] of spatial targetting parameter (0 evenly spatial distributed, 1 proportional to vulnerable biomass)
#' @param SRrelx Integer vector [stock] the form of the stock recruitment relationship (1 = Beverton-Holt, 2= Ricker)
#' @param M_agecur Matrix [stock, age] of Natural mortality rate at age
#' @param Mat_agecur Matrix [stock, age] of maturity (spawning fraction) age age
#' @param Kx Vector [stock] of von B growth parameter K
#' @param Linf Vector [stock] of von B asymptotic length parameter Linf
#' @param t0 Vector [stock] of von B theoretical age at zero length (t0)
#' @param Mx Vector [stock] mature natural mortality rate
#' @param R0x Vector [stock] unfished recruitment
#' @param R0ax Matrix [stock, area] unfished recruitment by area
#' @param SSBpRx Matrix [stock, area] spawning biomass per recruit by area
#' @param ax Vector [stock] weight-length parameter a W=aL^b
#' @param bx Vector [stock] weight-length parameter b W=aL^b
#' @param herm Character are there two stocks that are hermaphroditic pg = protogynous (F-M), pn = protandrous (M-F)?
#' @param sexspecific Logical are the two stocks two sexes (female is stock 1, male is stock 2)
#' @param Rel A list of inter-stock relationships see slot Rel of MOM object class
#' @author T.Carruthers
#' @keywords internal
#' @export
popdynOneMICE<-function(np,nf,nareas, maxage, Ncur, Vcur, Retcur, Fcur, PerrYrp, hsx, aRx, bRx, movx,Spat_targ,
                        SRrelx,M_agecur,Mat_agecur,
                        Kx,Linfx,t0x,Mx,R0x,R0ax,SSBpRx,ax,bx,herm=NULL, sexspecific=NULL, Rel){

  # FMarray.subcube(0,0, A, maxage-1, 0, A) =  (Effind(0) * Qc * fishdist(A) * Vuln.col(0))/Asize_c(A);
  # FMretarray.subcube(0,0, A, maxage-1, 0, A) =  (Effind(0) * Qc * fishdist(A) * Retc.col(0))/Asize_c(A);

  # Initial Bcur calc (before any weight at age recalculation change)
  # Bcalc ---------------------------------------------------------------------------
  Bcur<-SSBcur<-array(NA,dim(Ncur))
  Nind<-TEG(dim(Ncur)) # p, age, area
  Len_age<-matrix(Linfx*(1-exp(-(rep(1:maxage,each=np)-t0x)*(Kx))),nrow=np)
  Wt_age<-ax*Len_age^bx
  Bcur[Nind]<-Ncur[Nind]*Wt_age[Nind[,1:2]]

  # old surv
  # surv <- matrix(1, nsim, maxage)
  #surv[, 2:maxage] <- t(exp(-apply(StockPars[[1]]$M_ageArray[,,1], 1, cumsum)))[, 1:(maxage-1)]  # Survival array
  surv <- array(c(rep(1,np),t(exp(-apply(M_agecur, 1, cumsum)))[, 1:(maxage-1)]),c(np,maxage))  # Survival array
  oldM<-apply(surv*M_agecur*Mat_agecur,1,sum)/apply(surv*Mat_agecur,1,sum)

  if(np>1 & length(Rel)>0){
    Responses<-ResFromRel(Rel,Bcur,SSBcur,Ncur,seed=1) # ------------------------------------
    for(rr in 1:nrow(Responses))  eval(parse(text=paste0(Responses[rr,4],"[",Responses[rr,3],"]<-",as.numeric(Responses[rr,1]))))
  }
  # Parameters that could have changed: M, K, Linf, t0, a, b, hs
  # Recalculate growth
  Len_age<-matrix(Linfx*(1-exp(-(rep(1:maxage,each=np)-t0x)*(Kx))),nrow=np)
  Wt_age<-ax*Len_age^bx

  # Recalc SSBpR, SSB0
  M_agecurx<-M_agecur*Mx/oldM  # updated M

  # This is redundant code for updating parameters when R0 changes
  #surv <- cbind(rep(1,np),t(exp(-apply(M_agecurx, 1, cumsum)))[, 1:(maxage-1)])  # Survival array
  #SSB0x<-apply(R0x*surv*Mat_agecur*Wt_age,1,sum)
  #SSBpRx<-SSB0x/R0x
  #SSBpRax<-SSBpRx*distx
  #SSB0ax<-distx*SSB0x
  #R0ax<-distx*R0x
  #R0recalc thus aR bR recalc ---------------
  #bRx <- matrix(log(5 * hsx)/(0.8 * SSB0ax), nrow=np)  # Ricker SR params
  #aRx <- matrix(exp(bRx * SSB0ax)/SSBpRx, nrow=np)  # Ricker SR params

  # Bcalc ---------------------------------------------------------------------------

  Bcur[Nind]<-Ncur[Nind]*Wt_age[Nind[,1:2]]
  SSBcur[Nind]<-Bcur[Nind]*Mat_agecur[Nind[,1:2]]

  # Vulnerable biomass calculation --------------------------------------------------
  VBx<-Fdist<-FMx<-FMretx<-Zx<-array(NA,c(np,nf,maxage,nareas))
  VBind<-TEG(dim(VBx))
  VBx[VBind]<-Vcur[VBind[,1:3]]*Bcur[VBind[,c(1,3:4)]]
  Fdist[VBind]<-VBx[VBind]^Spat_targ[VBind[,1:2]]
  VBagg<-apply(Fdist,1:3,sum)
  Fdist[VBind]<-Fdist[VBind]/VBagg[VBind[,1:3]]

  FMx[VBind]<-Fdist[VBind]*Fcur[VBind[,1:2]]*Vcur[VBind[,1:3]]
  FMretx[VBind]<-Fdist[VBind]*Fcur[VBind[,1:2]]*Retcur[VBind[,1:3]]
  Ft<-apply(FMx,c(1,3,4),sum)#FMx[VBind]+M_agecur[VBind[,c(1,3)]]
  Zcur<-Ft+array(rep(M_agecur[VBind[,c(1,3)]],nareas),c(np,maxage,nareas))

  Nnext<-array(NA,c(np,maxage,nareas))

  for(p in 1:np){
    #NextYrN<-popdynOneTScpp(nareas, maxage, SSBcurr=colSums(SSB[x,,nyears, ]), Ncurr=N[x,,nyears,],
    #                        Zcurr=Zcur[p,,], PerrYr=PerrYrp[x, nyears+maxage-1], hs=hs[x],
    #                        R0a=R0a[x,], SSBpR=SSBpRax[p,], aR=aRx[p,], bR=bRx[p,],
    #                        mov=mov[x,,,], SRrel=SRrelx[p])

    NextYrN<-popdynOneTScpp(nareas, maxage, SSBcurr=colSums(SSBcur[p,,]), Ncurr=Ncur[p,,],
                            Zcurr=Zcur[p,,], PerrYr=PerrYrp[p], hs=hsx[p],
                            R0a=R0ax[p,], SSBpR=SSBpRx[p,], aR=aRx[p,], bR=bRx[p,],
                            mov=movx[p,,,], SRrel=SRrelx[p])

    Nnext[p,,]<-NextYrN

  }

  # returns new N and any updated parameters:
  list(Nnext=Nnext,M_agecurx=M_agecurx,R0x=R0x,R0ax=R0ax,hsx=hsx,
       aRx=aRx,bRx=bRx,Linfx=Linfx,Kx=Kx,t0x=t0x,Mx=Mx,ax=ax,bx=bx,
       Len_age=Len_age,Wt_age=Wt_age,surv=surv,FMx=FMx)

}


#' Returns Results of a set of MICE relationships
#'
#' Predicts stock-specific parameters from another stocks biomass, spawning biomass or numbers
#'
#' @param Rel A list of inter-stock relationships see slot Rel of MOM object class
#' @param Bcur An array of current stock biomass [stock, age, area]
#' @param SSBcur An array of current spawning stock biomass [stock, age, area]
#' @param Ncur An array of current stock numbers [stock, age, area]
#' @author T.Carruthers
#' @keywords internal
#' @export
ResFromRel<-function(Rel,Bcur,SSBcur,Ncur,seed){

  IVnams<-c("B","SSB","N")
  IVcode<-c("Bcur","SSBcur","Ncur")

  DVnam<-c("M","a", "b", "R0", "hs", "K", "Linf", "t0")
  modnam<-c("Mx","ax","bx","R0x","hsx","Kx","Linfx","t0x")

  nRel<-length(Rel)
  out<-array(NA,c(nRel,4))

  for(r in 1:nRel){

    fnams<-names(attr(Rel[[r]]$terms,"dataClasses"))
    DV<-fnams[1]
    Dp<-unlist(strsplit(DV,"_"))[2]
    Dnam<-unlist(strsplit(DV,"_"))[1]
    IV<-fnams[2:length(fnams)]
    nIV<-length(IV)
    IVs<-matrix(unlist(strsplit(IV,"_")),ncol=nIV)
    newdat<-NULL
    for(iv in 1:nIV){
      p<-as.numeric(IVs[2,iv])
      if(IVs[1,iv]=="B")newdat<-rbind(newdat,sum(Bcur[p,,]))
      if(IVs[1,iv]=="SSB")newdat<-rbind(newdat,sum(SSBcur[p,,]))
      if(IVs[1,iv]=="N")newdat<-rbind(newdat,sum(Ncur[p,,]))
    }
    newdat<-as.data.frame(newdat)
    names(newdat)<-IV
    ys<-predict(Rel[[r]],newdat=newdat)
    templm<-Rel[[r]]
    templm$fitted.values <- ys
    out[r,1]<-unlist(simulate(templm, nsim=1, seed=seed))
    out[r,2]<-DV
    out[r,3]<-Dp
    out[r,4]<-modnam[match(Dnam,DVnam)]

  }

  out

}


