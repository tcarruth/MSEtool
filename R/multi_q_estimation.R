


getq_multi_MICE<-function(x,StockPars, FleetPars, nsim, np,nf, nareas,maxage, nyears, N, VF, FretA, maxF=MOM@maxF, MPA,CatchFrac, bounds,Rel){

  x<-1
  Nx<-N[x,,,,]
  VFx<-VF[x,,,,]
  FretAx<-FretA[x,,,,]
  NIL(StockPars,"K")

  Kx<-matrix(unlist(lapply(StockPars,function(dat)dat['K'])),ncol=np)[x,]
  Linfx<-matrix(unlist(lapply(StockPars,function(dat)dat['Linf'])),ncol=np)[x,]
  t0x<-matrix(unlist(lapply(StockPars,function(dat)dat['t0'])),ncol=np)[x,]
  Mx<-matrix(unlist(lapply(StockPars,function(dat)dat['M'])),ncol=np)[x,]
  R0x<-matrix(unlist(lapply(StockPars,function(dat)dat['R0'])),ncol=np)[x,]
  hsx<-matrix(unlist(lapply(StockPars,function(dat)dat['hs'])),ncol=np)[x,]
  ax<-matrix(unlist(lapply(StockPars,function(dat)dat['a'])),ncol=np)[x,]
  bx<-matrix(unlist(lapply(StockPars,function(dat)dat['b'])),ncol=np)[x,]
  SRrelx<-matrix(unlist(lapply(StockPars,function(dat)dat['SRrel'])),ncol=np)[x,]

  distx<-array(NA,c(np,nareas))
  Perrx<-array(NA,c(np,nyears+maxage))
  movx<-array(NA,c(np,maxage,nareas,nareas))
  for(p in 1:np){
    distx[p,]<-StockPars[[p]]$R0a[x,]/sum(StockPars[[p]]$R0a[x,])
    Perrx[p,]<-StockPars[[p]]$Perr_y[x,1:(nyears+maxage)]
    movx[p,,,]<-StockPars[[p]]$mov[x,,,]
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

  CF<-t(matrix(unlist(lapply(CatchFrac,function(dat)dat[x,])),nrow=nf))
  Fdist<-CF/Effind[,,nyears] # Catch divided by effort (q proxy)
  Fdist<-Fdist/apply(Fdist,1,sum)    # q ratio proxy (real space)
  par<-c(rep(-5,np),log(Fdist[,2:nf]/(1-Fdist[,2:nf]))) # low initial F followed by logit guess at fraction based on Fdist according to catch fraction in recent year

  q_estMICe(par,np,nf,nyears,SSBx,Nx,VFx,FretAx,Effind,distx,movx,Spat_targ,M_ageArrayx,Mat_agex,Perrx,Kx,Linfx,t0x,Mx,R0x,hsx,ax,bx,SRrelx,Rel)

  optim(par,qestMICE)





  }



q_estMICe<-function(par,np,nf,nyears,Nx,VFx,FretAx,Effind,distx,movx,Spat_targ,M_ageArrayx,Mat_agex,Kx,Linfx,t0x,Mx,R0x,hsx,ax,bx,SRrelx,Rel){

  qs<-exp(par[1:np])
  qlogit<-array(0,c(np,nf))
  qlogit[,2:nf]<-par[(np+1):length(par)]
  qfrac<-exp(qlogit)/apply(exp(qlogit),1,sum)
  # Effind [np,nf,nyears]
  # Effdist<-qfrac*Effind
  # Efftot<-apply(Effdist,2,sum)

  # make all arrays for reallocation here (each year) and make this 'popdyn_MICE'

  popdynMICE<-function(qs,qfrac,np,nf,nyears,Nx,VFx,FretAx,Effind,distx,movx,Spat_targ,M_ageArrayx,Mat_agex,Kx,Linfx,t0x,Mx,R0x,hsx,ax,bx,SRrelx,Rel){

    Bx<-SSBx<-array(NA,dim(Nx))
    Wt_agey<-array(NA,c(np,maxage,nyears))
    Ky<-Linfy<-t0y<-My<-R0y<-hsy<-ay <-by<-array(NA,c(np,nyears))
    Ky[,1]<-Kx; Linfy[,1]<-Linfx; t0y[,1]<-t0x; My[,1]<-Mx; R0y[,1]<-R0x; hsy[,1]<-hsx; ay[,1]<-ax; by[,1]<-bx

    Len_age<-matrix(Linfx*(1-exp(-(rep(1:maxage,each=np)-t0x)*(Kx))),nrow=np)
    Wt_agey[,,1]<-ax*Len_age^bx

    for(y in 2:nyears){

      #y<-2; Fcur<-apply(Effind[,,y-1]*qs*qfrac,1,sum) ; Ncur<-Nx[,,y-1,];  Vcur<-VFx[,,,y-1];  Retcur<-FretAx[,,,y-1];  M_agecur=M_ageArrayx[,,y-1];Mat_agecur=Mat_agex[,,y-1];    PerrYrp=Perrx[,y+maxage-2]

      Fcur<-apply(Effind[,,y-1]*qs*qfrac,1,sum)

      out<-popdynOneMICE(np,nf,nareas, maxage, Ncur=Nx[,,y-1,],
                            Vcur=VFx[,,,y-1], Retcur=FretAx[,,,y-1], Fcur=Fcur, PerrYrp=Perrx[,y+maxage-2], hsx=hsy[,y-1],
                            distx=distx, movx=movx, SRrelx=SRrelx, M_agecur=M_ageArrayx[,,y-1], Mat_agecur=Mat_agex[,,y-1],
                            Kx=Ky[,y-1], Linfx=Linfy[,y-1], t0x=t0y[,y-1], Mx=My[,y-1], R0x=R0y[,y-1], ax=ay[,y-1],
                            bx=by[,y-1],Rel=Rel)

      Nx[,,y,]<-out$Nnext
      Wt_agey[,,y]<-out$Wt_age
      M_ageArrayx[,,y]<-out$M_agecurx
      Ky[,y]<-out$Kx; Linfy[,y]<-out$Linfx; t0y[,y]<-out$t0x; My[,y]<-out$Mx; R0y[,y]<-out$R0x; hsy[,y]<-out$hsx; ay[,y]<-out$ax; by[,y]<-out$bx

    }

    Nind<-TEG(dim(Nx))
    Bx[Nind]<-Nx[Nind]*Wt_agey[Nind[,1:3]]
    SSBx[Nind]<-Bx[Nind]*Mat_agex[Nind[,1:3]]
    # matplot(t(apply(SSBx,c(1,3),sum)))
  }


}


popdynOneMICE<-function(np,nf,nareas, maxage, Ncur, Vcur, Retcur, Fcur, PerrYrp, hsx,R0a, distx, movx, SRrelx,M_agecur,Mat_agecur,
                        Kx,Linfx,t0x,Mx,R0x,ax,bx,herm=NULL, Rel){

  # FMarray.subcube(0,0, A, maxage-1, 0, A) =  (Effind(0) * Qc * fishdist(A) * Vuln.col(0))/Asize_c(A);
  # FMretarray.subcube(0,0, A, maxage-1, 0, A) =  (Effind(0) * Qc * fishdist(A) * Retc.col(0))/Asize_c(A);

  # Initial Bcur calc (before any weight at age recalculation change)
  # Bcalc ---------------------------------------------------------------------------
  Bcur<-SSBcur<-Ncur
  Nind<-TEG(dim(Ncur)) # p, age, area
  Len_age<-matrix(Linfx*(1-exp(-(rep(1:maxage,each=np)-t0x)*(Kx))),nrow=np)
  Wt_age<-ax*Len_age^bx
  Bcur[Nind]<-Ncur[Nind]*Wt_age[Nind[,1:2]]

  # old surv
  surv <- cbind(rep(1,np),t(exp(-apply(M_agecur, 1, cumsum)))[, 1:(maxage-1)])  # Survival array
  oldM<-apply(surv*M_agecur*Mat_agecur,1,sum)/apply(surv*Mat_agecur,1,sum)

  Responses<-ResFromRel(Rel,Bcur,SSBcur,Ncur,seed=1) # ------------------------------------
  for(rr in 1:nrow(Responses))  eval(parse(text=paste0(Responses[rr,4],"[",Responses[rr,3],"]<-",as.numeric(Responses[rr,1]))))
  # Parameters that could have changed R0, M, K, Linf, t0, a, b, hs, bR, aR
  # Recalculate growth
  Len_age<-matrix(Linfx*(1-exp(-(rep(1:maxage,each=np)-t0x)*(Kx))),nrow=np)
  Wt_age<-ax*Len_age^bx

  # Recalc SSBpR, SSB0
  M_agecurx<-M_agecur*Mx/oldM  # updated M
  surv <- cbind(rep(1,np),t(exp(-apply(M_agecurx, 1, cumsum)))[, 1:(maxage-1)])  # Survival array
  SSB0x<-apply(R0x*surv*Mat_agecur*Wt_age,1,sum)
  SSBpRx<-SSB0x/R0x
  SSBpRax<-SSBpRx*dist
  SSB0ax<-dist*SSB0x
  R0ax<-dist*R0x

  # R0recalc thus aR bR recalc ---------------
  bRx <- matrix(log(5 * hsx)/(0.8 * SSB0ax), nrow=np)  # Ricker SR params
  aRx <- matrix(exp(bRx * SSB0ax)/SSBpRx, nrow=np)  # Ricker SR params

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

  FMx[VBind]<-Fdist[VBind]*Fcur[VBind[,1]]*Vcur[VBind[,1:3]]
  FMretx[VBind]<-Fdist[VBind]*Fcur[VBind[,1]]*Retcur[VBind[,1:3]]
  Zx[VBind]<-FMx[VBind]+M_agecur[VBind[,c(1,3)]]
  Zcur<-apply(Zx,c(1,3,4),sum)

  Nnext<-array(NA,c(np,maxage,nareas))

  for(p in 1:np){
    #NextYrN<-popdynOneTScpp(nareas, maxage, SSBcurr=colSums(SSB[x,,nyears, ]), Ncurr=N[x,,nyears,],
    #                        Zcurr=Zcur[p,,], PerrYr=PerrYrp[x, nyears+maxage-1], hs=hs[x],
    #                        R0a=R0a[x,], SSBpR=SSBpRax[p,], aR=aRx[p,], bR=bRx[p,],
    #                        mov=mov[x,,,], SRrel=SRrelx[p])

    NextYrN<-popdynOneTScpp(nareas, maxage, SSBcurr=SSBcur[p,,], Ncurr=Ncur[p,,],
                            Zcurr=Zcur[p,,], PerrYr=PerrYrp[p], hs=hsx[p],
                            R0a=R0ax[p,], SSBpR=SSBpRax[p,], aR=aR[p,], bR=bR[p,],
                            mov=movx[p,,,], SRrel=SRrelx[p])

    Nnext[p,,]<-NextYrN

  }

  # returns new N and any updated parameters:
  list(Nnext=Nnext,M_agecurx=M_agecurx,R0x=R0x,R0ax=R0ax,hsx=hsx,
       aRx=aRx,bRx=bRx,Linfx=Linfx,Kx=Kx,t0x=t0x,Mx=Mx,ax=ax,bx=bx,
       Len_age=Len_age,Wt_age=Wt_age,surv=surv)

}



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




#' optimize for catchability (q)
#'
#' Function optimizes catchability (q, where F=qE) required to get to user-specified stock
#' depletion
#'
#' @param x Integer, the simulation number
#' @param D A numeric vector nsim long of sampled depletion
#' @param SSB0 A numeric vector nsim long of total unfished spawning biomass
#' @param nareas The number of spatial areas
#' @param maxage The maximum age
#' @param N Array of the numbers-at-age in population. Dimensions are nsim, maxage, nyears, nareas.
#' Only values from the first year (i.e `N[,,1,]`) are used, which is the current N-at-age.
#' @param pyears The number of years to project forward. Equal to 'nyears' for optimizing for q.
#' @param M_ageArray An array (dimensions nsim, maxage, nyears+proyears) with the natural mortality-at-age and year
#' @param Mat_age An array (dimensions nsim, maxage, proyears+nyears) with the proportion mature for each age-class
#' @param Asize A matrix (dimensions nsim, nareas) with size of each area
#' @param Wt_age An array (dimensions nsim, maxage, nyears+proyears) with the weight-at-age and year
#' @param FleetP A list of Fleets containing entries for V, retA, Find and Spat_targ
#' @param mov An array (dimensions nsim, nareas, nareas) with the movement matrix
#' @param SRrel A numeric vector nsim long specifying the recruitment curve to use
#' @param hs A numeric vector nsim long with the steepness values for each simulation
#' @param R0a A matrix (dimensions nsim, nareas) with the unfished recruitment by area
#' @param SSBpR A matrix (dimensions nsim, nareas) with the unfished spawning-per-recruit by area
#' @param aR A numeric vector nareas long with the Ricker SRR a values
#' @param bR A numeric vector nareas long with the Ricker SRR b values
#' @param bounds A numeric vector of length 2 with bounds for the optimizer
#' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
#' @param MPA A matrix of spatial closures by year
#' @param useCPP logical - use the CPP code? For testing purposes only
#' @author T.Carruthers and A. Hordyk
#' @keywords internal
getq_multi <- function(x, D, SSB0, nareas, maxage, Np, pyears, M_ageArray, Mat_age, Asize, Wt_age,
                  FleetP, Perr, mov, SRrel, hs, R0a, SSBpR, aR, bR,
                  bounds = c(1e-05, 15),maxF, MPAc, CFp, useCPP=TRUE) {

  # Prep multi-fleet arrays
  nf<-length(FleetP)
  ns<-dim(FleetP[[1]]$V)[1]
  ay<-dim(FleetP[[1]]$V)[3]
  Vuln=retA=array(NA,c(nf,maxage,ay))
  ny<-pyears

  Effind<-t(matrix(unlist(lapply(FleetP,function(dat,x)dat['Find'][[1]][x,],x=x)),ncol=nf))
  Spat_targc<-unlist(lapply(FleetP,function(dat,x)dat['Spat_targ'][[1]][x],x=x))

  for(ff in 1:nf){ # kind of lazy but not THAT slow

    Vuln[ff,,]<-FleetP[[ff]]$V[x,,]
    retA[ff,,]<-FleetP[[ff]]$retA[x,,]

  }

  Fdist<-CFp[x,]/Effind[,ny] # Catch divided by effort (q proxy)
  Fdist<-Fdist/sum(Fdist)    # q ratio proxy (real space)
  par<-c(-5,log(Fdist[2:nf]/(1-Fdist[2:nf]))) # low initial F followed by logit guess at fraction based on Fdist according to catch fraction in recent year

  opt <- optim(par,optQ_multi, method="L-BFGS-B",lower=c(log(bounds[1]),rep(-5,nf-1)),upper=c(log(bounds[2]),rep(5,nf-1)),depc=D[x], SSB0c=SSB0[x], nareas=nareas, maxage=maxage, Ncurr=Np[x,,1,],
                  pyears=pyears, M_age=M_ageArray[x,,], MatAge=Mat_age[x,,], Asize_c=Asize[x,], WtAge=Wt_age[x,,],
                  Vuln=Vuln, Retc=retA, Prec=Perr[x,], movc=mov[x,,,], SRrelc=SRrel[x],
                  Effind=Effind,  Spat_targc=Spat_targc, hc=hs[x], R0c=R0a[x,],
                  SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], maxF=maxF, MPAc=MPAc, CFc=CFp[x,],useCPP=useCPP,mode='opt')

  out<-optQ_multi(opt$par,depc=D[x], SSB0c=SSB0[x], nareas=nareas, maxage=maxage, Ncurr=Np[x,,1,],
                  pyears=pyears, M_age=M_ageArray[x,,], MatAge=Mat_age[x,,], Asize_c=Asize[x,], WtAge=Wt_age[x,,],
                  Vuln=Vuln, Retc=retA, Prec=Perr[x,], movc=mov[x,,,], SRrelc=SRrel[x],
                  Effind=Effind,  Spat_targc=Spat_targc, hc=hs[x], R0c=R0a[x,], SSBpRc=SSBpR[x,],
                  aRc=aR[x,], bRc=bR[x,], maxF=maxF, MPAc=MPAc, CFc=CFp[x,],useCPP=useCPP,mode='out')

  return(out)

}


#' Optimize q for a single simulation
#'
#' @param logQ log q
#' @param depc Depletion value
#' @param SSB0c Unfished spawning biomass
#' @param nareas Number of areas
#' @param maxage Maximum age
#' @param Ncurr Current N-at-age
#' @param pyears Number of years to project population dynamics
#' @param M_age M-at-age
#' @param Asize_c Numeric vector (length nareas) with size of each area
#' @param MatAge Maturity-at-age
#' @param WtAge Weight-at-age
#' @param Vuln Vulnerability-at-age
#' @param Retc Retention-at-age
#' @param Prec Recruitment error by year
#' @param movc movement matrix
#' @param SRrelc SR parameter
#' @param Effind Historical fishing effort
#' @param Spat_targc Spatial targetting
#' @param hc Steepness
#' @param R0c Unfished recruitment by area
#' @param SSBpRc Unfished spawning biomass per recruit by area
#' @param aRc Ricker aR
#' @param bRc Ricker bR
#' @param maxF maximum F
#' @param MPA A matrix of spatial closures by year
#' @param useCPP Logical. Use the CPP code?
#' @param model Character 'opt' is optimization and returns an objective funciton to be minimized, anything else returns fitted values
#' @author A. Hordyk
#' @keywords internal
optQ_multi <- function(par, depc, SSB0c, nareas, maxage, Ncurr, pyears, M_age,
                 MatAge, Asize_c, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc,
                 R0c, SSBpRc, aRc, bRc, maxF, MPAc, CFc, useCPP, mode = 'opt') {


    nf<-dim(Effind)[1]
    ny<-dim(Effind)[2]
    ay<-dim(Vuln)[3]

    qtot<-exp(par[1])
    qlogit<-c(0,par[2:nf])
    qfrac<-exp(qlogit)/(sum(exp(qlogit)))
    Effdist<-qfrac*Effind
    Efftot<-apply(Effdist,2,sum)

    MPAind<-TEG(c(nf,ny,nareas))
    MPAtemp<-array(1/nf,dim(MPAc))
    MPAtemp[MPAind]=(MPAc[MPAind]*Effdist[MPAind[,1:2]])/Efftot[MPAind[,2]] # weighted by effort and fleet exposure
    MPAf<-apply(MPAtemp,2:3,sum)

    Vulnf<-Retf<-array(NA,c(maxage,ay))
    Vind<-TEG(c(nf,maxage,ny))
    Vtemp<-Vuln[,,1:ny]
    Vtemp[Vind]<-(Vuln[Vind]*Effdist[Vind[,c(1,3)]])/Efftot[Vind[,3]]
    Vulnf[,1:ny]<-apply(Vtemp,2:3,sum)
    Vulnf[,(ny+1):ay]<-Vulnf[,ny]  # Future vulnerability is the same
    Vulnf<-nlz(Vulnf,2,"max")       # normalize to max 1

    Rtemp<-Retc[,,1:ny]
    Rtemp[Vind]<-(Retc[Vind]*Effdist[Vind[,c(1,3)]])/Efftot[Vind[,3]]
    Retf[,1:ny]<-apply(Rtemp,2:3,sum)
    Retf[,(ny+1):ay]<-Retf[,ny]   # Future retention is the same

    Spat_targf<-sum(apply(Effdist,1,sum)*Spat_targc)/sum(Effdist) # Approximation according to historical F by fleet

    simpop <- popdynCPP(nareas, maxage, Ncurr, pyears, M_age, Asize_c,
                        MatAge, WtAge, Vuln=Vulnf, Retc=Retf, Prec, movc, SRrelc, Effind=Efftot, Spat_targc=Spat_targf, hc,
                        R0c=R0c, SSBpRc=SSBpRc, aRc=aRc, bRc=bRc, Qc=qtot, Fapic=0,
                        maxF=maxF, MPA=MPAf, control=1,  SSB0c=SSB0c)



    # 1:Narray 2:Barray 3:SSNarra 4:SBarray 5:VBarray 6:FMarray 7:FMretarray 8:Zarray
    ssb <- sum(simpop[[4]][,pyears,])
    BB<- apply(simpop[[2]][,pyears,],1,sum)  # age area
    Cf<-array(rep(BB,each=nf)*Vuln[,,ny]*Effdist[,ny],c(nf,maxage))
    Cpred<-apply(Cf,1,sum)/sum(Cf)

    depOBJ<-50*(log(depc) - log(ssb/SSB0c))^2
    cOBJ<-sum(log(Cpred[2:nf]/CFc[2:nf])^2)

    if(mode=='opt'){
      return(depOBJ+cOBJ)
    }else{
      return(list(qtot=qtot,qfrac=qfrac,CFc=CFc,Cpred=Cpred,depc=depc,deppred=ssb/SSB0c,Vulnf=Vulnf,Retf=Retf,MPAf=MPAf))
    }

}



